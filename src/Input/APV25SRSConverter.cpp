#include "Input/APV25SRSConverter.h"
#include "Terminal.h"

#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile2D.h>
#include <TTree.h>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <tuple>
#include <vector>

namespace fs = std::filesystem;
using json = nlohmann::json;

namespace {
struct APV25SRSChannelInfo {
   int boardID{0};
   int chipID{0};
   int channelID{0};
   int detectorID{0};
   int planeType{0};
   int id0{0};
   int id1{-1};
   double noiseSigma{-1.0};
};

struct GeometryChannelMapping {
   int detectorID{0};
   int chipID{0};
   int planeID{0};
   std::vector<int> id0;
   std::optional<std::vector<int>> id1;
};

using ChannelKey = std::tuple<int, int, int>;
using DetectorPlaneKey = std::pair<int, int>;
using GeometryMap = std::map<int, GeometryChannelMapping>;

struct PlaneQAHistograms {
   bool isPad{false};
   TDirectory* directory{nullptr};
   TH1D* stripOccupancy{nullptr};
   TH2D* padOccupancy{nullptr};
   TH1D* maximum{nullptr};
   TH1D* peakSample{nullptr};
   TH2D* maximumVsStrip{nullptr};
   TProfile2D* meanMaximumVsPad{nullptr};
};

struct PlaneQALayout {
   bool isPad{false};
   int minimumID0{0};
   int maximumID0{-1};
   int minimumID1{0};
   int maximumID1{-1};
};

std::string FormatCount(ULong64_t value) {
   std::string text = std::to_string(value);
   for (int position = static_cast<int>(text.size()) - 3; position > 0;
        position -= 3) {
      text.insert(static_cast<size_t>(position), ",");
   }
   return text;
}

bool LoadGeometry(const fs::path& path,
                  GeometryMap& geometry,
                  std::string& error) {
   std::ifstream input(path);
   if (!input) {
      error = "cannot open Geometry file: " + path.string();
      return false;
   }

   try {
      const json document = json::parse(input);
      if (!document.is_array()) {
         error = "Geometry root must be an array";
         return false;
      }

      geometry.clear();
      for (size_t index = 0; index < document.size(); ++index) {
         const auto& item = document[index];
         if (!item.is_object() ||
             !item.contains("detectorId") || !item["detectorId"].is_number_integer() ||
             !item.contains("chipId") || !item["chipId"].is_number_integer() ||
             !item.contains("planeId") || !item["planeId"].is_number_integer() ||
             !item.contains("id0") || !item["id0"].is_array()) {
            error = "invalid Geometry entry at index " + std::to_string(index);
            return false;
         }

         GeometryChannelMapping mapping;
         mapping.detectorID = item["detectorId"].get<int>();
         mapping.chipID = item["chipId"].get<int>();
         mapping.planeID = item["planeId"].get<int>();
         mapping.id0 = item["id0"].get<std::vector<int>>();
         if (mapping.id0.size() != 128) {
            error = "Geometry chip " + std::to_string(mapping.chipID) +
                    " has id0 length " + std::to_string(mapping.id0.size()) +
                    "; expected 128";
            return false;
         }

         if (item.contains("id1")) {
            if (!item["id1"].is_array()) {
               error = "Geometry chip " + std::to_string(mapping.chipID) +
                       " has non-array id1";
               return false;
            }
            mapping.id1 = item["id1"].get<std::vector<int>>();
            if (mapping.id1->size() != 128) {
               error = "Geometry chip " + std::to_string(mapping.chipID) +
                       " has id1 length " + std::to_string(mapping.id1->size()) +
                       "; expected 128";
               return false;
            }
         }

         if (!geometry.emplace(mapping.chipID, std::move(mapping)).second) {
            error = "duplicate Geometry chipId " +
                    std::to_string(item["chipId"].get<int>());
            return false;
         }
      }
   } catch (const std::exception& exception) {
      error = "invalid Geometry JSON: " + std::string(exception.what());
      return false;
   }

   if (geometry.empty()) {
      error = "Geometry contains no chip mappings";
      return false;
   }
   return true;
}
}  // namespace

bool APV25SRSConverter::AcquireRawData(const fs::path& rawDir,
                                       const std::string& runID,
                                       std::string& error) {
   const fs::path rawInput = rawDir / ("run" + runID + ".root");
   const fs::path geometryInput = rawDir.parent_path() / "Geometry.json";

   if (!fs::is_regular_file(rawInput)) {
      error = "conversion input does not exist: " + rawInput.string();
      return false;
   }
   if (!fs::is_regular_file(geometryInput)) {
      error = "Geometry input does not exist: " + geometryInput.string();
      return false;
   }

   m_inputPath = rawInput;
   m_geometryPath = geometryInput;
   return true;
}

bool APV25SRSConverter::Convert(const std::string& outputPath) {
   const std::string inputPath = m_inputPath.string();
   if (inputPath.empty()) {
      std::cerr << "[APV25SRSConverter] input is required\n";
      return false;
   }

   GeometryMap geometry;
   std::string geometryError;
   if (!LoadGeometry(m_geometryPath, geometry, geometryError)) {
      std::cerr << "[APV25SRSConverter] " << geometryError << '\n';
      return false;
   }

   std::map<DetectorPlaneKey, PlaneQALayout> qaLayouts;
   for (const auto& [chipID, mapping] : geometry) {
      const DetectorPlaneKey key{mapping.detectorID, mapping.planeID};
      const bool isPad = mapping.id1.has_value();
      auto [layoutIt, inserted] = qaLayouts.try_emplace(key);
      auto& layout = layoutIt->second;
      if (inserted) {
         layout.isPad = isPad;
      } else if (layout.isPad != isPad) {
         std::cerr << "[APV25SRSConverter] Detector " << mapping.detectorID
                   << " plane " << mapping.planeID
                   << " mixes strip and pad Geometry mappings\n";
         return false;
      }

      if (!isPad) continue;
      for (size_t channel = 0; channel < mapping.id0.size(); ++channel) {
         const int column = mapping.id0[channel];
         const int row = mapping.id1->at(channel);
         if (column < 0 || row < 0) continue;
         if (layout.maximumID0 < layout.minimumID0) {
            layout.minimumID0 = layout.maximumID0 = column;
            layout.minimumID1 = layout.maximumID1 = row;
         } else {
            layout.minimumID0 = std::min(layout.minimumID0, column);
            layout.maximumID0 = std::max(layout.maximumID0, column);
            layout.minimumID1 = std::min(layout.minimumID1, row);
            layout.maximumID1 = std::max(layout.maximumID1, row);
         }
      }
   }
   for (const auto& [key, layout] : qaLayouts) {
      if (layout.isPad && layout.maximumID0 < layout.minimumID0) {
         std::cerr << "[APV25SRSConverter] Detector " << key.first
                   << " plane " << key.second
                   << " has no valid pad coordinates in Geometry\n";
         return false;
      }
   }

   TFile inputFile(inputPath.c_str(), "READ");
   if (inputFile.IsZombie()) {
      std::cerr << "[APV25SRSConverter] Cannot open " << inputPath << '\n';
      return false;
   }
   TTree* rawTree = static_cast<TTree*>(inputFile.Get("raw"));
   if (!rawTree) {
      std::cerr << "[APV25SRSConverter] raw tree not found\n";
      return false;
   }

   const char* required[] = {"apv_evt", "apv_id", "apv_ch", "apv_q"};
   for (const char* branch : required) {
      if (!rawTree->GetBranch(branch)) {
         std::cerr << "[APV25SRSConverter] Missing raw branch: " << branch << '\n';
         return false;
      }
   }

   const bool hasTime = rawTree->GetBranch("time_s") && rawTree->GetBranch("time_us");
   unsigned int apvEvent = 0;
   int timeSeconds = 0, timeMicroseconds = 0;
   std::vector<unsigned int>* apvIDs = nullptr;
   std::vector<unsigned int>* apvChannels = nullptr;
   std::vector<std::vector<short>>* waveformsRaw = nullptr;

   rawTree->SetBranchAddress("apv_evt", &apvEvent);
   rawTree->SetBranchAddress("apv_id", &apvIDs);
   rawTree->SetBranchAddress("apv_ch", &apvChannels);
   rawTree->SetBranchAddress("apv_q", &waveformsRaw);
   if (hasTime) {
      rawTree->SetBranchAddress("time_s", &timeSeconds);
      rawTree->SetBranchAddress("time_us", &timeMicroseconds);
   }

   fs::path output(outputPath);
   if (output.has_parent_path()) fs::create_directories(output.parent_path());
   TFile outputFile(outputPath.c_str(), "RECREATE");
   if (outputFile.IsZombie()) {
      std::cerr << "[APV25SRSConverter] Cannot create " << outputPath << '\n';
      return false;
   }

   ULong64_t eventID = 0, timestamp = 0;
   std::vector<int> detectorIDs, planeTypes, id0s, id1s;
   std::vector<int> boardIDs, chipIDs, channelIDs;
   std::vector<short> waveformMaximums;
   std::vector<int> peakSamples;
   std::vector<std::vector<short>> waveforms;
   std::vector<unsigned int> channelFlags;

   TTree eventTree("Events", "BeamAnalysis canonical raw events");
   eventTree.Branch("event_id", &eventID);
   eventTree.Branch("timestamp", &timestamp);
   eventTree.Branch("detector_id", &detectorIDs);
   eventTree.Branch("plane_type", &planeTypes);
   eventTree.Branch("id0", &id0s);
   eventTree.Branch("id1", &id1s);
   eventTree.Branch("board_id", &boardIDs);
   eventTree.Branch("chip_id", &chipIDs);
   eventTree.Branch("channel_id", &channelIDs);
   eventTree.Branch("waveform", &waveforms);
   eventTree.Branch("waveform_max", &waveformMaximums);
   eventTree.Branch("peak_sample", &peakSamples);
   eventTree.Branch("channel_flags", &channelFlags);

   TDirectory* qaDirectory = outputFile.mkdir("QA");
   TDirectory* globalDirectory = qaDirectory->mkdir("Global");
   globalDirectory->cd();
   TH1D globalMultiplicity(
       "hEventChannelMultiplicity",
       "Channels per event;Channels/event;Events", 100, -0.5, 99.5);
   TH1D globalMaximum(
       "hWaveformMaximum",
       "Raw waveform maximum;Maximum [ADC];Channel waveforms", 2048, 0, 2048);
   TH1D globalPeakSample(
       "hPeakSample",
       "Waveform peak sample;Sample index;Channel waveforms", 30, 0, 30);

   globalMultiplicity.SetDirectory(nullptr);
   globalMaximum.SetDirectory(nullptr);
   globalPeakSample.SetDirectory(nullptr);
   outputFile.cd();

   std::map<DetectorPlaneKey, PlaneQAHistograms> planeQA;
   std::map<ChannelKey, APV25SRSChannelInfo> observedChannels;
   const auto entries = rawTree->GetEntries();
   const bool hasPedestals = inputFile.Get("pedestals") != nullptr;
   constexpr Long64_t progressInterval = 10000;
   const auto conversionStarted = std::chrono::steady_clock::now();
   ULong64_t waveformCount = 0;
   ULong64_t emptyWaveformCount = 0;
   size_t waveformLength = 0;
   short minimumObservedMaximum = std::numeric_limits<short>::max();
   short maximumObservedMaximum = std::numeric_limits<short>::min();
   std::set<int> observedDetectors;

   if (Terminal::Verbose()) {
      Terminal::Detail("input " + inputPath);
      Terminal::Detail("output " + outputPath);
   }

   for (Long64_t entry = 0; entry < entries; ++entry) {
      rawTree->GetEntry(entry);
      const size_t count = apvIDs ? apvIDs->size() : 0;
      if (!apvChannels || !waveformsRaw ||
          apvChannels->size() != count || waveformsRaw->size() != count) {
         std::cerr << "[APV25SRSConverter] Inconsistent vectors at entry " << entry << '\n';
         return false;
      }
      eventID = apvEvent;
      timestamp = hasTime
                      ? static_cast<ULong64_t>(timeSeconds) * 1000000ULL +
                            static_cast<ULong64_t>(timeMicroseconds)
                      : static_cast<ULong64_t>(entry);
      detectorIDs.clear();
      planeTypes.clear();
      id0s.clear();
      id1s.clear();
      boardIDs.clear();
      chipIDs.clear();
      channelIDs.clear();
      waveformMaximums.clear();
      peakSamples.clear();
      waveforms.clear();
      channelFlags.clear();

      globalMultiplicity.Fill(static_cast<double>(count));

      for (size_t i = 0; i < count; ++i) {
         const int boardID = static_cast<int>((*apvIDs)[i]);
         const int channelID = static_cast<int>((*apvChannels)[i]);
         const auto geometryIt = geometry.find(boardID);
         if (geometryIt == geometry.end()) {
            std::cerr << "[APV25SRSConverter] Geometry is missing chipId "
                      << boardID << " at raw entry " << entry << '\n';
            return false;
         }
         if (channelID < 0 || channelID >= 128) {
            std::cerr << "[APV25SRSConverter] APV channel " << channelID
                      << " is outside [0, 127] for chipId " << boardID
                      << " at raw entry " << entry << '\n';
            return false;
         }
         const auto& mapping = geometryIt->second;
         const int detectorID = mapping.detectorID;
         const int planeType = mapping.planeID;
         const int id0 = mapping.id0[static_cast<size_t>(channelID)];
         const int id1 = mapping.id1
                             ? mapping.id1->at(static_cast<size_t>(channelID))
                             : -1;
         detectorIDs.push_back(detectorID);
         planeTypes.push_back(planeType);
         id0s.push_back(id0);
         id1s.push_back(id1);
         boardIDs.push_back(boardID);
         chipIDs.push_back(0);
         channelIDs.push_back(channelID);
         const auto& waveform = (*waveformsRaw)[i];
         waveforms.push_back(waveform);

         short waveformMaximum = 0;
         int peakSample = -1;
         if (waveform.empty()) {
            ++emptyWaveformCount;
         } else {
            if (waveformLength == 0) waveformLength = waveform.size();
            const auto maximum = std::max_element(waveform.begin(), waveform.end());
            waveformMaximum = *maximum;
            peakSample = static_cast<int>(std::distance(waveform.begin(), maximum));
            globalMaximum.Fill(waveformMaximum);
            globalPeakSample.Fill(peakSample);
            minimumObservedMaximum = std::min(minimumObservedMaximum, waveformMaximum);
            maximumObservedMaximum = std::max(maximumObservedMaximum, waveformMaximum);
         }
         waveformMaximums.push_back(waveformMaximum);
         peakSamples.push_back(peakSample);
         ++waveformCount;

         const DetectorPlaneKey qaKey{detectorID, planeType};
         const auto layoutIt = qaLayouts.find(qaKey);
         const bool isPad = layoutIt != qaLayouts.end() && layoutIt->second.isPad;
         const bool hasPhysicalCoordinate = id0 >= 0 && (!isPad || id1 >= 0);
         if (hasPhysicalCoordinate) {
            auto [qaIt, inserted] = planeQA.try_emplace(qaKey);
            if (inserted) {
               qaDirectory->cd();
               const std::string detectorName = "Detector_" + std::to_string(detectorID);
               TDirectory* detectorDirectory = qaDirectory->GetDirectory(detectorName.c_str());
               if (!detectorDirectory)
                  detectorDirectory = qaDirectory->mkdir(detectorName.c_str());
               const std::string planeName = "Plane_" + std::to_string(planeType);
               TDirectory* planeDirectory = detectorDirectory->mkdir(planeName.c_str());
               planeDirectory->cd();

               auto& qa = qaIt->second;
               qa.isPad = isPad;
               qa.directory = planeDirectory;
               if (isPad) {
                  const auto& layout = layoutIt->second;
                  const int columns = layout.maximumID0 - layout.minimumID0 + 1;
                  const int rows = layout.maximumID1 - layout.minimumID1 + 1;
                  qa.padOccupancy = new TH2D(
                      "hOccupancyPad", "Pad occupancy;Column;Row;Entries",
                      columns, layout.minimumID0 - 0.5, layout.maximumID0 + 0.5,
                      rows, layout.minimumID1 - 0.5, layout.maximumID1 + 0.5);
                  qa.meanMaximumVsPad = new TProfile2D(
                      "hMeanWaveformMaximumPad",
                      "Mean raw waveform maximum per pad;Column;Row;Mean maximum [ADC]",
                      columns, layout.minimumID0 - 0.5, layout.maximumID0 + 0.5,
                      rows, layout.minimumID1 - 0.5, layout.maximumID1 + 0.5);
               } else {
                  qa.stripOccupancy = new TH1D(
                      "hOccupancy", "Channel occupancy;Strip ID;Entries",
                      256, -0.5, 255.5);
                  qa.maximumVsStrip = new TH2D(
                      "hMaximumVsStrip",
                      "Raw maximum versus strip;Strip ID;Maximum [ADC]", 256, -0.5,
                      255.5, 2048, 0, 2048);
               }
               qa.maximum = new TH1D(
                   "hWaveformMaximum",
                   "Raw waveform maximum;Maximum [ADC];Channel waveforms", 2048, 0, 2048);
               qa.peakSample = new TH1D(
                   "hPeakSample", "Waveform peak sample;Sample index;Channel waveforms",
                   30, 0, 30);
               outputFile.cd();
            }
            auto& qa = qaIt->second;
            if (qa.isPad)
               qa.padOccupancy->Fill(id0, id1);
            else
               qa.stripOccupancy->Fill(id0);
            if (!waveform.empty()) {
               qa.maximum->Fill(waveformMaximum);
               qa.peakSample->Fill(peakSample);
               if (qa.isPad)
                  qa.meanMaximumVsPad->Fill(id0, id1, waveformMaximum);
               else
                  qa.maximumVsStrip->Fill(id0, waveformMaximum);
            }
         }
         observedDetectors.insert(detectorID);
         channelFlags.push_back(0);
         observedChannels[{boardID, channelID, id0}] =
             {boardID, 0, channelID, detectorID, planeType, id0, id1, -1.0};
      }
      eventTree.Fill();

      const Long64_t processed = entry + 1;
      if (Terminal::Interactive() &&
          (processed % progressInterval == 0 || processed == entries)) {
         const double elapsed = std::chrono::duration<double>(
                                    std::chrono::steady_clock::now() - conversionStarted)
                                    .count();
         const double rate = elapsed > 0.0 ? processed / elapsed : 0.0;
         const double eta = rate > 0.0 ? (entries - processed) / rate : 0.0;
         const double percent = entries > 0 ? 100.0 * processed / entries : 100.0;
         std::ostringstream progress;
         progress << "      " << FormatCount(processed) << " / "
                  << FormatCount(entries) << " (" << std::fixed
                  << std::setprecision(1) << percent << "%)  "
                  << std::setprecision(0) << rate << " evt/s  ETA "
                  << std::setprecision(1) << eta << " s";
         std::cout << '\r' << std::left << std::setw(88) << progress.str()
                   << std::right << std::flush;
      }
   }
   Terminal::ClearProgress();

   TTree* pedestalTree = static_cast<TTree*>(inputFile.Get("pedestals"));
   if (pedestalTree && pedestalTree->GetEntries() > 0) {
      std::vector<unsigned int>* pedAPVIDs = nullptr;
      std::vector<unsigned int>* pedChannels = nullptr;
      std::vector<double>* pedStd = nullptr;
      const char* pedestalRequired[] = {"apv_id", "apv_ch", "apv_pedstd"};
      for (const char* branch : pedestalRequired) {
         if (!pedestalTree->GetBranch(branch)) {
            std::cerr << "[APV25SRSConverter] Missing pedestal branch: "
                      << branch << '\n';
            return false;
         }
      }
      pedestalTree->SetBranchAddress("apv_id", &pedAPVIDs);
      pedestalTree->SetBranchAddress("apv_ch", &pedChannels);
      pedestalTree->SetBranchAddress("apv_pedstd", &pedStd);
      pedestalTree->GetEntry(0);
      if (!pedAPVIDs || !pedChannels || !pedStd ||
          pedChannels->size() != pedAPVIDs->size() ||
          pedStd->size() != pedAPVIDs->size()) {
         std::cerr << "[APV25SRSConverter] Inconsistent pedestal vectors\n";
         return false;
      }
      for (size_t i = 0; i < pedAPVIDs->size(); ++i) {
         const int boardID = static_cast<int>((*pedAPVIDs)[i]);
         const int channelID = static_cast<int>((*pedChannels)[i]);
         const auto geometryIt = geometry.find(boardID);
         if (geometryIt == geometry.end()) {
            std::cerr << "[APV25SRSConverter] Geometry is missing pedestal chipId "
                      << boardID << '\n';
            return false;
         }
         if (channelID < 0 || channelID >= 128) {
            std::cerr << "[APV25SRSConverter] Pedestal APV channel "
                      << channelID << " is outside [0, 127] for chipId "
                      << boardID << '\n';
            return false;
         }
         const auto& mapping = geometryIt->second;
         const int id0 = mapping.id0[static_cast<size_t>(channelID)];
         const int id1 = mapping.id1
                             ? mapping.id1->at(static_cast<size_t>(channelID))
                             : -1;
         auto& info = observedChannels[{boardID, channelID, id0}];
         info.boardID = boardID;
         info.channelID = channelID;
         info.id0 = id0;
         info.id1 = id1;
         info.noiseSigma = (*pedStd)[i];
         info.detectorID = mapping.detectorID;
         info.planeType = mapping.planeID;
      }
   }

   int mapBoardID, mapChipID, mapChannelID;
   int detectorID, planeType, id0, id1, polarity = 1;
   double pedestal = 0.0, noiseSigma, gain = 1.0;
   unsigned int status = 0;
   TTree channelTree("Channels", "BeamAnalysis channel map and calibration");
   channelTree.Branch("board_id", &mapBoardID);
   channelTree.Branch("chip_id", &mapChipID);
   channelTree.Branch("channel_id", &mapChannelID);
   channelTree.Branch("detector_id", &detectorID);
   channelTree.Branch("plane_type", &planeType);
   channelTree.Branch("id0", &id0);
   channelTree.Branch("id1", &id1);
   channelTree.Branch("pedestal", &pedestal);
   channelTree.Branch("noise_sigma", &noiseSigma);
   channelTree.Branch("gain", &gain);
   channelTree.Branch("polarity", &polarity);
   channelTree.Branch("status", &status);
   for (const auto& [key, info] : observedChannels) {
      mapBoardID = info.boardID;
      mapChipID = info.chipID;
      mapChannelID = info.channelID;
      detectorID = info.detectorID;
      planeType = info.planeType;
      id0 = info.id0;
      id1 = info.id1;
      noiseSigma = info.noiseSigma;
      channelTree.Fill();
   }

   int schemaVersion = 2;
   std::string sourceFormat = "APV25_SRS_ROOT";
   std::string converterVersion = "2.0.0";
   std::string geometrySource = m_geometryPath.string();
   TTree metadataTree("Metadata", "BeamAnalysis canonical format metadata");
   metadataTree.Branch("schema_version", &schemaVersion);
   metadataTree.Branch("source_format", &sourceFormat);
   metadataTree.Branch("converter_version", &converterVersion);
   metadataTree.Branch("geometry_source", &geometrySource);
   metadataTree.Fill();

   eventTree.Write();
   channelTree.Write();
   metadataTree.Write();
   globalDirectory->cd();
   globalMultiplicity.Write();
   globalMaximum.Write();
   globalPeakSample.Write();
   for (const auto& [key, qa] : planeQA) {
      qa.directory->cd();
      if (qa.isPad) {
         qa.padOccupancy->Write();
         qa.meanMaximumVsPad->Write();
      } else {
         qa.stripOccupancy->Write();
         qa.maximumVsStrip->Write();
      }
      qa.maximum->Write();
      qa.peakSample->Write();
   }
   outputFile.Close();
   inputFile.Close();

   const double elapsed = std::chrono::duration<double>(
                              std::chrono::steady_clock::now() - conversionStarted)
                              .count();
   std::ostringstream summary;
   summary << FormatCount(entries) << " events · "
           << FormatCount(waveformCount) << " waveforms · "
           << FormatCount(static_cast<Long64_t>(observedChannels.size()))
           << " channels · " << std::fixed << std::setprecision(0)
           << (elapsed > 0.0 ? entries / elapsed : 0.0) << " evt/s";
   Terminal::Detail(summary.str());
   if (emptyWaveformCount > 0) {
      Terminal::Note(FormatCount(emptyWaveformCount) + " empty waveforms");
   }
   if (Terminal::Verbose()) {
      std::ostringstream detail;
      detail << observedDetectors.size() << " detectors · "
             << waveformLength << " samples/waveform · pedestals "
             << (hasPedestals ? "available" : "unavailable");
      if (waveformCount > emptyWaveformCount) {
         detail << " · ADC max " << minimumObservedMaximum << ".."
                << maximumObservedMaximum;
      }
      Terminal::Detail(Terminal::Muted(detail.str()));
   }
   return true;
}
