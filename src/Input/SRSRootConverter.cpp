#include "Input/SRSRootConverter.h"

#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <tuple>
#include <vector>

namespace fs = std::filesystem;

namespace {
struct SRSChannelInfo {
   int boardID{0};
   int chipID{0};
   int channelID{0};
   int detectorID{0};
   int planeType{0};
   int stripID{0};
   double noiseSigma{-1.0};
};

using ChannelKey = std::tuple<int, int, int>;
using DetectorPlaneKey = std::pair<int, int>;

struct PlaneQAHistograms {
   TH1D* occupancy{nullptr};
   TH1D* maximum{nullptr};
   TH1D* peakSample{nullptr};
   TH2D* maximumVsStrip{nullptr};
};

std::string FormatCount(ULong64_t value) {
   std::string text = std::to_string(value);
   for (int position = static_cast<int>(text.size()) - 3; position > 0;
        position -= 3) {
      text.insert(static_cast<size_t>(position), ",");
   }
   return text;
}

std::tuple<int, int, int> MapSRSChannel(unsigned int boardID,
                                        unsigned int channelID,
                                        unsigned int mmStrip) {
   constexpr std::array<int, 16> boardToRawIndex = {
       0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7};
   const int rawDataIndex = boardID < boardToRawIndex.size()
                                ? boardToRawIndex[boardID]
                                : static_cast<int>(boardID) / 2;
   const int planeType = rawDataIndex % 2 == 0 ? 0 : 1;
   const int detectorID = rawDataIndex / 2 + 1;
   int stripID = static_cast<int>(mmStrip);

   if (boardID == 12)
      stripID = 256 - static_cast<int>(channelID);
   else if (boardID == 13)
      stripID = 128 - static_cast<int>(channelID);

   if (boardID == 14)
      stripID = 256 - static_cast<int>(channelID);
   else if (boardID == 15)
      stripID = 128 - static_cast<int>(channelID);

   return {detectorID, stripID, planeType};
}
}  // namespace

bool SRSRootConverter::AcquireRawData(const fs::path& rawDir,
                                      const std::string& runID,
                                      std::string& error) {
   const fs::path rawInput = rawDir / ("run" + runID + ".root");

   if (!fs::is_regular_file(rawInput)) {
      error = "conversion input does not exist: " + rawInput.string();
      return false;
   }

   m_inputPath = rawInput;
   return true;
}

bool SRSRootConverter::Convert(const std::string& outputPath) {
   const std::string inputPath = m_inputPath.string();
   if (inputPath.empty()) {
      std::cerr << "[SRSRootConverter] input is required\n";
      return false;
   }

   TFile inputFile(inputPath.c_str(), "READ");
   if (inputFile.IsZombie()) {
      std::cerr << "[SRSRootConverter] Cannot open " << inputPath << '\n';
      return false;
   }
   TTree* rawTree = static_cast<TTree*>(inputFile.Get("raw"));
   if (!rawTree) {
      std::cerr << "[SRSRootConverter] raw tree not found\n";
      return false;
   }

   const char* required[] = {"apv_evt", "apv_id", "apv_ch", "mm_strip", "apv_q"};
   for (const char* branch : required) {
      if (!rawTree->GetBranch(branch)) {
         std::cerr << "[SRSRootConverter] Missing raw branch: " << branch << '\n';
         return false;
      }
   }

   const bool hasTime = rawTree->GetBranch("time_s") && rawTree->GetBranch("time_us");
   unsigned int apvEvent = 0;
   int timeSeconds = 0, timeMicroseconds = 0;
   std::vector<unsigned int>* apvIDs = nullptr;
   std::vector<unsigned int>* apvChannels = nullptr;
   std::vector<unsigned int>* strips = nullptr;
   std::vector<std::vector<short>>* waveformsRaw = nullptr;

   rawTree->SetBranchAddress("apv_evt", &apvEvent);
   rawTree->SetBranchAddress("apv_id", &apvIDs);
   rawTree->SetBranchAddress("apv_ch", &apvChannels);
   rawTree->SetBranchAddress("mm_strip", &strips);
   rawTree->SetBranchAddress("apv_q", &waveformsRaw);
   if (hasTime) {
      rawTree->SetBranchAddress("time_s", &timeSeconds);
      rawTree->SetBranchAddress("time_us", &timeMicroseconds);
   }

   fs::path output(outputPath);
   if (output.has_parent_path()) fs::create_directories(output.parent_path());
   TFile outputFile(outputPath.c_str(), "RECREATE");
   if (outputFile.IsZombie()) {
      std::cerr << "[SRSRootConverter] Cannot create " << outputPath << '\n';
      return false;
   }

   ULong64_t eventID = 0, timestamp = 0;
   std::vector<int> detectorIDs, planeTypes, stripIDs;
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
   eventTree.Branch("id", &stripIDs);
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
   std::map<ChannelKey, SRSChannelInfo> observedChannels;
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

   std::cout << "  Input:             " << inputPath << '\n'
             << "  Output:            " << outputPath << '\n'
             << "  Events:            " << FormatCount(entries) << '\n'
             << "  Pedestals:         " << (hasPedestals ? "available" : "unavailable")
             << "\n------------------------------------------------------------\n";

   for (Long64_t entry = 0; entry < entries; ++entry) {
      rawTree->GetEntry(entry);
      const size_t count = apvIDs ? apvIDs->size() : 0;
      if (!apvChannels || !strips || !waveformsRaw ||
          apvChannels->size() != count || strips->size() != count || waveformsRaw->size() != count) {
         std::cerr << "[SRSRootConverter] Inconsistent vectors at entry " << entry << '\n';
         return false;
      }
      eventID = apvEvent;
      timestamp = hasTime
                      ? static_cast<ULong64_t>(timeSeconds) * 1000000ULL +
                            static_cast<ULong64_t>(timeMicroseconds)
                      : static_cast<ULong64_t>(entry);
      detectorIDs.clear();
      planeTypes.clear();
      stripIDs.clear();
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
         auto [detectorID, stripID, planeType] =
             MapSRSChannel((*apvIDs)[i], (*apvChannels)[i], (*strips)[i]);
         detectorIDs.push_back(detectorID);
         planeTypes.push_back(planeType);
         stripIDs.push_back(stripID);
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
            qaIt->second.occupancy = new TH1D(
                "hOccupancy", "Channel occupancy;Strip ID;Entries", 256, -0.5, 255.5);
            qaIt->second.maximum = new TH1D(
                "hWaveformMaximum",
                "Raw waveform maximum;Maximum [ADC];Channel waveforms", 2048, 0, 2048);
            qaIt->second.peakSample = new TH1D(
                "hPeakSample", "Waveform peak sample;Sample index;Channel waveforms", 30, 0, 30);
            qaIt->second.maximumVsStrip = new TH2D(
                "hMaximumVsStrip",
                "Raw maximum versus strip;Strip ID;Maximum [ADC]", 256, -0.5,
                255.5, 2048, 0, 2048);
            outputFile.cd();
         }
         auto& qa = qaIt->second;
         qa.occupancy->Fill(stripID);
         if (!waveform.empty()) {
            qa.maximum->Fill(waveformMaximum);
            qa.peakSample->Fill(peakSample);
            qa.maximumVsStrip->Fill(stripID, waveformMaximum);
         }
         observedDetectors.insert(detectorID);
         channelFlags.push_back(0);
         observedChannels[{boardID, channelID, stripID}] =
             {boardID, 0, channelID, detectorID, planeType, stripID, -1.0};
      }
      eventTree.Fill();

      const Long64_t processed = entry + 1;
      if (processed % progressInterval == 0 || processed == entries) {
         const double elapsed = std::chrono::duration<double>(
                                    std::chrono::steady_clock::now() - conversionStarted)
                                    .count();
         const double rate = elapsed > 0.0 ? processed / elapsed : 0.0;
         const double eta = rate > 0.0 ? (entries - processed) / rate : 0.0;
         const double percent = entries > 0 ? 100.0 * processed / entries : 100.0;
         std::ostringstream progress;
         progress << "  Progress: " << FormatCount(processed) << " / "
                  << FormatCount(entries) << " (" << std::fixed
                  << std::setprecision(1) << percent << "%)  "
                  << std::setprecision(0) << rate << " evt/s  ETA "
                  << std::setprecision(1) << eta << " s";
         std::cout << '\r' << std::left << std::setw(88) << progress.str()
                   << std::right;
         if (processed == entries)
            std::cout << '\n';
         else
            std::cout << std::flush;
      }
   }

   TTree* pedestalTree = static_cast<TTree*>(inputFile.Get("pedestals"));
   if (pedestalTree && pedestalTree->GetEntries() > 0) {
      std::vector<unsigned int>* pedAPVIDs = nullptr;
      std::vector<unsigned int>* pedChannels = nullptr;
      std::vector<unsigned int>* pedStrips = nullptr;
      std::vector<double>* pedStd = nullptr;
      pedestalTree->SetBranchAddress("apv_id", &pedAPVIDs);
      pedestalTree->SetBranchAddress("apv_ch", &pedChannels);
      pedestalTree->SetBranchAddress("mm_strip", &pedStrips);
      pedestalTree->SetBranchAddress("apv_pedstd", &pedStd);
      pedestalTree->GetEntry(0);
      if (pedAPVIDs && pedChannels && pedStrips && pedStd) {
         for (size_t i = 0; i < pedAPVIDs->size(); ++i) {
            const int boardID = static_cast<int>((*pedAPVIDs)[i]);
            const int channelID = static_cast<int>((*pedChannels)[i]);
            auto [detectorID, stripID, planeType] =
                MapSRSChannel((*pedAPVIDs)[i], (*pedChannels)[i], (*pedStrips)[i]);
            auto& info = observedChannels[{boardID, channelID, stripID}];
            info.boardID = boardID;
            info.channelID = channelID;
            info.stripID = stripID;
            info.noiseSigma = (*pedStd)[i];
            info.detectorID = detectorID;
            info.planeType = planeType;
         }
      }
   }

   int mapBoardID, mapChipID, mapChannelID;
   int detectorID, planeType, stripID, polarity = 1;
   double pedestal = 0.0, noiseSigma, gain = 1.0;
   unsigned int status = 0;
   TTree channelTree("Channels", "BeamAnalysis channel map and calibration");
   channelTree.Branch("board_id", &mapBoardID);
   channelTree.Branch("chip_id", &mapChipID);
   channelTree.Branch("channel_id", &mapChannelID);
   channelTree.Branch("detector_id", &detectorID);
   channelTree.Branch("plane_type", &planeType);
   channelTree.Branch("id", &stripID);
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
      stripID = info.stripID;
      noiseSigma = info.noiseSigma;
      channelTree.Fill();
   }

   eventTree.Write();
   channelTree.Write();
   globalDirectory->cd();
   globalMultiplicity.Write();
   globalMaximum.Write();
   globalPeakSample.Write();
   for (const auto& [key, qa] : planeQA) {
      qa.occupancy->GetDirectory()->cd();
      qa.occupancy->Write();
      qa.maximum->Write();
      qa.peakSample->Write();
      qa.maximumVsStrip->Write();
   }
   outputFile.Close();
   inputFile.Close();

   const double elapsed = std::chrono::duration<double>(
                              std::chrono::steady_clock::now() - conversionStarted)
                              .count();
   const size_t qaHistogramCount = 3 + 4 * planeQA.size();
   std::cout << "------------------------------------------------------------\n"
             << "  Conversion summary\n"
             << "  Events written:      " << FormatCount(entries) << '\n'
             << "  Channel waveforms:   " << FormatCount(waveformCount) << '\n'
             << "  Detectors observed:  " << observedDetectors.size() << '\n'
             << "  Active channels:     " << observedChannels.size() << '\n';
   if (waveformLength > 0)
      std::cout << "  Waveform length:     " << waveformLength << " samples\n";
   if (waveformCount > emptyWaveformCount) {
      std::cout << "  Raw maximum range:   " << minimumObservedMaximum << " .. "
                << maximumObservedMaximum << " ADC\n";
   }
   if (emptyWaveformCount > 0)
      std::cout << "  Empty waveforms:     " << FormatCount(emptyWaveformCount) << '\n';
   std::cout << "  QA histograms:       " << qaHistogramCount << '\n'
             << "  Elapsed:             " << std::fixed << std::setprecision(1)
             << elapsed << " s\n"
             << "  Average rate:        " << std::setprecision(0)
             << (elapsed > 0.0 ? entries / elapsed : 0.0) << " events/s\n";
   return true;
}
