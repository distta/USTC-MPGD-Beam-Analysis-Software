#include "Input/VMM3aConverter.h"

#include <TFile.h>
#include <TError.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

namespace fs = std::filesystem;
using json = nlohmann::json;

namespace {

using GeometryKey = std::tuple<int, int, int>;

ErrorHandlerFunc_t previousRootErrorHandler = nullptr;

void VMMInputErrorHandler(int level, Bool_t abort,
                          const char* location, const char* message) {
   const std::string_view where =
       location ? std::string_view(location) : std::string_view();
   const std::string_view text =
       message ? std::string_view(message) : std::string_view();
   const bool unusedClassWarning =
       level == kWarning && where == "TClass::Init" &&
       (text == "no dictionary for class Hit is available" ||
        text == "no dictionary for class ClusterPlane is available" ||
        text == "no dictionary for class ClusterDetector is available");
   if (unusedClassWarning) return;
   if (previousRootErrorHandler)
      previousRootErrorHandler(level, abort, location, message);
}

struct GeometryEntry {
   int detector{0};
   int plane{0};
   int fec{0};
   int chip{0};
   std::vector<int> id0;
   std::vector<int> id1;
};

struct VMMHit {
   unsigned int rawHitID{0};
   int detector{0};
   int plane{0};
   int id0{-1};
   int id1{-1};
   int chip{-1};
   int channel{-1};
   unsigned short adc{0};
   double timeNs{0.0};
};

struct TriggerCandidate {
   unsigned int rawHitID{0};
   unsigned short adc{0};
   double rawTimeNs{0.0};
   double correctedTimeNs{0.0};
};

struct Trigger {
   unsigned long long id{0};
   double timeNs{0.0};
};

struct ChannelRow {
   int chip{-1};
   int channel{-1};
   int detector{-1};
   int plane{-1};
   int id0{-1};
   int id1{-1};

   auto Key() const {
      return std::tie(detector, plane, id0, id1, chip, channel);
   }
};

std::string ReplaceRunID(std::string value, const std::string& runID) {
   constexpr const char* token = "{run_id}";
   size_t position = 0;
   while ((position = value.find(token, position)) != std::string::npos) {
      value.replace(position, std::char_traits<char>::length(token), runID);
      position += runID.size();
   }
   return value;
}

fs::path ResolvePath(const fs::path& rawDir, const std::string& value) {
   fs::path path(value);
   return path.is_absolute() ? path.lexically_normal()
                            : (rawDir / path).lexically_normal();
}

bool LoadGeometry(const fs::path& path,
                  std::map<GeometryKey, GeometryEntry>& geometry,
                  std::string& error) {
   std::ifstream input(path);
   if (!input) {
      error = "cannot open VMM geometry: " + path.string();
      return false;
   }

   try {
      const json document = json::parse(input);
      if (!document.is_object() ||
          !document.contains("vmm_geometry") ||
          !document["vmm_geometry"].is_array()) {
         error = "geometry.json must contain a vmm_geometry array";
         return false;
      }

      for (size_t index = 0; index < document["vmm_geometry"].size();
           ++index) {
         const auto& item = document["vmm_geometry"][index];
         if (!item.is_object() ||
             !item.contains("detector") ||
             !item.contains("plane") ||
             !item.contains("fec") ||
             !item.contains("vmm") ||
             !item.contains("id0")) {
            error = "invalid vmm_geometry entry " + std::to_string(index);
            return false;
         }

         GeometryEntry entry;
         entry.detector = item["detector"].get<int>();
         entry.plane = item["plane"].get<int>();
         entry.fec = item["fec"].get<int>();
         entry.chip = item["vmm"].get<int>();
         entry.id0 = item["id0"].get<std::vector<int>>();
         if (entry.id0.size() != 64) {
            error = "vmm_geometry entry " + std::to_string(index) +
                    " has id0 size " + std::to_string(entry.id0.size()) +
                    "; expected 64";
            return false;
         }
         entry.id1.assign(64, -1);
         if (item.contains("id1")) {
            entry.id1 = item["id1"].get<std::vector<int>>();
            if (entry.id1.size() != 64) {
               error = "vmm_geometry entry " + std::to_string(index) +
                       " has id1 size " +
                       std::to_string(entry.id1.size()) +
                       "; expected 64";
               return false;
            }
         }

         const GeometryKey key{
             entry.detector, entry.fec, entry.chip};
         if (!geometry.emplace(key, std::move(entry)).second) {
            error = "duplicate (detector,fec,vmm) VMM geometry key at "
                    "entry " + std::to_string(index);
            return false;
         }
      }
   } catch (const std::exception& exception) {
      error = "invalid VMM geometry JSON: " +
              std::string(exception.what());
      return false;
   }

   return true;
}

}  // namespace

bool VMM3aConverter::AcquireRawData(const fs::path& rawDir,
                                    const std::string& runID,
                                    std::string& error) {
   const std::string inputName = ReplaceRunID(
       m_config.value("inputFile", "run{run_id}_.root"), runID);
   const std::string geometryName =
       m_config.value("geometryFile", "geometry.json");
   m_inputPath = ResolvePath(rawDir, inputName);
   m_geometryPath = ResolvePath(rawDir.parent_path(), geometryName);
   if (!fs::path(geometryName).is_absolute() &&
       !fs::is_regular_file(m_geometryPath)) {
      const fs::path legacyPath = ResolvePath(rawDir, geometryName);
      if (fs::is_regular_file(legacyPath))
         m_geometryPath = legacyPath;
   }

   if (!fs::is_regular_file(m_inputPath)) {
      error = "VMM3a input does not exist: " + m_inputPath.string();
      return false;
   }
   if (!fs::is_regular_file(m_geometryPath)) {
      error = "VMM3a geometry does not exist: " +
              m_geometryPath.string();
      return false;
   }
   return true;
}

bool VMM3aConverter::Convert(const std::string& outputPath) {
   std::map<GeometryKey, GeometryEntry> geometry;
   std::string error;
   if (!LoadGeometry(m_geometryPath, geometry, error)) {
      std::cerr << "[VMM3aConverter] " << error << '\n';
      return false;
   }

   previousRootErrorHandler = GetErrorHandler();
   SetErrorHandler(VMMInputErrorHandler);
   TFile input(m_inputPath.c_str(), "READ");
   SetErrorHandler(previousRootErrorHandler);
   previousRootErrorHandler = nullptr;
   if (input.IsZombie()) {
      std::cerr << "[VMM3aConverter] cannot open " << m_inputPath << '\n';
      return false;
   }

   const std::string treeName = m_config.value("inputTree", "hits");
   const std::string timeBranch = m_config.value("timeBranch", "time");
   auto* source = static_cast<TTree*>(input.Get(treeName.c_str()));
   if (!source) {
      std::cerr << "[VMM3aConverter] tree " << treeName
                << " was not found\n";
      return false;
   }
   const std::vector<std::string> requiredBranches{
       "id", "det", "plane", "fec", "vmm", "ch", "adc", timeBranch};
   for (const auto& branch : requiredBranches) {
      if (!source->GetBranch(branch.c_str())) {
         std::cerr << "[VMM3aConverter] missing branch " << branch << '\n';
         return false;
      }
   }

   const json builder =
       m_config.value("eventBuilder", json::object());
   const json triggerConfig =
       builder.value("trigger", json::object());
   const json windowConfig =
       builder.value("timeWindow", json::object());

   const int triggerDetector =
       triggerConfig.value("detector", 0);
   const auto triggerChannels =
       triggerConfig.value("channels", std::vector<int>{3, 7});
   const auto triggerOffsets =
       triggerConfig.value("timeOffsetsNs",
                           std::vector<double>{0.0, -3.958});
   if (triggerChannels.size() != 2 ||
       triggerOffsets.size() != triggerChannels.size() ||
       triggerChannels[0] == triggerChannels[1]) {
      std::cerr << "[VMM3aConverter] trigger channels/timeOffsetsNs must "
                   "contain two entries\n";
      return false;
   }

   const unsigned int triggerADC =
       triggerConfig.value("adcThreshold", 300U);
   const double coincidenceWindow =
       triggerConfig.value("coincidenceWindowNs", 50.0);
   const double deadTime =
       triggerConfig.value("deadTimeNs", 0.0);
   const double windowOffset =
       windowConfig.value("offsetNs", 0.0);
   const double windowBefore =
       windowConfig.value("beforeNs", 250.0);
   const double windowAfter =
       windowConfig.value("afterNs", 250.0);
   const unsigned int hitADCThreshold =
       builder.value("hitADCThreshold", 0U);
   const auto includedDetectorList =
       builder.value("detectors", std::vector<int>{1, 2, 3, 10});
   const std::set<int> includedDetectors(
       includedDetectorList.begin(), includedDetectorList.end());

   std::ostringstream detectorList;
   for (size_t index = 0; index < includedDetectorList.size(); ++index) {
      if (index) detectorList << ',';
      detectorList << includedDetectorList[index];
   }
   std::cout << "[VMM3aConverter] input=" << m_inputPath
             << ", tree=" << treeName << ", time=" << timeBranch << '\n'
             << "[VMM3aConverter] geometry=" << m_geometryPath
             << ", mappings=" << geometry.size()
             << ", detectors=[" << detectorList.str() << "]\n"
             << "[VMM3aConverter] trigger det=" << triggerDetector
             << ", channels=" << triggerChannels[0] << '/'
             << triggerChannels[1] << ", adc>" << triggerADC
             << ", coincidence=" << coincidenceWindow
             << " ns, dead-time=" << deadTime << " ns\n"
             << "[VMM3aConverter] event window offset=" << windowOffset
             << " ns, before/after=" << windowBefore << '/'
             << windowAfter << " ns, hit adc>=" << hitADCThreshold << '\n';

   if (coincidenceWindow < 0.0 || deadTime < 0.0 ||
       windowBefore < 0.0 || windowAfter < 0.0) {
      std::cerr << "[VMM3aConverter] time windows must be non-negative\n";
      return false;
   }
   if (timeBranch != "time") {
      std::cerr << "[VMM3aConverter] warning: event building with '"
                << timeBranch
                << "'; VMM3a global coincidence normally requires 'time'\n";
   }

   UInt_t rawHitID = 0;
   UChar_t detector = 0;
   UChar_t plane = 0;
   UChar_t fec = 0;
   UChar_t chip = 0;
   UChar_t channel = 0;
   UShort_t adcValue = 0;
   Double_t rawTimeNs = 0.0;
   source->SetBranchAddress("id", &rawHitID);
   source->SetBranchAddress("det", &detector);
   source->SetBranchAddress("plane", &plane);
   source->SetBranchAddress("fec", &fec);
   source->SetBranchAddress("vmm", &chip);
   source->SetBranchAddress("ch", &channel);
   source->SetBranchAddress("adc", &adcValue);
   source->SetBranchAddress(timeBranch.c_str(), &rawTimeNs);

   std::array<std::vector<TriggerCandidate>, 2> triggerCandidates;
   std::vector<VMMHit> hits;
   size_t unmappedHits = 0;
   size_t relabeledPlaneHits = 0;
   std::map<std::pair<int, int>, size_t> mappedHitCounts;
   const Long64_t entries = source->GetEntries();
   hits.reserve(static_cast<size_t>(entries / 2));

   for (Long64_t entryIndex = 0; entryIndex < entries; ++entryIndex) {
      source->GetEntry(entryIndex);

      if (static_cast<int>(detector) == triggerDetector &&
          adcValue > triggerADC) {
         for (size_t slot = 0; slot < triggerChannels.size(); ++slot) {
            if (static_cast<int>(channel) == triggerChannels[slot]) {
               triggerCandidates[slot].push_back(
                   {rawHitID, adcValue, rawTimeNs,
                    rawTimeNs + triggerOffsets[slot]});
            }
         }
      }

      if (!includedDetectors.count(static_cast<int>(detector)) ||
          adcValue < hitADCThreshold) {
         continue;
      }

      const auto mapping = geometry.find(
          {static_cast<int>(detector), static_cast<int>(fec),
           static_cast<int>(chip)});
      if (mapping == geometry.end() || channel >= 64) {
         ++unmappedHits;
         continue;
      }
      const int mappedID0 = mapping->second.id0[channel];
      const int mappedID1 = mapping->second.id1[channel];
      if (mappedID0 < 0) continue;
      if (static_cast<int>(plane) != mapping->second.plane)
         ++relabeledPlaneHits;
      ++mappedHitCounts[
          {mapping->second.detector, mapping->second.plane}];

      hits.push_back(
          {rawHitID,
           mapping->second.detector,
           mapping->second.plane,
           mappedID0,
           mappedID1,
           static_cast<int>(chip),
           static_cast<int>(channel),
           adcValue,
           rawTimeNs});
   }

   const auto byTimeAndID = [](const auto& left, const auto& right) {
      return std::tie(left.correctedTimeNs, left.rawHitID) <
             std::tie(right.correctedTimeNs, right.rawHitID);
   };
   for (auto& candidates : triggerCandidates)
      std::sort(candidates.begin(), candidates.end(), byTimeAndID);
   std::sort(hits.begin(), hits.end(),
             [](const VMMHit& left, const VMMHit& right) {
                return std::tie(left.timeNs, left.rawHitID) <
                       std::tie(right.timeNs, right.rawHitID);
             });

   std::vector<Trigger> triggers;
   std::array<size_t, 2> candidateIndex{};
   double previousTrigger = -std::numeric_limits<double>::infinity();
   while (candidateIndex[0] < triggerCandidates[0].size() &&
          candidateIndex[1] < triggerCandidates[1].size()) {
      const auto& first =
          triggerCandidates[0][candidateIndex[0]];
      const auto& second =
          triggerCandidates[1][candidateIndex[1]];
      const double difference =
          first.correctedTimeNs - second.correctedTimeNs;
      if (std::abs(difference) <= coincidenceWindow) {
         const double triggerTime =
             0.5 * (first.correctedTimeNs + second.correctedTimeNs);
         if (triggerTime - previousTrigger >= deadTime) {
            triggers.push_back(
                {static_cast<unsigned long long>(triggers.size()),
                 triggerTime});
            previousTrigger = triggerTime;
         }
         ++candidateIndex[0];
         ++candidateIndex[1];
      } else if (difference < 0.0) {
         ++candidateIndex[0];
      } else {
         ++candidateIndex[1];
      }
   }

   const fs::path outputPathObject(outputPath);
   if (outputPathObject.has_parent_path())
      fs::create_directories(outputPathObject.parent_path());
   TFile output(outputPath.c_str(), "RECREATE");
   if (output.IsZombie()) {
      std::cerr << "[VMM3aConverter] cannot create " << outputPath << '\n';
      return false;
   }

   ULong64_t evt = 0;
   Double_t eventTimeNs = 0.0;
   std::vector<int> detectors;
   std::vector<int> planes;
   std::vector<int> id0s;
   std::vector<int> id1s;
   std::vector<int> chips;
   std::vector<int> channels;
   std::vector<std::vector<short>> waveforms;
   std::vector<unsigned short> adcs;
   std::vector<unsigned int> rawHitIDs;
   std::vector<double> hitTimesNs;

   TTree eventTree("Events", "BeamAnalysis schema v3 events");
   eventTree.Branch("evt", &evt);
   eventTree.Branch("time", &eventTimeNs);
   eventTree.Branch("det", &detectors);
   eventTree.Branch("plane", &planes);
   eventTree.Branch("id0", &id0s);
   eventTree.Branch("id1", &id1s);
   eventTree.Branch("chip", &chips);
   eventTree.Branch("channel", &channels);
   eventTree.Branch("waveform", &waveforms);
   eventTree.Branch("adc", &adcs);
   eventTree.Branch("rawhitid", &rawHitIDs);
   eventTree.Branch("hit_time_ns", &hitTimesNs);

   std::map<std::tuple<int, int, int, int, int, int>, ChannelRow>
       observedChannels;
   unsigned long long eventHitCount = 0;
   for (const auto& trigger : triggers) {
      evt = trigger.id;
      eventTimeNs = trigger.timeNs;
      detectors.clear();
      planes.clear();
      id0s.clear();
      id1s.clear();
      chips.clear();
      channels.clear();
      waveforms.clear();
      adcs.clear();
      rawHitIDs.clear();
      hitTimesNs.clear();

      const double start =
          trigger.timeNs + windowOffset - windowBefore;
      const double stop =
          trigger.timeNs + windowOffset + windowAfter;
      const auto begin = std::lower_bound(
          hits.begin(), hits.end(), start,
          [](const VMMHit& hit, double value) {
             return hit.timeNs < value;
          });
      const auto end = std::lower_bound(
          begin, hits.end(), stop,
          [](const VMMHit& hit, double value) {
             return hit.timeNs < value;
          });

      const size_t count =
          static_cast<size_t>(std::distance(begin, end));
      detectors.reserve(count);
      planes.reserve(count);
      id0s.reserve(count);
      id1s.reserve(count);
      chips.reserve(count);
      channels.reserve(count);
      waveforms.reserve(count);
      adcs.reserve(count);
      rawHitIDs.reserve(count);
      hitTimesNs.reserve(count);

      for (auto hit = begin; hit != end; ++hit) {
         detectors.push_back(hit->detector);
         planes.push_back(hit->plane);
         id0s.push_back(hit->id0);
         id1s.push_back(hit->id1);
         chips.push_back(hit->chip);
         channels.push_back(hit->channel);
         waveforms.emplace_back();
         adcs.push_back(hit->adc);
         rawHitIDs.push_back(hit->rawHitID);
         hitTimesNs.push_back(hit->timeNs - trigger.timeNs);

         ChannelRow row{hit->chip, hit->channel, hit->detector,
                        hit->plane, hit->id0, hit->id1};
         observedChannels.emplace(row.Key(), row);
      }
      eventHitCount += count;
      eventTree.Fill();
   }

   int mapChip = -1;
   int mapChannel = -1;
   int mapDetector = -1;
   int mapPlane = -1;
   int mapID0 = -1;
   int mapID1 = -1;
   double pedestal = 0.0;
   double noiseSigma = -1.0;
   double gain = 1.0;
   int polarity = 1;
   TTree channelTree("Channels", "BeamAnalysis schema v3 channels");
   channelTree.Branch("chip", &mapChip);
   channelTree.Branch("channel", &mapChannel);
   channelTree.Branch("det", &mapDetector);
   channelTree.Branch("plane", &mapPlane);
   channelTree.Branch("id0", &mapID0);
   channelTree.Branch("id1", &mapID1);
   channelTree.Branch("pedestal", &pedestal);
   channelTree.Branch("noise_sigma", &noiseSigma);
   channelTree.Branch("gain", &gain);
   channelTree.Branch("polarity", &polarity);
   for (const auto& [key, row] : observedChannels) {
      (void)key;
      mapChip = row.chip;
      mapChannel = row.channel;
      mapDetector = row.detector;
      mapPlane = row.plane;
      mapID0 = row.id0;
      mapID1 = row.id1;
      channelTree.Fill();
   }

   int schemaVersion = 3;
   TTree metadataTree("Metadata",
                      "BeamAnalysis canonical format metadata");
   metadataTree.Branch("schema_version", &schemaVersion);
   metadataTree.Fill();

   output.cd();
   eventTree.Write();
   channelTree.Write();
   metadataTree.Write();
   output.Close();
   input.Close();

   const double meanMultiplicity =
       triggers.empty()
           ? 0.0
           : static_cast<double>(eventHitCount) / triggers.size();
   std::cout << "[VMM3aConverter] input hits=" << entries
             << ", trigger candidates=" << triggerCandidates[0].size()
             << '/' << triggerCandidates[1].size()
             << ", events=" << triggers.size()
             << ", event hits=" << eventHitCount
             << ", mean multiplicity=" << std::fixed
             << std::setprecision(2) << meanMultiplicity
             << ", active channels=" << observedChannels.size() << '\n';
   std::ostringstream mappedSummary;
   for (const auto& [key, count] : mappedHitCounts) {
      if (mappedSummary.tellp() > 0) mappedSummary << ", ";
      mappedSummary << "det" << key.first << ".plane" << key.second
                    << '=' << count;
   }
   std::cout << "[VMM3aConverter] mapped hits: "
             << mappedSummary.str();
   if (relabeledPlaneHits)
      std::cout << " (plane relabeled=" << relabeledPlaneHits << ')';
   std::cout << '\n'
             << "[VMM3aConverter] output=" << outputPath
             << ", schema=v3, Events=" << triggers.size()
             << ", Channels=" << observedChannels.size() << '\n';
   if (unmappedHits > 0) {
      std::cerr << "[VMM3aConverter] skipped " << unmappedHits
                << " included-detector hits without geometry mapping\n";
   }
   return true;
}
