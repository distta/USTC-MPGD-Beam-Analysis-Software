#include "Script/Base/RawDataParser.h"
#include "Detector/DetectorFactory.h"

#include <TDirectory.h>
#include <TH1D.h>
#include <TH2D.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <utility>

RawDataParser::RawDataParser(const std::string& rawFile)
    : m_rawFile(rawFile) {
}

RawDataParser::~RawDataParser() {
   if (m_file) {
      m_file->Close();
      delete m_file;
      m_file = nullptr;
   }
}

bool RawDataParser::Initialize() {
   // 打开ROOT文件
   m_file = TFile::Open(m_rawFile.c_str(), "READ");
   if (!m_file || m_file->IsZombie()) {
      std::cerr << "[RawDataParser] Failed to open file: "
                << m_rawFile << std::endl;
      return false;
   }

   m_tree = static_cast<TTree*>(m_file->Get("Events"));
   if (!m_tree) {
      std::cerr << "[RawDataParser] Events tree not found in "
                << m_rawFile
                << "; convert APV25SRS/BT input to the standard ROOT format first"
                << std::endl;
      return false;
   }

   TTree* metadataTree = static_cast<TTree*>(m_file->Get("Metadata"));
   if (metadataTree && metadataTree->GetBranch("schema_version")) {
      metadataTree->SetBranchAddress("schema_version", &m_schemaVersion);
      metadataTree->GetEntry(0);
   } else {
      m_schemaVersion = 0;
   }
   if (m_schemaVersion != 2 && m_schemaVersion != 3) {
      std::cerr << "[RawDataParser] Unsupported schema version: "
                << m_schemaVersion
                << "; supported versions are 2 and 3"
                << std::endl;
      return false;
   }

   const std::vector<std::string> required =
       m_schemaVersion == 3
           ? std::vector<std::string>{
                 "evt", "time", "det", "plane", "id0", "id1",
                 "chip", "channel", "waveform", "adc", "rawhitid",
                 "hit_time_ns"}
           : std::vector<std::string>{
                 "event_id", "timestamp", "detector_id", "plane_type",
                 "id0", "id1", "waveform"};

   for (const auto& branch : required) {
      if (!m_tree->GetBranch(branch.c_str())) {
         std::cerr << "[RawDataParser] Missing branch: "
                   << branch << std::endl;
         return false;
      }
   }

   if (m_schemaVersion == 3) {
      m_tree->SetBranchAddress("evt", &m_eventID);
      m_tree->SetBranchAddress("time", &m_eventTimeNs);
      m_tree->SetBranchAddress("det", &m_detectorIDs);
      m_tree->SetBranchAddress("plane", &m_planeTypes);
      m_tree->SetBranchAddress("id0", &m_id0s);
      m_tree->SetBranchAddress("id1", &m_id1s);
      m_tree->SetBranchAddress("chip", &m_chipIDs);
      m_tree->SetBranchAddress("channel", &m_channelIDs);
      m_tree->SetBranchAddress("waveform", &m_waveforms);
      m_tree->SetBranchAddress("adc", &m_adcs);
      m_tree->SetBranchAddress("rawhitid", &m_rawHitIDs);
      m_tree->SetBranchAddress("hit_time_ns", &m_hitTimesNs);
   } else {
      m_tree->SetBranchAddress("event_id", &m_eventID);
      m_tree->SetBranchAddress("timestamp", &m_timestamp);
      m_tree->SetBranchAddress("detector_id", &m_detectorIDs);
      m_tree->SetBranchAddress("plane_type", &m_planeTypes);
      m_tree->SetBranchAddress("id0", &m_id0s);
      m_tree->SetBranchAddress("id1", &m_id1s);
      m_tree->SetBranchAddress("waveform", &m_waveforms);
   }

   m_numOfEvents = m_tree->GetEntries();

   if (!LoadCanonicalChannelData()) {
      return false;
   }

   return true;
}

bool RawDataParser::LoadCanonicalChannelData() {
   TTree* channelTree = static_cast<TTree*>(m_file->Get("Channels"));
   if (!channelTree) {
      std::cerr << "[RawDataParser] Channels tree not found; "
                << "noise sigma is unavailable"
                << std::endl;
      return false;
   }

   int detectorID = 0;
   int planeType = 0;
   int id0 = 0;
   int id1 = -1;

   double pedestal = 0.0;
   double noiseSigma = -1.0;
   double gain = 1.0;

   int polarity = 1;

   const char* detectorBranch =
       m_schemaVersion == 3 ? "det" : "detector_id";
   const char* planeBranch =
       m_schemaVersion == 3 ? "plane" : "plane_type";
   const char* required[] = {
       detectorBranch, planeBranch, "id0", "id1"};
   for (const char* branch : required) {
      if (!channelTree->GetBranch(branch)) {
         std::cerr << "[RawDataParser] Channels tree is missing branch: "
                   << branch << std::endl;
         return false;
      }
   }

   channelTree->SetBranchAddress(detectorBranch, &detectorID);
   channelTree->SetBranchAddress(planeBranch, &planeType);
   channelTree->SetBranchAddress("id0", &id0);
   channelTree->SetBranchAddress("id1", &id1);

   if (channelTree->GetBranch("pedestal")) {
      channelTree->SetBranchAddress("pedestal", &pedestal);
   }

   if (channelTree->GetBranch("noise_sigma")) {
      channelTree->SetBranchAddress("noise_sigma", &noiseSigma);
   }

   if (channelTree->GetBranch("gain")) {
      channelTree->SetBranchAddress("gain", &gain);
   }

   if (channelTree->GetBranch("polarity")) {
      channelTree->SetBranchAddress("polarity", &polarity);
   }

   const auto entries = channelTree->GetEntries();

   for (Long64_t i = 0; i < entries; ++i) {
      channelTree->GetEntry(i);
      const ChannelCoordinate coordinate{id0, id1};

      if (noiseSigma >= 0) {
         m_pedSigmaMap[detectorID][planeType][coordinate] = noiseSigma;
      }

      m_calibrationMap[detectorID][planeType][coordinate] = {
          pedestal,
          gain,
          polarity};
   }

   return true;
}

double RawDataParser::GetSigma(int detID, int type, int id0, int id1) const {
   auto it1 = m_pedSigmaMap.find(detID);
   if (it1 == m_pedSigmaMap.end()) return -1;

   auto it2 = it1->second.find(type);
   if (it2 == it1->second.end()) return -1;

   auto it3 = it2->second.find({id0, id1});
   if (it3 == it2->second.end()) return -1;

   return it3->second;
}

bool RawDataParser::WriteDebugRoot(const std::string& outputFile) {
   if (!m_tree) {
      std::cerr << "[RawDataParser] ERROR: TTree not initialized\n";
      return false;
   }

   TFile debugFile(outputFile.c_str(), "RECREATE");
   if (debugFile.IsZombie()) {
      std::cerr << "[RawDataParser] ERROR: cannot create debug file: "
                << outputFile << '\n';
      return false;
   }

   const Long64_t nEvents = m_tree->GetEntries();

   std::cout << "[RawDataParser] Writing debug ROOT: "
             << outputFile << std::endl;

   for (Long64_t entry = 0; entry < nEvents; ++entry) {
      m_tree->GetEntry(entry);

      if (!m_detectorIDs || !m_planeTypes || !m_id0s || !m_id1s || !m_waveforms) {
         continue;
      }

      const size_t nHits = m_detectorIDs->size();

      if (m_planeTypes->size() != nHits ||
          m_id0s->size() != nHits ||
          m_id1s->size() != nHits ||
          m_waveforms->size() != nHits) {
         std::cerr << "[RawDataParser] Inconsistent vectors at entry "
                   << entry << '\n';
         continue;
      }

      std::ostringstream dirName;
      dirName << "event_" << entry;

      TDirectory* eventDir = debugFile.mkdir(dirName.str().c_str());
      if (!eventDir) {
         eventDir = debugFile.GetDirectory(dirName.str().c_str());
      }

      if (!eventDir) {
         continue;
      }

      eventDir->cd();

      std::map<std::pair<int, int>, TH1D*> hMaxMap;
      std::map<std::pair<int, int>, TH2D*> hWaveMap;

      for (size_t i = 0; i < nHits; ++i) {
         const int detID = (*m_detectorIDs)[i];
         const int planeType = (*m_planeTypes)[i];
         const int id0 = (*m_id0s)[i];
         const int id1 = (*m_id1s)[i];

         if (detID < 0) continue;
         if (planeType < 0) continue;
         if (id0 < 0) continue;

         const auto& waveform = (*m_waveforms)[i];
         if (waveform.empty()) continue;

         const auto key = std::make_pair(detID, planeType);

         if (hMaxMap.find(key) == hMaxMap.end()) {
            std::ostringstream hMaxName;
            hMaxName << "h_maxA_event_" << entry
                     << "_det_" << detID
                     << "_plane_" << planeType;

            std::ostringstream hMaxTitle;
            hMaxTitle << "Event " << entry
                      << ", detID " << detID
                      << ", plane " << planeType
                      << ";id0;maxA";

            hMaxMap[key] = new TH1D(
                hMaxName.str().c_str(),
                hMaxTitle.str().c_str(),
                256,
                -0.5,
                255.5);

            const int nSamples = static_cast<int>(waveform.size());

            std::ostringstream hWaveName;
            hWaveName << "h_waveform_event_" << entry
                      << "_det_" << detID
                      << "_plane_" << planeType;

            std::ostringstream hWaveTitle;
            hWaveTitle << "Event " << entry
                       << ", detID " << detID
                       << ", plane " << planeType
                       << ";sample;id0;ADC";

            hWaveMap[key] = new TH2D(
                hWaveName.str().c_str(),
                hWaveTitle.str().c_str(),
                nSamples,
                -0.5,
                nSamples - 0.5,
                256,
                -0.5,
                255.5);
         }

         double maxA = -std::numeric_limits<double>::max();

         const auto detIt = m_calibrationMap.find(detID);

         for (size_t s = 0; s < waveform.size(); ++s) {
            double value = static_cast<double>(waveform[s]);

            if (detIt != m_calibrationMap.end()) {
               const auto planeIt = detIt->second.find(planeType);

               if (planeIt != detIt->second.end()) {
                  const auto stripIt = planeIt->second.find({id0, id1});

                  if (stripIt != planeIt->second.end()) {
                     const auto& calibration = stripIt->second;

                     value = (value - calibration.pedestal) * calibration.gain * calibration.polarity;
                  }
               }
            }

            if (value > maxA) {
               maxA = value;
            }

            TH2D* hWave = hWaveMap[key];

            const int xbin = hWave->GetXaxis()->FindBin(
                static_cast<double>(s));
            const int ybin = hWave->GetYaxis()->FindBin(id0);

            hWave->SetBinContent(xbin, ybin, value);
         }

         TH1D* hMax = hMaxMap[key];

         const int bin = hMax->GetXaxis()->FindBin(id0);
         const double oldValue = hMax->GetBinContent(bin);

         if (maxA > oldValue) {
            hMax->SetBinContent(bin, maxA);
         }
      }

      for (auto& item : hMaxMap) {
         item.second->Write();
         delete item.second;
      }

      for (auto& item : hWaveMap) {
         item.second->Write();
         delete item.second;
      }

      debugFile.cd();

      if ((entry + 1) % 1000 == 0) {
         std::cout << "[RawDataParser] debug events written: "
                   << (entry + 1) << " / " << nEvents << std::endl;
      }
   }

   debugFile.Close();

   std::cout << "[RawDataParser] Debug ROOT written: "
             << outputFile << std::endl;

   return true;
}

std::unordered_map<int, std::vector<RawData>>
    RawDataParser::LoadEvent(int eventID) {
   std::unordered_map<int, std::vector<RawData>> result;

   if (!m_tree) {
      std::cerr << "[RawDataParser] ERROR: TTree not initialized\n";
      return result;
   }

   if (eventID < 0 || eventID >= GetTotalEvents()) {
      std::cerr << "[RawDataParser] ERROR: Invalid event index "
                << eventID << "\n";
      return result;
   }

   // ---- Load TTree Entry ----
   m_tree->GetEntry(eventID);

   auto& factory = DetectorFactory::GetInstance();
   const auto& detectors = factory.GetAllDetectors();

   if (!m_detectorIDs || !m_planeTypes || !m_id0s || !m_id1s ||
       !m_waveforms) {
      return result;
   }

   const size_t nHits = m_detectorIDs->size();

   if (m_planeTypes->size() != nHits ||
       m_id0s->size() != nHits ||
       m_id1s->size() != nHits ||
       m_waveforms->size() != nHits) {
      std::cerr << "[RawDataParser] Inconsistent vectors at entry "
                << eventID << '\n';
      return result;
   }
   if (m_schemaVersion == 3 &&
       (!m_chipIDs || !m_channelIDs || !m_adcs || !m_rawHitIDs ||
        !m_hitTimesNs ||
        m_chipIDs->size() != nHits ||
        m_channelIDs->size() != nHits ||
        m_adcs->size() != nHits ||
        m_rawHitIDs->size() != nHits ||
        m_hitTimesNs->size() != nHits)) {
      std::cerr << "[RawDataParser] Inconsistent schema v3 vectors at entry "
                << eventID << '\n';
      return result;
   }

   result.reserve(16);

   for (size_t i = 0; i < nHits; ++i) {
      const int detID = (*m_detectorIDs)[i];

      if (detectors.find(detID) == detectors.end()) {
         continue;
      }

      const int planeType = (*m_planeTypes)[i];
      const int id0 = (*m_id0s)[i];
      const int id1Value = (*m_id1s)[i];
      if (id0 < 0) continue;

      if (m_schemaVersion == 3) {
         double calibratedADC = static_cast<double>((*m_adcs)[i]);
         const auto detIt = m_calibrationMap.find(detID);
         if (detIt != m_calibrationMap.end()) {
            const auto planeIt = detIt->second.find(planeType);
            if (planeIt != detIt->second.end()) {
               const auto channelIt =
                   planeIt->second.find({id0, id1Value});
               if (channelIt != planeIt->second.end()) {
                  const auto& calibration = channelIt->second;
                  calibratedADC =
                      (calibratedADC - calibration.pedestal) *
                      calibration.gain * calibration.polarity;
               }
            }
         }

         RawData decoded{id0, id1Value, planeType, {}};
         decoded.chip = (*m_chipIDs)[i];
         decoded.channel = (*m_channelIDs)[i];
         decoded.rawHitID = (*m_rawHitIDs)[i];
         decoded.adcValue = calibratedADC;
         decoded.hitTimeNs = (*m_hitTimesNs)[i];
         result[detID].push_back(std::move(decoded));
         continue;
      }

      std::vector<short> calibrated = (*m_waveforms)[i];
      double gain = 1.0;

      const auto detIt = m_calibrationMap.find(detID);
      if (detIt != m_calibrationMap.end()) {
         const auto planeIt = detIt->second.find(planeType);
         if (planeIt != detIt->second.end()) {
            const auto stripIt = planeIt->second.find({id0, id1Value});
            if (stripIt != planeIt->second.end()) {
               const auto& calibration = stripIt->second;
               gain = calibration.gain;

               for (auto& sample : calibrated) {
                  const double value =
                      (sample - calibration.pedestal) * calibration.gain * calibration.polarity;

                  sample = static_cast<short>(std::clamp(
                      std::lround(value),
                      static_cast<long>(
                          std::numeric_limits<short>::min()),
                      static_cast<long>(
                          std::numeric_limits<short>::max())));
               }
            }
         }
      }

      const double sigma = GetSigma(detID, planeType, id0, id1Value);
      if (sigma >= 0) {
         const auto maxIt = std::max_element(
             calibrated.begin(),
             calibrated.end());

         if (maxIt == calibrated.end()) {
            continue;
         }

         const double threshold = 3.0 * sigma * std::abs(gain);
         if (static_cast<double>(*maxIt) <= threshold) {
            continue;
         }
      }

      result[detID].push_back(
          RawData{id0, id1Value, planeType, std::move(calibrated)});
   }

   return result;
}
