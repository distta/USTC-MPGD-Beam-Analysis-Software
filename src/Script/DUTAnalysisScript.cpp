#include "Script/DUTAnalysisScript.h"
#include "Algorithm/AnalysisUtils.h"
#include "Detector/DetectorFactory.h"
#include "Event/DataModel.h"
#include "Event/DetectorFrame.h"
#include "Script/Base/RawDataParser.h"
#include "Script/Base/ScriptFactory.h"

#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include <TCanvas.h>
#include <TDirectory.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TMath.h>
#include <TPad.h>
#include <TStyle.h>
#include <TTree.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <unordered_set>

using namespace std;

void RunDUTAlign(const std::vector<Event>& events, std::shared_ptr<Detector> detector, int detID);
LocalHit CalcuDutResidual(std::shared_ptr<Detector> detector, const std::vector<Cluster>& clusters, const TVector3& predL, double& residualX, double& residualY);
double DUTChi2Objective(const double* par, const std::vector<Event>& events, std::shared_ptr<Detector> detector, int detID);
struct ResidualAnalysisResult {
   double sigma68 = 0;
   double sigmaFit = 0;
   double sigmaFitError = 0;
   double mean = 0;
   double efficiency = 0;  // valid reconstruction / total events
   double eff3S = 0;
   double eff5S = 0;
   double eff1mm = 0;  // |residual| < 1 mm
};

struct DUTAlignmentQAPoint {
   double predX = 0;
   double predY = 0;
   double hitX = 0;
   double hitY = 0;
   double resX = 0;
   double resY = 0;
   bool validX = false;
   bool validY = false;
};

ResidualAnalysisResult AnalyzeResidualSequence(const std::vector<double>& residuals, TFile* outputFile, const std::string& histName);
void WriteDUTAlignmentQA(TFile* outputFile, int detID,
                         const std::vector<DUTAlignmentQAPoint>& points,
                         double predXMin, double predXMax,
                         double predYMin, double predYMax);


void WriteDUTAlignmentQA(TFile* outputFile, int detID,
                         const std::vector<DUTAlignmentQAPoint>& points,
                         double predXMin, double predXMax,
                         double predYMin, double predYMax) {
   if (!outputFile || points.empty()) return;

   std::vector<double> absResX, absResY;
   double hitXMin = predXMin, hitXMax = predXMax;
   double hitYMin = predYMin, hitYMax = predYMax;
   for (const auto& point : points) {
      if (point.validX) {
         absResX.push_back(std::abs(point.resX));
         hitXMin = std::min(hitXMin, point.hitX);
         hitXMax = std::max(hitXMax, point.hitX);
      }
      if (point.validY) {
         absResY.push_back(std::abs(point.resY));
         hitYMin = std::min(hitYMin, point.hitY);
         hitYMax = std::max(hitYMax, point.hitY);
      }
   }

   auto residualRange = [](std::vector<double> values) {
      if (values.empty()) return 1.0;
      const size_t index = std::min(values.size() - 1,
                                    static_cast<size_t>(0.99 * values.size()));
      std::nth_element(values.begin(), values.begin() + index, values.end());
      return std::clamp(1.2 * values[index], 1.0, 100.0);
   };
   const double resXRange = residualRange(std::move(absResX));
   const double resYRange = residualRange(std::move(absResY));

   const double hitXPadding = std::max(0.5, 0.02 * (hitXMax - hitXMin));
   const double hitYPadding = std::max(0.5, 0.02 * (hitYMax - hitYMin));
   hitXMin -= hitXPadding;
   hitXMax += hitXPadding;
   hitYMin -= hitYPadding;
   hitYMax += hitYPadding;

   auto* qaDirectory = outputFile->GetDirectory("AlignmentQA");
   if (!qaDirectory) qaDirectory = outputFile->mkdir("AlignmentQA");
   const std::string detectorDirectoryName = "DUT_" + std::to_string(detID);
   auto* detectorDirectory = qaDirectory->GetDirectory(detectorDirectoryName.c_str());
   if (!detectorDirectory) detectorDirectory = qaDirectory->mkdir(detectorDirectoryName.c_str());
   TDirectory::TContext context(detectorDirectory);

   TH2D hitXVsPredX(Form("hHitXVsPredX_DUT%d", detID),
                    Form("DUT %d measured X vs track intercept X;Track intercept X [mm];Measured hit X [mm]", detID),
                    160, predXMin, predXMax, 160, hitXMin, hitXMax);
   TH2D resXVsPredX(Form("hResXVsPredX_DUT%d", detID),
                    Form("DUT %d residual X vs track intercept X;Track intercept X [mm];Residual X [mm]", detID),
                    160, predXMin, predXMax, 200, -resXRange, resXRange);
   TH2D resXVsHitX(Form("hResXVsHitX_DUT%d", detID),
                   Form("DUT %d residual X vs measured X;Measured hit X [mm];Residual X [mm]", detID),
                   160, hitXMin, hitXMax, 200, -resXRange, resXRange);
   TH2D hitYVsPredY(Form("hHitYVsPredY_DUT%d", detID),
                    Form("DUT %d measured Y vs track intercept Y;Track intercept Y [mm];Measured hit Y [mm]", detID),
                    160, predYMin, predYMax, 160, hitYMin, hitYMax);
   TH2D resYVsPredY(Form("hResYVsPredY_DUT%d", detID),
                    Form("DUT %d residual Y vs track intercept Y;Track intercept Y [mm];Residual Y [mm]", detID),
                    160, predYMin, predYMax, 200, -resYRange, resYRange);
   TH2D resYVsHitY(Form("hResYVsHitY_DUT%d", detID),
                   Form("DUT %d residual Y vs measured Y;Measured hit Y [mm];Residual Y [mm]", detID),
                   160, hitYMin, hitYMax, 200, -resYRange, resYRange);

   for (const auto& point : points) {
      if (point.validX) {
         hitXVsPredX.Fill(point.predX, point.hitX);
         resXVsPredX.Fill(point.predX, point.resX);
         resXVsHitX.Fill(point.hitX, point.resX);
      }
      if (point.validY) {
         hitYVsPredY.Fill(point.predY, point.hitY);
         resYVsPredY.Fill(point.predY, point.resY);
         resYVsHitY.Fill(point.hitY, point.resY);
      }
   }

   for (TH2D* histogram : {&hitXVsPredX, &resXVsPredX, &resXVsHitX,
                           &hitYVsPredY, &resYVsPredY, &resYVsHitY}) {
      histogram->SetStats(false);
      histogram->Write();
   }

   TCanvas canvasX(Form("cAlignmentX_DUT%d", detID),
                   Form("DUT %d X alignment QA", detID), 1800, 600);
   canvasX.Divide(3, 1);
   canvasX.cd(1);
   hitXVsPredX.Draw("COLZ");
   canvasX.cd(2);
   resXVsPredX.Draw("COLZ");
   canvasX.cd(3);
   resXVsHitX.Draw("COLZ");
   canvasX.Write();

   TCanvas canvasY(Form("cAlignmentY_DUT%d", detID),
                   Form("DUT %d Y alignment QA", detID), 1800, 600);
   canvasY.Divide(3, 1);
   canvasY.cd(1);
   hitYVsPredY.Draw("COLZ");
   canvasY.cd(2);
   resYVsPredY.Draw("COLZ");
   canvasY.cd(3);
   resYVsHitY.Draw("COLZ");
   canvasY.Write();
}

void DUTAnalysisScript::LoadConfig(const json& config) {
   m_runAlignment = config.value("runAlignment", true);
   m_progressInterval = config.value("progressInterval", 1000);
   m_maxEvents = config.value("maxEvents", -1);

   const json effCfg = config.value("efficiencyMap", json::object());
   m_effXMin = effCfg.value("xMin", m_effXMin);
   m_effXMax = effCfg.value("xMax", m_effXMax);
   m_effYMin = effCfg.value("yMin", m_effYMin);
   m_effYMax = effCfg.value("yMax", m_effYMax);
   m_effXBins = effCfg.value("xBins", m_effXBins);
   m_effYBins = effCfg.value("yBins", m_effYBins);
   m_effMinEntriesPerBin = effCfg.value(
       "minEntriesPerBin", m_effMinEntriesPerBin);
   m_effExcludedXBins = effCfg.value("excludeXBins", std::vector<int>{});
   m_effExcludedYBins = effCfg.value("excludeYBins", std::vector<int>{});

   if (m_effXBins <= 0) m_effXBins = 1;
   if (m_effYBins <= 0) m_effYBins = 1;
   if (m_effMinEntriesPerBin <= 0) m_effMinEntriesPerBin = 1;
   if (m_effXMax < m_effXMin) std::swap(m_effXMax, m_effXMin);
   if (m_effYMax < m_effYMin) std::swap(m_effYMax, m_effYMin);
   if (m_effXMax == m_effXMin) m_effXMax = m_effXMin + 1.0;
   if (m_effYMax == m_effYMin) m_effYMax = m_effYMin + 1.0;

   auto cleanExcludedBins = [](std::vector<int>& bins, int maxBin) {
      bins.erase(std::remove_if(bins.begin(), bins.end(),
                                [maxBin](int bin) { return bin < 1 || bin > maxBin; }),
                 bins.end());
      std::sort(bins.begin(), bins.end());
      bins.erase(std::unique(bins.begin(), bins.end()), bins.end());
   };
   cleanExcludedBins(m_effExcludedXBins, m_effXBins);
   cleanExcludedBins(m_effExcludedYBins, m_effYBins);

   m_efficiencyConfig.xBins = m_effXBins;
   m_efficiencyConfig.yBins = m_effYBins;
   m_efficiencyConfig.xMin = m_effXMin;
   m_efficiencyConfig.xMax = m_effXMax;
   m_efficiencyConfig.yMin = m_effYMin;
   m_efficiencyConfig.yMax = m_effYMax;
   m_efficiencyConfig.minEntriesPerBin = m_effMinEntriesPerBin;
   m_efficiencyConfig.excludedXBins = m_effExcludedXBins;
   m_efficiencyConfig.excludedYBins = m_effExcludedYBins;
   m_efficiencyConfig.margin = effCfg.value("margin", m_efficiencyConfig.margin);
   const json marginScan = effCfg.value("marginScan", json::object());
   m_efficiencyConfig.envelopeScanMin = marginScan.value(
       "min", m_efficiencyConfig.envelopeScanMin);
   m_efficiencyConfig.envelopeScanMax = marginScan.value(
       "max", m_efficiencyConfig.envelopeScanMax);
   m_efficiencyConfig.envelopeScanStep = marginScan.value(
       "step", m_efficiencyConfig.envelopeScanStep);

   const json fake = effCfg.value("fake", json::object());
   m_efficiencyConfig.enableFakeEfficiency = fake.value(
       "enabled", m_efficiencyConfig.enableFakeEfficiency);
   m_efficiencyConfig.fakeSeed = fake.value(
       "seed", m_efficiencyConfig.fakeSeed);
   m_efficiencyConfig.fakePartnersPerEvent = std::max(
       1, fake.value("partnersPerEvent", m_efficiencyConfig.fakePartnersPerEvent));

}

void DUTAnalysisScript::Print() const {
   cout << "DUTAnalysisScript Configuration:" << endl;
   cout << "  Run Alignment: " << (m_runAlignment ? "Yes" : "No") << endl;
   cout << "  Efficiency Region X: [" << m_effXMin << ", " << m_effXMax << "] mm, bins=" << m_effXBins << endl;
   cout << "  Efficiency Region Y: [" << m_effYMin << ", " << m_effYMax << "] mm, bins=" << m_effYBins << endl;
   cout << "  Efficiency Minimum Entries/Bin: "
        << m_effMinEntriesPerBin << endl;
   cout << "  Efficiency Averaging: Event-weighted" << endl;
   cout << "  Cluster Envelope Margin: "
        << m_efficiencyConfig.margin << " mm" << endl;
   cout << "  Fake Efficiency: "
        << (m_efficiencyConfig.enableFakeEfficiency ? "enabled" : "disabled")
        << ", random partners/event="
        << m_efficiencyConfig.fakePartnersPerEvent << endl;
   if (!m_effExcludedXBins.empty() || !m_effExcludedYBins.empty()) {
      cout << "  Efficiency Excluded X bins:";
      for (int bin : m_effExcludedXBins) cout << " " << bin;
      if (m_effExcludedXBins.empty()) cout << " none";
      cout << endl;
      cout << "  Efficiency Excluded Y bins:";
      for (int bin : m_effExcludedYBins) cout << " " << bin;
      if (m_effExcludedYBins.empty()) cout << " none";
      cout << endl;
   }
   cout << "  Progress Interval: " << m_progressInterval << endl;
   cout << "  Max Events: " << (m_maxEvents > 0 ? to_string(m_maxEvents) : "All") << endl;
}

bool DUTAnalysisScript::Execute() {
   auto t0 = chrono::high_resolution_clock::now();

   auto parser = GetParser();
   if (!parser) {
      cerr << "Error: Parser not set!" << endl;
      return false;
   }

   auto& factory = DetectorFactory::GetInstance();

   Print();

   // 加载track信息
   string trackFile = GetOutputDir() + "TrackInfo.root";
   cout << "Loading track info..." << endl;
   cout << "File: " << trackFile << endl;

   TFile* f = TFile::Open(trackFile.c_str(), "READ");
   if (!f || f->IsZombie()) {
      cerr << "Error: Cannot open " << trackFile << endl;
      cerr << "Please run Track Analysis first!" << endl;
      return false;
   }

   TTree* trackTree = (TTree*)f->Get("Tracks");
   if (!trackTree) {
      cerr << "Error: No Tracks tree!" << endl;
      f->Close();
      return false;
   }

   Int_t eventID;
   Track* track = nullptr;
   double sigTime = 0.0;

   // Tracks contains one row per reconstructed track.  Count event IDs using
   // only the lightweight scalar branch first, then analyze exactly-one-track
   // events in the full pass below.
   trackTree->SetBranchStatus("*", false);
   trackTree->SetBranchStatus("eventID", true);
   trackTree->SetBranchAddress("eventID", &eventID);
   const Long64_t totalTrackEntries = trackTree->GetEntries();
   std::unordered_map<int, int> tracksPerEvent;
   tracksPerEvent.reserve(static_cast<size_t>(totalTrackEntries));
   for (Long64_t i = 0; i < totalTrackEntries; ++i) {
      trackTree->GetEntry(i);
      ++tracksPerEvent[eventID];
   }

   Long64_t singleTrackEvents = 0;
   for (const auto& [id, count] : tracksPerEvent) {
      (void)id;
      if (count == 1) ++singleTrackEvents;
   }

   trackTree->SetBranchStatus("*", true);
   trackTree->SetBranchAddress("track", &track);
   if (trackTree->GetBranch("t0")) trackTree->SetBranchAddress("t0", &sigTime);

   cout << "\nProcessing DUT data..." << endl;
   std::vector<Event> events;  // Script本地数据

   int processed = 0;

   cout << "[DUTAnalysis] track entries=" << totalTrackEntries
        << ", unique events: " << tracksPerEvent.size()
        << ", single-track events: " << singleTrackEvents << endl;

   const auto& duts = factory.GetDetectorsByRole(Detector::Role::DUT);
   const auto& trackers = factory.GetDetectorsByRole(Detector::Role::Tracker);
   std::unordered_set<int> strictSingleHitTrackerEvents;

   string dutFile = GetOutputDir() + "DUTInfo.root";
   TFile* fDut = new TFile(dutFile.c_str(), "RECREATE");
   TTree* tDut = new TTree("DUTTree", "DUT data");

   // 加载事件并处理DUT数据
   const Long64_t processingTarget =
       m_maxEvents > 0
           ? std::min<Long64_t>(m_maxEvents, singleTrackEvents)
           : singleTrackEvents;
   events.reserve(static_cast<size_t>(processingTarget));
   strictSingleHitTrackerEvents.reserve(static_cast<size_t>(processingTarget));

   auto printProgress = [&](Long64_t completed) {
      const double percent = processingTarget > 0
                                 ? 100.0 * completed / processingTarget
                                 : 100.0;
      cout << "\r  Processing single-track events: " << completed
           << "/" << processingTarget << " ("
           << fixed << setprecision(1) << percent << "%)" << flush;
   };
   printProgress(0);

   Long64_t lastProgress = 0;
   for (Long64_t i = 0;
        i < totalTrackEntries && processed < processingTarget; ++i) {
      trackTree->GetEntry(i);

      const auto multiplicity = tracksPerEvent.find(eventID);
      if (multiplicity == tracksPerEvent.end() || multiplicity->second != 1) continue;

      auto rawHits = parser->LoadEvent(eventID);
      if (rawHits.empty()) continue;

      if (m_efficiencyConfig.enableFakeEfficiency) {
         bool strictSingleHit = true;
         for (const auto& trackerDetector : trackers) {
            auto trackerFrame = make_shared<DetectorFrame>(*trackerDetector);
            const auto raw = rawHits.find(trackerDetector->GetID());
            if (raw != rawHits.end()) trackerFrame->SetRawData(raw->second);
            trackerFrame->Process();
            if (trackerFrame->LocalHits().size() != 1) {
               strictSingleHit = false;
               break;
            }
         }
         if (strictSingleHit) strictSingleHitTrackerEvents.insert(eventID);
      }

      Event evt{.eventID = int(eventID), .track = *track, .t0 = sigTime};

      for (auto& det : duts) {
         const int id = det->GetID();
         auto detEvt = make_shared<DetectorFrame>(*det);
         detEvt->SetRawData(rawHits[det->GetID()]);
         detEvt->Process(evt.t0);
         evt.detectorFramesMap[id] = std::move(detEvt);
      }

      events.push_back(move(evt));
      processed++;
      if ((m_progressInterval > 0 && processed % m_progressInterval == 0) ||
          processed == processingTarget) {
         printProgress(processed);
         lastProgress = processed;
      }
   }
   if (lastProgress != processed) printProgress(processed);
   cout << "\n[DUTAnalysis] processed " << events.size()
        << " single-track DUT events" << endl;

   // 运行DUT对齐
   if (m_runAlignment) {
      for (auto& det : duts) {
         RunDUTAlign(events, det, det->GetID());
      }
   }

   cout << "\nSaving DUT results..." << endl;

   // DUT数据分支
   Int_t dutID;
   Double_t resX, resY, predX, predY;
   Double_t hitX, hitY;
   vector<Int_t> clusterIndices;
   vector<ChannelHit> channelHits;
   vector<Cluster> clusters;
   Int_t hitFlag;
   Cluster clusterX, clusterY;
   vector<ChannelHit> channelHitsX, channelHitsY;

   tDut->Branch("eventID", &eventID);
   tDut->Branch("dutID", &dutID);
   tDut->Branch("predX", &predX);
   tDut->Branch("predY", &predY);
   tDut->Branch("resX", &resX);
   tDut->Branch("resY", &resY);
   tDut->Branch("hitX", &hitX);
   tDut->Branch("hitY", &hitY);
   tDut->Branch("clusterIndices", &clusterIndices);
   tDut->Branch("channelHits", &channelHits);
   tDut->Branch("clusters", &clusters);
   tDut->Branch("hitFlag", &hitFlag);
   tDut->Branch("clusterX", &clusterX);
   tDut->Branch("clusterY", &clusterY);
   tDut->Branch("channelHitsX", &channelHitsX);
   tDut->Branch("channelHitsY", &channelHitsY);

   for (auto& det : duts) {
      int id = det->GetID();
      for (auto& evt : events) {
         eventID = evt.eventID;
         dutID = id;
         hitX = DUTAnalysisConfig::kInvalidValue;
         hitY = DUTAnalysisConfig::kInvalidValue;
         resX = DUTAnalysisConfig::kInvalidValue;
         resY = DUTAnalysisConfig::kInvalidValue;
         hitFlag = 0;
         clusterIndices.clear();
         channelHits.clear();
         clusters.clear();
         clusterX = CreateInvalidCluster(DUTAnalysisConfig::kTypeX);
         clusterY = CreateInvalidCluster(DUTAnalysisConfig::kTypeY);
         channelHitsX.clear();
         channelHitsY.clear();

         TVector3 predG = det->CalcHitFromTrack(evt.track);
         TVector3 predL = det->GlobalToLocal(predG);
         predX = predL.X();
         predY = predL.Y();

         if (predY < m_effYMin || predY > m_effYMax) continue;
         // if (predY < 60 || predY > 66) continue;
         // if (predX > 70 || predX < 55) continue;
         // if (predY < 0 || predY > 20) continue;
         if (predX > m_effXMax || predX < m_effXMin) continue;

         auto frameIt = evt.detectorFramesMap.find(id);
         if (frameIt != evt.detectorFramesMap.end()) {
            const auto& detFrame = frameIt->second;

            channelHits = detFrame->ChannelHits();
            clusters = detFrame->Clusters();

            if (!clusters.empty()) {
               LocalHit localHit = CalcuDutResidual(det, clusters, predL, resX, resY);
               const auto& clusterIdx = localHit.clusterIndices;

               int idxX = clusterIdx[0];
               int idxY = clusterIdx[1];

               clusterX = (idxX >= 0) ? clusters[idxX] : CreateInvalidCluster(DUTAnalysisConfig::kTypeX);
               clusterY = (idxY >= 0) ? clusters[idxY] : CreateInvalidCluster(DUTAnalysisConfig::kTypeY);

               // 填充channelHitsX/Y
               channelHitsX.clear();
               channelHitsY.clear();
               if (idxX >= 0) {
                  for (int idx : clusterX.channelHitIndices) {
                     channelHitsX.push_back(channelHits[idx]);
                  }
               }
               if (idxY >= 0) {
                  for (int idx : clusterY.channelHitIndices) {
                     channelHitsY.push_back(channelHits[idx]);
                  }
               }

               hitX = localHit.localPos.X();
               hitY = localHit.localPos.Y();

               bool hasX = (idxX >= 0);
               bool hasY = (idxY >= 0);
               if (hasX && hasY)
                  hitFlag = 3;
               else if (hasX)
                  hitFlag = 1;
               else if (hasY)
                  hitFlag = 2;
               else
                  hitFlag = 0;

               clusterIndices = clusterIdx;

            } else {
               clusterX = CreateInvalidCluster(DUTAnalysisConfig::kTypeX);
               clusterY = CreateInvalidCluster(DUTAnalysisConfig::kTypeY);
               hitX = DUTAnalysisConfig::kInvalidValue;
               hitY = DUTAnalysisConfig::kInvalidValue;
               resX = DUTAnalysisConfig::kInvalidValue;
               resY = DUTAnalysisConfig::kInvalidValue;
               hitFlag = 0;
               clusterIndices.clear();
            }
         }

         tDut->Fill();
      }
   }

   tDut->Write();

   tDut->SetBranchAddress("dutID", &dutID);
   tDut->SetBranchAddress("predX", &predX);
   tDut->SetBranchAddress("predY", &predY);
   tDut->SetBranchAddress("hitX", &hitX);
   tDut->SetBranchAddress("hitY", &hitY);
   tDut->SetBranchAddress("resX", &resX);
   tDut->SetBranchAddress("resY", &resY);

   printf("\n%-5s %8s %10s %10s %11s %10s %11s %11s %12s\n",
          "DUT", "Events", "Xres[um]", "Yres[um]", "EventEff[%]",
          "Fake2D[%]", "NonUniX[%]", "NonUniY[%]", "NonUni2D[%]");

   for (const auto& det : duts) {
      const int id = det->GetID();
      vector<double> residualsX, residualsY;
      vector<DUTAlignmentQAPoint> alignmentPoints;
      for (Long64_t entry = 0; entry < tDut->GetEntries(); ++entry) {
         tDut->GetEntry(entry);
         if (dutID != id) continue;
         const bool validX = std::isfinite(hitX) && std::isfinite(resX) &&
                             hitX != DUTAnalysisConfig::kInvalidValue &&
                             resX != DUTAnalysisConfig::kInvalidValue;
         const bool validY = std::isfinite(hitY) && std::isfinite(resY) &&
                             hitY != DUTAnalysisConfig::kInvalidValue &&
                             resY != DUTAnalysisConfig::kInvalidValue;
         residualsX.push_back(resX);
         residualsY.push_back(resY);
         if (validX || validY) {
            alignmentPoints.push_back(
                {predX, predY, hitX, hitY, resX, resY, validX, validY});
         }
      }

      WriteDUTAlignmentQA(fDut, id, alignmentPoints,
                          m_effXMin, m_effXMax, m_effYMin, m_effYMax);
      const auto resolutionX = AnalyzeResidualSequence(
          residualsX, fDut, Form("hResX_DUT%d", id));
      const auto resolutionY = AnalyzeResidualSequence(
          residualsY, fDut, Form("hResY_DUT%d", id));
      const auto efficiency = DUTEfficiency::Analyze(
          events, strictSingleHitTrackerEvents, det,
          m_efficiencyConfig, fDut);

      printf("%-5d %8lld %10.1f %10.1f %11.3f %10.4f %11.3f %11.3f %12.3f\n",
             id, efficiency.eligibleEvents,
             1000.0 * resolutionX.sigmaFit, 1000.0 * resolutionY.sigmaFit,
             100.0 * efficiency.eventWeighted2D,
             100.0 * efficiency.fake2D,
             100.0 * efficiency.nonuniformityX,
             100.0 * efficiency.nonuniformityY,
             100.0 * efficiency.nonuniformity2D);
   }

   fDut->Close();
   delete fDut;

   f->Close();
   delete f;

   cout << "DUT data saved to: " << dutFile << endl;

   auto t1 = chrono::high_resolution_clock::now();
   double sec = chrono::duration<double>(t1 - t0).count();

   cout << "[DUTAnalysis] wrote DUTInfo.root in "
        << fixed << setprecision(2) << sec << " s" << endl;

   return true;
}

void RunDUTAlign(const std::vector<Event>& events, std::shared_ptr<Detector> detector, int detID) {

   if (events.empty()) {
      std::cerr << "[DUT Alignment] No events to analyze for DUT " << detID << "!" << std::endl;
      return;
   }

   std::cout << "[DUT " << detID << "] Aligning (6-parameter)..." << std::endl;
   auto minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
   minimizer->SetTolerance(0.005);
   minimizer->SetPrintLevel(0);

   UInt_t nPar = 6;

   // 使用lambda捕获this和参数
   auto chi2Func = [&events, &detector, detID](const double* par) -> double {
      return DUTChi2Objective(par, events, detector, detID);
   };

   ROOT::Math::Functor f(chi2Func, nPar);
   minimizer->SetFunction(f);

   for (UInt_t i = 0; i < nPar; ++i) {
      minimizer->SetVariable(i, Form("p%d", i), 0.0, 0.001);
   }

   minimizer->Minimize();

   // 应用结果
   const double* result = minimizer->X();
   double dx = result[0];
   double dy = result[1];
   double dz = result[2];
   double rotX = result[3];
   double rotY = result[4];
   double rotZ = result[5];

   detector->SetAlignment(dx, dy, dz, rotX, rotY, rotZ);

   TVector3 pos = detector->GetPos();
   TVector3 rot = detector->GetRot();

   std::cout << "DUT " << detID << " alignment: "
             << std::fixed << std::setprecision(5)
             << "\"position\": [" << pos.X() << "," << pos.Y() << "," << pos.Z() << "],"
             << "\"rotation\": [" << rot.X() << "," << rot.Y() << "," << rot.Z() << "]"
             << std::endl;
   std::cout << "DUT " << detID << " alignment corrections: "
             << "dx=" << dx << ", dy=" << dy << ", dz=" << dz << ", "
             << "rotX=" << rotX << ", rotY=" << rotY << ", rotZ=" << rotZ
             << std::endl;

   std::cout << "[DUT " << detID << "] Final chi2: " << minimizer->MinValue() << std::endl;

   delete minimizer;
}

// ========== DUT对齐私有方法 ==========

LocalHit CalcuDutResidual(std::shared_ptr<Detector> detector, const std::vector<Cluster>& clusters, const TVector3& predL, double& residualX, double& residualY) {

   double predX = predL.X();
   double predY = predL.Y();

   // 从探测器配置读取参数，消除魔法数字
   const auto* config = detector->GetPlanarConfig();
   const int typeX = DUTAnalysisConfig::kTypeX;
   const int typeY = DUTAnalysisConfig::kTypeY;

   LocalHit localHit;
   if (!config) {
      residualX = DUTAnalysisConfig::kInvalidValue;
      residualY = DUTAnalysisConfig::kInvalidValue;
      localHit.localPos.SetXYZ(DUTAnalysisConfig::kInvalidValue,
                               DUTAnalysisConfig::kInvalidValue, 0);
      localHit.clusterIndices = {-1, -1};
      return localHit;
   }

   // X方向处理：找到最优cluster
   int bestClusterXIndex = -1;
   double minResX = std::numeric_limits<double>::infinity();
   double bestPosX = DUTAnalysisConfig::kInvalidValue;
   if (config->readoutPlanePitch.find(typeX) != config->readoutPlanePitch.end()) {
      double pitchX = config->readoutPlanePitch.at(typeX);
      for (size_t i = 0; i < clusters.size(); ++i) {
         if (clusters[i].type == typeX) {
            double currentResX = std::abs(clusters[i].centroid * pitchX - predX);
            if (currentResX < minResX) {
               minResX = currentResX;
               bestClusterXIndex = static_cast<int>(i);
            }
         }
      }
   }

   int bestClusterYIndex = -1;
   double minResY = std::numeric_limits<double>::infinity();
   double bestPosY = DUTAnalysisConfig::kInvalidValue;
   if (config->readoutPlanePitch.find(typeY) != config->readoutPlanePitch.end()) {
      double pitchY = config->readoutPlanePitch.at(typeY);
      for (size_t i = 0; i < clusters.size(); ++i) {
         if (clusters[i].type == typeY) {
            double currentResY = std::abs(clusters[i].centroid * pitchY - predY);
            if (currentResY < minResY) {
               minResY = currentResY;
               bestClusterYIndex = static_cast<int>(i);
            }
         }
      }
   }

   // 构建LocalHit
   if (bestClusterXIndex != -1) {
      double pitchX = config->readoutPlanePitch.at(typeX);
      bestPosX = clusters[bestClusterXIndex].pos * pitchX;
      residualX = bestPosX - predX;
   } else {
      bestPosX = DUTAnalysisConfig::kInvalidValue;
      residualX = DUTAnalysisConfig::kInvalidValue;
   }

   if (bestClusterYIndex != -1) {
      double pitchY = config->readoutPlanePitch.at(typeY);
      bestPosY = clusters[bestClusterYIndex].pos * pitchY;
      residualY = bestPosY - predY;
   } else {
      bestPosY = DUTAnalysisConfig::kInvalidValue;
      residualY = DUTAnalysisConfig::kInvalidValue;
   }

   localHit.localPos.SetXYZ(bestPosX, bestPosY, 0);
   localHit.clusterIndices = {bestClusterXIndex, bestClusterYIndex};

   return localHit;
}

double DUTChi2Objective(const double* par, const std::vector<Event>& events, std::shared_ptr<Detector> detector, int detID) {

   const double dx = par[0];
   const double dy = par[1];
   const double dz = par[2];
   const double rotX = par[3];
   const double rotY = par[4];
   const double rotZ = par[5];

   detector->SetAlignment(dx, dy, dz, rotX, rotY, rotZ);

   // 对所有事件求平均 χ²
   double chi2 = 0.0;
   int nEvents = 0;

   for (const auto& evt : events) {
      double residualX = 0.0;
      double residualY = 0.0;

      // 使用detectorFramesMap获取Clusters
      auto frameIt = evt.detectorFramesMap.find(detID);
      if (frameIt == evt.detectorFramesMap.end() || frameIt->second->Clusters().empty()) {
         continue;
      }

      const auto& clusters = frameIt->second->Clusters();
      TVector3 predG = detector->CalcHitFromTrack(evt.track);
      TVector3 predL = detector->GlobalToLocal(predG);

      double predX = predL.X();
      double predY = predL.Y();

      // if (predY < 10 || predY > 50) continue;
      // if (predX < 60 & predX > 40) continue;
      // if (predX > 90) continue;
      // 调用新的CalcuDutResidual方法
      LocalHit localHit = CalcuDutResidual(detector, clusters, predL, residualX, residualY);

      double res = residualX * residualX + residualY * residualY;
      if (res > 3) continue;
      chi2 += res;
      nEvents++;
   }

   // 返回平均 χ²
   return (nEvents > 0) ? chi2 / nEvents : 1e9;
}

Cluster DUTAnalysisScript::CreateInvalidCluster(int type) {
   Cluster invalidCluster;
   invalidCluster.type = type;
   invalidCluster.size = DUTAnalysisConfig::kInvalidSize;
   invalidCluster.range = DUTAnalysisConfig::kInvalidSize;
   invalidCluster.charge = DUTAnalysisConfig::kInvalidValue;
   invalidCluster.maxAmp = DUTAnalysisConfig::kInvalidValue;
   invalidCluster.time = DUTAnalysisConfig::kInvalidValue;
   invalidCluster.centroid = DUTAnalysisConfig::kInvalidValue;
   invalidCluster.pos = DUTAnalysisConfig::kInvalidValue;
   invalidCluster.channelHitIndices.clear();
   return invalidCluster;
}

ResidualAnalysisResult AnalyzeResidualSequence(const std::vector<double>& rawRes, TFile* outputFile, const std::string& histName) {

   double rangeSigma = 3.0;
   int bins = 100;
   bool drawPlots = true;

   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);

   ResidualAnalysisResult result;
   if (rawRes.empty() || !outputFile) return result;

   std::vector<double> residuals;
   for (double res : rawRes) {
      if (abs(res) < 5)
         residuals.push_back(res);
   }
   if (residuals.empty()) return result;

   const int totalEvents = rawRes.size();
   const int N = residuals.size();
   /* ================= Units ================= */
   const double unit = 1.0;  //

   /* ================= Basic statistics ================= */
   double mean = 0.0;
   for (double r : residuals) mean += r * unit;
   mean /= N;

   double rms = 0.0;
   for (double r : residuals) {
      const double x = r * unit - mean;
      rms += x * x;
   }
   rms = std::sqrt(rms / N);

   std::sort(residuals.begin(), residuals.end());

   auto quantile = [&](double q) {
      int idx = std::lround(q * (residuals.size() - 1));
      return residuals[idx];
   };

   double q16 = quantile(0.16);
   double q50 = quantile(0.50);
   double q84 = quantile(0.84);

   double sigma68 = 0.5 * (q84 - q16);

   const double xmin = mean - rangeSigma * rms;
   const double xmax = mean + rangeSigma * rms;

   outputFile->cd();

   /* ================= Histogram ================= */
   TH1D* hist = new TH1D(histName.c_str(), (histName + ";Residual[mm];Events").c_str(), bins, xmin, xmax);

   for (double r : residuals)
      hist->Fill(r * unit);

   hist->SetLineColor(kBlack);
   hist->SetMarkerStyle(20);
   hist->SetMarkerSize(0.8);
   // hist->Scale(1.0 / hist->Integral("width"));

   TLatex* title = new TLatex(0.12, 0.82, "**Experiment**Preliminary");
   title->SetNDC();
   title->SetTextSize(0.04);
   title->SetTextFont(62);
   hist->GetListOfFunctions()->Add(title);

   /* ====================================================
    *  Single Gaussian
    * ==================================================== */
   TF1* fSingle = new TF1("fSingle", "gaus", xmin, xmax);
   fSingle->SetLineColor(kRed);
   fSingle->SetLineWidth(2);

   hist->Fit(fSingle, "QR0");

   const double muS = fSingle->GetParameter(1);
   const double sigS = fSingle->GetParameter(2);

   int n3S = 0, n5S = 0, n3mm = 0, n1mm = 0;
   for (double r : residuals) {
      const double x = r * unit;
      if (std::fabs(x - muS) < 3 * sigS) n3S++;
      if (std::fabs(x - muS) < 5 * sigS) n5S++;
      if (std::fabs(x - muS) < 0.5) n3mm++;
      if (std::fabs(x - muS) < 1.0) n1mm++;
   }

   const double eff3S = double(n3S) / totalEvents;
   const double eff5S = double(n5S) / totalEvents;
   const double eff3mm = double(n3mm) / totalEvents;
   const double eff1mm = double(n1mm) / totalEvents;

   TLegend* legSingle = new TLegend(0.55, 0.55, 0.88, 0.88);
   legSingle->SetBorderSize(1);
   legSingle->SetFillColor(kWhite);
   legSingle->SetTextFont(42);
   legSingle->SetTextSize(0.027);
   legSingle->SetTextAlign(12);

   legSingle->AddEntry(fSingle, "Single Gaussian Fit", "l");
   legSingle->AddEntry((TObject*)0, Form("#sigma_{68} = %.4f mm", sigma68), "");
   legSingle->AddEntry((TObject*)0, Form("RMS (|x|<5 mm) = %.4f mm", rms), "");
   legSingle->AddEntry((TObject*)0, Form("#mu = %.4f #pm %.4f mm", muS, fSingle->GetParError(1)), "");
   legSingle->AddEntry((TObject*)0, Form("#sigma = %.4f #pm %.4f mm", sigS, fSingle->GetParError(2)), "");
   legSingle->AddEntry((TObject*)0, Form("Eff(|x|<3#sigma) = %.2f%%", eff3S * 100), "");
   legSingle->AddEntry((TObject*)0, Form("Eff(|x|<5#sigma) = %.2f%%", eff5S * 100), "");
   legSingle->AddEntry((TObject*)0, Form("Eff(|x|<1 mm) = %.2f%%", eff1mm * 100), "");

   /* ---------- Single Gaussian canvas ---------- */
   if (drawPlots) {
      TCanvas* cSingle =
          new TCanvas((histName + "Single").c_str(), ("Residual - " + histName + " Single Gaussian").c_str(), 800, 600);

      hist->Draw("E");
      fSingle->Draw("same");
      legSingle->Draw();

      cSingle->Write();
   }

   /* ====================================================
    *  Double Gaussian
    * ==================================================== */
   TF1* fDouble = new TF1(
       "fDouble",
       "[0]*TMath::Gaus(x,[1],[2],true) + "
       "[3]*TMath::Gaus(x,[4],[5],true)",
       xmin, xmax);

   fDouble->SetParameters(
       0.7 * hist->GetMaximum(), muS, 0.8 * sigS, 0.3 * hist->GetMaximum(), muS, 2.5 * sigS);

   // fDouble->SetParLimits(2, 0.2 * sigS, 2.0 * sigS);
   // fDouble->SetParLimits(4, 1.5 * sigS, 10.0 * sigS);

   fDouble->SetLineColor(kRed);
   fDouble->SetLineWidth(2);

   hist->Fit(fDouble, "RQ0");

   const double A1 = fDouble->GetParameter(0);
   const double mu1 = fDouble->GetParameter(1);
   const double sig1 = fDouble->GetParameter(2);

   const double A2 = fDouble->GetParameter(3);
   const double mu2 = fDouble->GetParameter(4);
   const double sig2 = fDouble->GetParameter(5);

   const double w1 = A1 / (A1 + A2);
   const double w2 = A2 / (A1 + A2);

   const double sigmaEff =
       sqrt(w1 * sig1 * sig1 +
            w2 * sig2 * sig2 +
            w1 * w2 * (mu1 - mu2) * (mu1 - mu2));

   int n5Tail = 0;
   for (double r : residuals) {
      const double x = r * unit;
      if (std::fabs(x - mu1) < 5 * sig2) n5Tail++;
   }

   const double eff5Tail = double(n5Tail) / totalEvents;

   TLegend* legDouble = new TLegend(0.60, 0.55, 0.88, 0.88);
   legDouble->SetBorderSize(1);
   legDouble->SetFillColor(kWhite);
   legDouble->SetTextFont(42);
   legDouble->SetTextSize(0.027);
   legDouble->SetTextAlign(12);

   legDouble->AddEntry(fDouble, "Double Gaussian Fit", "l");
   legDouble->AddEntry((TObject*)0, Form("Std Dev = %.4f mm", rms), "");
   legDouble->AddEntry((TObject*)0, Form("#mu_{core} = %.4f #pm %.4f mm", muS, fSingle->GetParError(1)), "");
   legDouble->AddEntry((TObject*)0,
                       Form("#sigma_{core} = %.4f #pm %.4f mm", sig1, fDouble->GetParError(2)), "");
   legDouble->AddEntry((TObject*)0, Form("#mu_{tail} = %.4f #pm %.4f mm", mu2, fDouble->GetParError(4)), "");
   legDouble->AddEntry((TObject*)0,
                       Form("#sigma_{tail} = %.4f #pm %.4f mm", sig2, fDouble->GetParError(5)), "");
   legDouble->AddEntry((TObject*)0,
                       Form("Tail fraction = %.3f", w2), "");
   legDouble->AddEntry((TObject*)0,
                       Form("#sigma_{eff} = %.4f mm", sigmaEff), "");
   legDouble->AddEntry((TObject*)0,
                       Form("Eff(|x|<5#sigma_{tail}) = %.2f%%",
                            eff5Tail * 100),
                       "");

   /* ---------- Double Gaussian canvas ---------- */
   if (drawPlots) {
      TCanvas* cDouble =
          new TCanvas((histName + "Double").c_str(), ("Residual - " + histName + " Double Gaussian").c_str(), 800, 600);

      hist->Draw("E");
      fDouble->Draw("same");
      legDouble->Draw();
      cDouble->Write();
   }

   // Fill result struct
   result.sigma68 = sigma68;
   result.sigmaFit = sigS;
   result.sigmaFitError = fSingle->GetParError(2);
   result.mean = muS;
   result.efficiency = double(N) / totalEvents;
   result.eff3S = eff3S;
   result.eff5S = eff5S;
   result.eff1mm = eff1mm;

   return result;
}

REGISTER_SCRIPT("DUTAnalysis", DUTAnalysisScript);
