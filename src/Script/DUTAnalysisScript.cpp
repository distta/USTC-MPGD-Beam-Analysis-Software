#include "Script/DUTAnalysisScript.h"
#include "Terminal.h"
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
#include <memory>
#include <numeric>
#include <unordered_map>
#include <unordered_set>

using namespace std;

struct DUTAlignmentRegion {
   double xMin;
   double xMax;
   double yMin;
   double yMax;
};

void RunDUTAlign(const std::vector<Event>& events,
                 std::shared_ptr<Detector> detector, int detID,
                 const DUTAlignmentRegion& region);
LocalHit CalcuDutResidual(std::shared_ptr<Detector> detector,
                          const std::vector<Cluster>& clusters,
                          const TVector3& predL, double& residualX,
                          double& residualY, bool useCentroid = false);
double DUTChi2Objective(const double* par, const std::vector<Event>& events,
                        std::shared_ptr<Detector> detector, int detID,
                        const DUTAlignmentRegion& region,
                        double matchingSigma, double minMatchDominance);
double DUTCCResidualSlopeY(const double* par,
                           const std::vector<Event>& events,
                           std::shared_ptr<Detector> detector, int detID,
                           const DUTAlignmentRegion& region,
                           double matchingSigma);

std::pair<int, int> PlanarAxisTypes(const planarConfig& config) {
   int typeX = -1;
   int typeY = -1;
   double bestX = std::numeric_limits<double>::infinity();
   double bestY = std::numeric_limits<double>::infinity();
   for (int type : config.readoutPlaneType) {
      const double angle =
          config.readoutPlaneAngle.at(type) * TMath::DegToRad();
      const double xScore = std::abs(std::sin(angle));
      const double yScore = std::abs(std::cos(angle));
      if (xScore < bestX) {
         bestX = xScore;
         typeX = type;
      }
      if (yScore < bestY) {
         bestY = yScore;
         typeY = type;
      }
   }
   return {typeX, typeY};
}

struct DUTAlignmentMatch {
   double residualX = DUTAnalysisConfig::kInvalidValue;
   double residualY = DUTAnalysisConfig::kInvalidValue;
   double dominanceX = 0.0;
   double dominanceY = 0.0;
};

// Smooth nearest-cluster association used only by DUT alignment.  All
// candidates are reconsidered at every geometry evaluation, so the fit never
// freezes a possibly wrong cluster set selected from the initial geometry.
DUTAlignmentMatch CalcuDutAlignmentMatch(
    std::shared_ptr<Detector> detector,
    const std::vector<Cluster>& clusters, const TVector3& predL,
    double matchingSigma) {
   DUTAlignmentMatch match;
   const auto* config = detector->GetPlanarConfig();
   if (!config || !(matchingSigma > 0.0) || !std::isfinite(matchingSigma))
      return match;

   const auto [typeX, typeY] = PlanarAxisTypes(*config);
   const auto axisResidual = [&](int type, double prediction,
                                 double& dominance) {
      const auto pitchEntry = config->readoutPlanePitch.find(type);
      if (pitchEntry == config->readoutPlanePitch.end())
         return DUTAnalysisConfig::kInvalidValue;

      std::vector<double> residuals;
      residuals.reserve(clusters.size());
      double bestResidual2 = std::numeric_limits<double>::infinity();
      for (const auto& cluster : clusters) {
         if (cluster.type != type || !std::isfinite(cluster.centroid))
            continue;
         const double residual =
             cluster.centroid * pitchEntry->second - prediction;
         if (!std::isfinite(residual)) continue;
         residuals.push_back(residual);
         bestResidual2 = std::min(bestResidual2, residual * residual);
      }
      if (residuals.empty()) return DUTAnalysisConfig::kInvalidValue;

      const double inverseTwoSigma2 =
          0.5 / (matchingSigma * matchingSigma);
      double sumWeight = 0.0;
      double weightedResidual = 0.0;
      double bestWeight = 0.0;
      for (double residual : residuals) {
         // Subtracting the best exponent keeps the calculation stable even
         // when the configured geometry is still far from the final one.
         const double weight = std::exp(
             -(residual * residual - bestResidual2) * inverseTwoSigma2);
         sumWeight += weight;
         weightedResidual += weight * residual;
         bestWeight = std::max(bestWeight, weight);
      }
      if (!(sumWeight > 0.0) || !std::isfinite(sumWeight))
         return DUTAnalysisConfig::kInvalidValue;
      dominance = bestWeight / sumWeight;
      return weightedResidual / sumWeight;
   };

   match.residualX = axisResidual(typeX, predL.X(), match.dominanceX);
   match.residualY = axisResidual(typeY, predL.Y(), match.dominanceY);
   return match;
}
struct ResidualAnalysisResult {
   double sigma68 = 0;
   double sigmaFit = 0;
   double sigmaFitError = 0;
   double mean = 0;
   double rms = 0;
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

std::pair<double, double> DUTDistributionRange(
    std::vector<double> values) {
   values.erase(
       std::remove_if(
           values.begin(), values.end(),
           [](double value) { return !std::isfinite(value); }),
       values.end());
   if (values.empty()) return {0.0, 1.0};
   std::sort(values.begin(), values.end());
   const size_t lowIndex = static_cast<size_t>(
       0.005 * static_cast<double>(values.size() - 1));
   const size_t highIndex = static_cast<size_t>(
       0.995 * static_cast<double>(values.size() - 1));
   double minimum = values[lowIndex];
   double maximum = values[highIndex];

   double low = std::min(0.0, minimum);
   double high = std::max(0.0, maximum);
   if (!(high > low)) return {low - 0.5, high + 0.5};
   const double padding = 0.05 * (high - low);
   if (low < 0.0) low -= padding;
   high += padding;
   return {low, high};
}

void StyleDUTDistribution(TH1D& histogram, Color_t color) {
   histogram.SetDirectory(nullptr);
   histogram.SetStats(false);
   histogram.SetLineColor(color);
   histogram.SetLineWidth(2);
}

void WritePlanarDUTHitProperties(
    TDirectory* detectorDirectory, int dutID,
    const std::vector<double>& clusterChargesX,
    const std::vector<double>& clusterChargesY,
    const std::vector<int>& clusterSizesX,
    const std::vector<int>& clusterSizesY,
    const std::vector<double>& hitADCsX,
    const std::vector<double>& hitADCsY) {
   if (!detectorDirectory) return;
   auto* hitProperties = detectorDirectory->GetDirectory("HitProperties");
   if (!hitProperties)
      hitProperties = detectorDirectory->mkdir("HitProperties");
   if (!hitProperties) return;
   TDirectory::TContext context(hitProperties);

   std::vector<double> combinedCharges = clusterChargesX;
   combinedCharges.insert(
       combinedCharges.end(), clusterChargesY.begin(),
       clusterChargesY.end());
   const auto chargeRange = DUTDistributionRange(combinedCharges);
   TH1D chargeX(
       "hClusterChargeX",
       ("DUT " + std::to_string(dutID) +
        " track-matched cluster ADC sum;"
        "Cluster charge = #Sigma strip ADC [ADC];Entries")
           .c_str(),
       200, chargeRange.first, chargeRange.second);
   TH1D chargeY(
       "hClusterChargeY", "", 200,
       chargeRange.first, chargeRange.second);
   StyleDUTDistribution(chargeX, kBlue + 1);
   StyleDUTDistribution(chargeY, kRed + 1);
   for (double value : clusterChargesX) chargeX.Fill(value);
   for (double value : clusterChargesY) chargeY.Fill(value);
   chargeX.SetMaximum(
       1.15 * std::max(chargeX.GetMaximum(), chargeY.GetMaximum()));
   TCanvas chargeCanvas(
       "cClusterCharge",
       "Track-matched DUT cluster ADC sum", 900, 700);
   chargeX.Draw("HIST");
   chargeY.Draw("HIST SAME");
   TLegend chargeLegend(0.68, 0.75, 0.88, 0.88);
   chargeLegend.SetBorderSize(0);
   chargeLegend.AddEntry(&chargeX, "X strips", "l");
   chargeLegend.AddEntry(&chargeY, "Y strips", "l");
   chargeLegend.Draw();
   chargeCanvas.Write();

   std::vector<int> combinedClusterSizes = clusterSizesX;
   combinedClusterSizes.insert(
       combinedClusterSizes.end(), clusterSizesY.begin(),
       clusterSizesY.end());
   std::sort(combinedClusterSizes.begin(), combinedClusterSizes.end());
   int maximumClusterSize = 1;
   if (!combinedClusterSizes.empty()) {
      const size_t highIndex = static_cast<size_t>(
          0.995 *
          static_cast<double>(combinedClusterSizes.size() - 1));
      maximumClusterSize =
          std::max(1, combinedClusterSizes[highIndex]);
   }
   TH1D sizeX(
       "hClusterSizeX",
       ("DUT " + std::to_string(dutID) +
        " track-matched cluster size;Cluster size [strips];Entries")
           .c_str(),
       maximumClusterSize, 0.5, maximumClusterSize + 0.5);
   TH1D sizeY(
       "hClusterSizeY", "", maximumClusterSize,
       0.5, maximumClusterSize + 0.5);
   StyleDUTDistribution(sizeX, kBlue + 1);
   StyleDUTDistribution(sizeY, kRed + 1);
   for (int value : clusterSizesX) sizeX.Fill(value);
   for (int value : clusterSizesY) sizeY.Fill(value);
   sizeX.SetMaximum(
       1.15 * std::max(sizeX.GetMaximum(), sizeY.GetMaximum()));
   TCanvas sizeCanvas(
       "cClusterSize",
       "Track-matched DUT cluster size", 900, 700);
   sizeX.Draw("HIST");
   sizeY.Draw("HIST SAME");
   TLegend sizeLegend(0.68, 0.75, 0.88, 0.88);
   sizeLegend.SetBorderSize(0);
   sizeLegend.AddEntry(&sizeX, "X strips", "l");
   sizeLegend.AddEntry(&sizeY, "Y strips", "l");
   sizeLegend.Draw();
   sizeCanvas.Write();

   std::vector<double> combinedADCs = hitADCsX;
   combinedADCs.insert(
       combinedADCs.end(), hitADCsY.begin(), hitADCsY.end());
   const auto adcRange = DUTDistributionRange(combinedADCs);
   TH1D adcX(
       "hHitADCX",
       ("DUT " + std::to_string(dutID) +
        " track-matched hit amplitude;Hit amplitude [ADC];Entries")
           .c_str(),
       200, adcRange.first, adcRange.second);
   TH1D adcY(
       "hHitADCY", "", 200, adcRange.first, adcRange.second);
   StyleDUTDistribution(adcX, kBlue + 1);
   StyleDUTDistribution(adcY, kRed + 1);
   for (double value : hitADCsX) adcX.Fill(value);
   for (double value : hitADCsY) adcY.Fill(value);
   adcX.SetMaximum(
       1.15 * std::max(adcX.GetMaximum(), adcY.GetMaximum()));
   TCanvas adcCanvas(
       "cHitADC",
       "Track-matched DUT hit amplitude", 900, 700);
   adcX.Draw("HIST");
   adcY.Draw("HIST SAME");
   TLegend adcLegend(0.68, 0.75, 0.88, 0.88);
   adcLegend.SetBorderSize(0);
   adcLegend.AddEntry(&adcX, "X strips", "l");
   adcLegend.AddEntry(&adcY, "Y strips", "l");
   adcLegend.Draw();
   adcCanvas.Write();
}

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

   const json alignmentCfg =
       config.value("alignmentRegion", json::object());
   m_alignmentXMin = alignmentCfg.value("xMin", m_alignmentXMin);
   m_alignmentXMax = alignmentCfg.value("xMax", m_alignmentXMax);
   m_alignmentYMin = alignmentCfg.value("yMin", m_alignmentYMin);
   m_alignmentYMax = alignmentCfg.value("yMax", m_alignmentYMax);
   if (m_alignmentXMax < m_alignmentXMin)
      std::swap(m_alignmentXMax, m_alignmentXMin);
   if (m_alignmentYMax < m_alignmentYMin)
      std::swap(m_alignmentYMax, m_alignmentYMin);

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
   double sigTime = numeric_limits<double>::quiet_NaN();
   Bool_t hasTrackTime = false;

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
   const bool hasT0Branch = trackTree->GetBranch("t0") != nullptr;
   const bool hasT0FlagBranch = trackTree->GetBranch("hasT0") != nullptr;
   if (hasT0Branch) trackTree->SetBranchAddress("t0", &sigTime);
   if (hasT0FlagBranch)
      trackTree->SetBranchAddress("hasT0", &hasTrackTime);

   std::vector<Event> events;  // Script本地数据

   int processed = 0;


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

      if (!hasT0FlagBranch)
         hasTrackTime = hasT0Branch && isfinite(sigTime);

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
         if (hasTrackTime)
            detEvt->Process(evt.t0);
         else
            detEvt->Process();
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
   cout << '\n';

   // 运行DUT对齐
   if (m_runAlignment) {
      for (auto& det : duts) {
         RunDUTAlign(events, det, det->GetID(),
                     {m_alignmentXMin, m_alignmentXMax,
                      m_alignmentYMin, m_alignmentYMax});
      }
   }


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
      const auto* planar = det->GetPlanarConfig();
      const auto [typeX, typeY] =
          planar ? PlanarAxisTypes(*planar)
                 : std::pair<int, int>{
                       DUTAnalysisConfig::kTypeX,
                       DUTAnalysisConfig::kTypeY};
      int hitMapXBins = m_effXBins;
      int hitMapYBins = m_effYBins;
      double hitMapXMin = m_effXMin;
      double hitMapXMax = m_effXMax;
      double hitMapYMin = m_effYMin;
      double hitMapYMax = m_effYMax;
      if (planar) {
         const auto setReadoutBinning =
             [&](int type, int& bins, double& minimum, double& maximum) {
                const auto pitch = planar->readoutPlanePitch.find(type);
                const auto strips =
                    planar->readoutPlaneStripNumber.find(type);
                if (pitch == planar->readoutPlanePitch.end() ||
                    strips == planar->readoutPlaneStripNumber.end() ||
                    !(pitch->second > 0.0) || strips->second <= 0) {
                   return;
                }
                bins = strips->second;
                minimum = 0.0;
                maximum = pitch->second * strips->second;
             };
         setReadoutBinning(
             typeX, hitMapXBins, hitMapXMin, hitMapXMax);
         setReadoutBinning(
             typeY, hitMapYBins, hitMapYMin, hitMapYMax);
      }
      TH2D reconstructedHitMap(
          "hReconstructedHitMap",
          "All reconstructed DUT hits;Reconstructed X [mm];"
          "Reconstructed Y [mm]",
          hitMapXBins, hitMapXMin, hitMapXMax,
          hitMapYBins, hitMapYMin, hitMapYMax);
      TH2D matchedHitMap(
          "hMatchedHitMap",
          "Track-selected reconstructed DUT hits;Reconstructed X [mm];"
          "Reconstructed Y [mm]",
          hitMapXBins, hitMapXMin, hitMapXMax,
          hitMapYBins, hitMapYMin, hitMapYMax);
      for (TH2D* histogram : {&reconstructedHitMap, &matchedHitMap}) {
         histogram->SetDirectory(nullptr);
         histogram->SetStats(false);
         histogram->SetOption("COLZ");
      }
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
         clusterX = CreateInvalidCluster(typeX);
         clusterY = CreateInvalidCluster(typeY);
         channelHitsX.clear();
         channelHitsY.clear();

         TVector3 predG = det->CalcHitFromTrack(evt.track);
         TVector3 predL = det->GlobalToLocal(predG);
         predX = predL.X();
         predY = predL.Y();

         auto frameIt = evt.detectorFramesMap.find(id);
         if (frameIt != evt.detectorFramesMap.end()) {
            for (const auto& localHit : frameIt->second->LocalHits()) {
               const double localX = localHit.localPos.X();
               const double localY = localHit.localPos.Y();
               if (std::isfinite(localX) && std::isfinite(localY))
                  reconstructedHitMap.Fill(localX, localY);
            }
         }

         if (predY < m_effYMin || predY > m_effYMax) continue;
         // if (predY < 60 || predY > 66) continue;
         // if (predX > 70 || predX < 55) continue;
         // if (predY < 0 || predY > 20) continue;
         if (predX > m_effXMax || predX < m_effXMin) continue;

         if (frameIt != evt.detectorFramesMap.end()) {
            const auto& detFrame = frameIt->second;

            channelHits = detFrame->ChannelHits();
            clusters = detFrame->Clusters();

            if (!clusters.empty()) {
               LocalHit localHit = CalcuDutResidual(det, clusters, predL, resX, resY);
               const auto& clusterIdx = localHit.clusterIndices;

               int idxX = clusterIdx[0];
               int idxY = clusterIdx[1];

               clusterX = (idxX >= 0) ? clusters[idxX] : CreateInvalidCluster(typeX);
               clusterY = (idxY >= 0) ? clusters[idxY] : CreateInvalidCluster(typeY);

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
               if (hasX && hasY && std::isfinite(hitX) &&
                   std::isfinite(hitY)) {
                  matchedHitMap.Fill(hitX, hitY);
               }

            } else {
               clusterX = CreateInvalidCluster(typeX);
               clusterY = CreateInvalidCluster(typeY);
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

      auto* detectorDirectory =
          fDut->GetDirectory(("DUT_" + std::to_string(id)).c_str());
      if (!detectorDirectory)
         detectorDirectory =
             fDut->mkdir(("DUT_" + std::to_string(id)).c_str());
      auto* hitMapDirectory = detectorDirectory->GetDirectory("HitMaps");
      if (!hitMapDirectory)
         hitMapDirectory = detectorDirectory->mkdir("HitMaps");
      if (hitMapDirectory) {
         TDirectory::TContext context(hitMapDirectory);
         reconstructedHitMap.Write();
         matchedHitMap.Write();

         std::unique_ptr<TH1D> reconstructedProfileX(
             reconstructedHitMap.ProjectionX(
                 "hReconstructedHitProfileX"));
         std::unique_ptr<TH1D> reconstructedProfileY(
             reconstructedHitMap.ProjectionY(
                 "hReconstructedHitProfileY"));
         std::unique_ptr<TH1D> matchedProfileX(
             matchedHitMap.ProjectionX("hMatchedHitProfileX"));
         std::unique_ptr<TH1D> matchedProfileY(
             matchedHitMap.ProjectionY("hMatchedHitProfileY"));

         reconstructedProfileX->SetTitle(
             "DUT hit profile X;Reconstructed X [mm];Entries");
         reconstructedProfileY->SetTitle(
             "DUT hit profile Y;Reconstructed Y [mm];Entries");
         for (TH1D* profile :
              {reconstructedProfileX.get(), reconstructedProfileY.get()}) {
            profile->SetDirectory(nullptr);
            profile->SetStats(false);
            profile->SetLineColor(kBlue + 1);
            profile->SetLineWidth(2);
            profile->Write();
         }
         for (TH1D* profile :
              {matchedProfileX.get(), matchedProfileY.get()}) {
            profile->SetDirectory(nullptr);
            profile->SetStats(false);
            profile->SetLineColor(kRed + 1);
            profile->SetLineWidth(2);
            profile->Fit("gaus", "Q");
            profile->Write();
         }

         TCanvas profileCanvas(
             "cHitProfiles", "DUT hit profiles", 1400, 650);
         profileCanvas.Divide(2, 1);
         profileCanvas.cd(1);
         reconstructedProfileX->SetMaximum(
             1.15 * std::max(reconstructedProfileX->GetMaximum(),
                             matchedProfileX->GetMaximum()));
         reconstructedProfileX->Draw("HIST");
         matchedProfileX->Draw("HIST SAME");
         TLegend legendX(0.62, 0.76, 0.88, 0.88);
         legendX.SetBorderSize(0);
         legendX.AddEntry(
             reconstructedProfileX.get(), "All reconstructed", "l");
         legendX.AddEntry(
             matchedProfileX.get(), "Track selected", "l");
         legendX.Draw();

         profileCanvas.cd(2);
         reconstructedProfileY->SetMaximum(
             1.15 * std::max(reconstructedProfileY->GetMaximum(),
                             matchedProfileY->GetMaximum()));
         reconstructedProfileY->Draw("HIST");
         matchedProfileY->Draw("HIST SAME");
         TLegend legendY(0.62, 0.76, 0.88, 0.88);
         legendY.SetBorderSize(0);
         legendY.AddEntry(
             reconstructedProfileY.get(), "All reconstructed", "l");
         legendY.AddEntry(
             matchedProfileY.get(), "Track selected", "l");
         legendY.Draw();
         profileCanvas.Write();
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

   for (const auto& det : duts) {
      const int id = det->GetID();
      vector<double> residualsX, residualsY;
      vector<DUTAlignmentQAPoint> alignmentPoints;
      vector<double> selectedClusterChargesX, selectedClusterChargesY;
      vector<int> selectedClusterSizesX, selectedClusterSizesY;
      vector<double> selectedHitADCsX, selectedHitADCsY;
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
         if (clusterX.size > 0) {
            if (std::isfinite(clusterX.charge))
               selectedClusterChargesX.push_back(clusterX.charge);
            selectedClusterSizesX.push_back(clusterX.size);
            for (const auto& hit : channelHitsX)
               if (hit.isValid && std::isfinite(hit.amp))
                  selectedHitADCsX.push_back(hit.amp);
         }
         if (clusterY.size > 0) {
            if (std::isfinite(clusterY.charge))
               selectedClusterChargesY.push_back(clusterY.charge);
            selectedClusterSizesY.push_back(clusterY.size);
            for (const auto& hit : channelHitsY)
               if (hit.isValid && std::isfinite(hit.amp))
                  selectedHitADCsY.push_back(hit.amp);
         }
         if (validX || validY) {
            alignmentPoints.push_back(
                {predX, predY, hitX, hitY, resX, resY, validX, validY});
         }
      }

      auto* detectorDirectory =
          fDut->GetDirectory(("DUT_" + std::to_string(id)).c_str());
      if (!detectorDirectory)
         detectorDirectory =
             fDut->mkdir(("DUT_" + std::to_string(id)).c_str());
      WritePlanarDUTHitProperties(
          detectorDirectory, id,
          selectedClusterChargesX, selectedClusterChargesY,
          selectedClusterSizesX, selectedClusterSizesY,
          selectedHitADCsX, selectedHitADCsY);
      WriteDUTAlignmentQA(fDut, id, alignmentPoints,
                          m_effXMin, m_effXMax, m_effYMin, m_effYMax);
      const auto resolutionX = AnalyzeResidualSequence(
          residualsX, fDut, Form("hResX_DUT%d", id));
      const auto resolutionY = AnalyzeResidualSequence(
          residualsY, fDut, Form("hResY_DUT%d", id));
      const auto efficiency = DUTEfficiency::Analyze(
          events, strictSingleHitTrackerEvents, det,
          m_efficiencyConfig, fDut);

      const auto mean = [](const auto& values) {
         if (values.empty()) return 0.0;
         return std::accumulate(values.begin(), values.end(), 0.0) /
                static_cast<double>(values.size());
      };
      const auto row = [](const string& label, const string& value) {
         cout << "  " << left << setw(38) << label << right << setw(39)
              << value << '\n';
      };
      const auto separator = [] { cout << string(100, '-') << '\n'; };
      const auto fixedValue = [](double value, int precision,
                                 const string& suffix = "") {
         ostringstream text;
         text << fixed << setprecision(precision) << value << suffix;
         return text.str();
      };
      const double efficiencyError =
          efficiency.eligibleEvents > 0
              ? 100.0 * sqrt(efficiency.eventWeighted2D *
                             (1.0 - efficiency.eventWeighted2D) /
                             efficiency.eligibleEvents)
              : 0.0;
      const double runtime = chrono::duration<double>(
                                 chrono::high_resolution_clock::now() - t0)
                                 .count();
      cout << string(100, '=') << '\n'
           << Terminal::Accent("Configuration") << '\n';
      row("Detector", det->GetName());
      row("Detector ID", to_string(id));
      row("Detector type",
          det->GetPlanarPadConfig() ? "planar_pad" : "planar");
      row("Alignment", m_runAlignment ? "enabled" : "disabled");
      row("Alignment parameters", "dx, dy, dz, rotX, rotY, rotZ");
      row("Matching method", "track inside cluster envelope");
      row("Matching margin",
          fixedValue(m_efficiencyConfig.margin, 2, " mm"));
      row("Fiducial region X",
          "[" + fixedValue(m_effXMin, 1) + ", " +
              fixedValue(m_effXMax, 1) + "] mm");
      row("Fiducial region Y",
          "[" + fixedValue(m_effYMin, 1) + ", " +
              fixedValue(m_effYMax, 1) + "] mm");
      row("Efficiency bins",
          to_string(m_effXBins) + " x " + to_string(m_effYBins));
      row("Fake-efficiency estimation",
          m_efficiencyConfig.enableFakeEfficiency ? "enabled" : "disabled");
      separator();
      cout << Terminal::Accent("Input") << '\n';
      row("Reconstructed events", Terminal::Count(tracksPerEvent.size()));
      row("Reconstructed tracks", Terminal::Count(totalTrackEntries));
      row("Tracks inside fiducial region",
          Terminal::Count(efficiency.eligibleEvents));
      separator();
      cout << Terminal::Accent("DUT Response") << '\n';
      cout << "  " << left << setw(38) << "Quantity" << right << setw(19)
           << "X" << setw(20) << "Y" << '\n';
      const auto responseRow = [](const string& label, const string& x,
                                  const string& y) {
         cout << "  " << left << setw(38) << label << right << setw(19) << x
              << setw(20) << y << '\n';
      };
      responseRow("Reconstructed clusters",
                  Terminal::Count(selectedClusterSizesX.size()),
                  Terminal::Count(selectedClusterSizesY.size()));
      responseRow("Channel hits",
                  Terminal::Count(selectedHitADCsX.size()),
                  Terminal::Count(selectedHitADCsY.size()));
      responseRow("Mean cluster size",
                  fixedValue(mean(selectedClusterSizesX), 3, " strips"),
                  fixedValue(mean(selectedClusterSizesY), 3, " strips"));
      responseRow("Mean hit ADC",
                  fixedValue(mean(selectedHitADCsX), 2),
                  fixedValue(mean(selectedHitADCsY), 2));
      separator();
      cout << Terminal::Accent("Efficiency") << '\n';
      row("Matched tracks", Terminal::Count(efficiency.matchedEvents));
      row("Eligible tracks", Terminal::Count(efficiency.eligibleEvents));
      {
         ostringstream value;
         value << Terminal::Count(efficiency.matchedEvents) << " / "
               << Terminal::Count(efficiency.eligibleEvents) << " = "
               << fixed << setprecision(2)
               << 100.0 * efficiency.eventWeighted2D << " ± "
               << efficiencyError << '%';
         row("Overall efficiency", value.str());
      }
      row("Efficiency X",
          fixedValue(100.0 * efficiency.eventWeightedX, 2, "%"));
      row("Efficiency Y",
          fixedValue(100.0 * efficiency.eventWeightedY, 2, "%"));
      row("Fake-match probability",
          fixedValue(100.0 * efficiency.fake2D, 3, "%"));
      separator();
      cout << Terminal::Accent("Spatial Performance") << '\n'
           << "  " << left << setw(38) << "Quantity" << right << setw(19)
           << "X" << setw(20) << "Y" << '\n';
      const auto spatialRow = [](const string& label, double x, double y) {
         ostringstream xValue, yValue;
         xValue << fixed << setprecision(2) << 1000.0 * x << " um";
         yValue << fixed << setprecision(2) << 1000.0 * y << " um";
         cout << "  " << left << setw(38) << label << right << setw(19)
              << xValue.str() << setw(20) << yValue.str() << '\n';
      };
      spatialRow("Residual mean", resolutionX.mean, resolutionY.mean);
      spatialRow("Residual sigma68", resolutionX.sigma68,
                 resolutionY.sigma68);
      spatialRow("Residual RMS", resolutionX.rms, resolutionY.rms);
      separator();
      cout << Terminal::Accent("Status") << '\n';
      row(Terminal::Success("[PASS] DUT Analysis completed"), "");
      row("Runtime", fixedValue(runtime, 1, " s"));
      cout << string(100, '=') << '\n';
   }

   fDut->Close();
   delete fDut;

   f->Close();
   delete f;

   return true;
}

void RunDUTAlign(const std::vector<Event>& events,
                 std::shared_ptr<Detector> detector, int detID,
                 const DUTAlignmentRegion& region) {

   if (events.empty()) {
      std::cerr << "[DUT Alignment] No events to analyze for DUT " << detID << "!" << std::endl;
      return;
   }

   std::cout << "[DUT " << detID << "] Aligning (6-parameter)..." << std::endl;
   UInt_t nPar = 6;

   // The broad CC residual distribution can hide a coherent Y scale trend
   // from a generic chi2 minimizer.  Use the zero of the observable
   // CC-residual-vs-predY slope as the rotX seed; rotX remains free in the
   // subsequent staged fit.
   double prefittedRotX = 0.0;
   double bestAbsSlope = std::numeric_limits<double>::infinity();
   double previousRotX = -0.05;
   double previousSlope = std::numeric_limits<double>::quiet_NaN();
   double scanParameters[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
   for (int step = 0; step <= 100; ++step) {
      const double candidateRotX = -0.05 + 0.001 * step;
      scanParameters[3] = candidateRotX;
      const double slope = DUTCCResidualSlopeY(
          scanParameters, events, detector, detID, region, 2.0);
      if (!std::isfinite(slope)) continue;
      if (std::abs(slope) < bestAbsSlope) {
         bestAbsSlope = std::abs(slope);
         prefittedRotX = candidateRotX;
      }
      if (std::isfinite(previousSlope) && slope * previousSlope < 0.0) {
         const double denominator = slope - previousSlope;
         if (std::abs(denominator) > 1e-12) {
            const double root = previousRotX -
                previousSlope * (candidateRotX - previousRotX) /
                    denominator;
            if (std::abs(root) <= 0.05) prefittedRotX = root;
         }
         break;
      }
      previousRotX = candidateRotX;
      previousSlope = slope;
   }
   detector->SetAlignment(0, 0, 0, 0, 0, 0);

   // Deterministic annealing: start with broad, continuously varying
   // associations and tighten them only after the geometry has moved.  In the
   // final stage, reject events for which no cluster clearly dominates.
   const double matchingSigmas[] = {2.0, 1.0, 0.5};
   const double minMatchDominances[] = {0.0, 0.55, 0.75};
   double fitResult[6] = {0.0, 0.0, 0.0, prefittedRotX, 0.0, 0.0};
   double finalMinimum = std::numeric_limits<double>::quiet_NaN();
   for (size_t stage = 0; stage < 3; ++stage) {
      auto minimizer =
          ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
      minimizer->SetTolerance(1e-6);
      minimizer->SetPrecision(1e-10);
      minimizer->SetStrategy(2);
      minimizer->SetMaxFunctionCalls(10000);
      minimizer->SetMaxIterations(2000);
      minimizer->SetPrintLevel(0);

      auto chi2Func = [&events, &detector, detID, &region,
                       stage, &matchingSigmas,
                       &minMatchDominances](const double* par) {
         return DUTChi2Objective(par, events, detector, detID, region,
                                 matchingSigmas[stage],
                                 minMatchDominances[stage]);
      };
      ROOT::Math::Functor f(chi2Func, nPar);
      minimizer->SetFunction(f);
      minimizer->SetLimitedVariable(0, "dx", fitResult[0], 0.01, -5.0, 5.0);
      minimizer->SetLimitedVariable(1, "dy", fitResult[1], 0.01, -5.0, 5.0);
      // This beam has too little angular spread to determine z independently
      // from transverse translations.
      minimizer->SetFixedVariable(2, "dz", 0.0);
      minimizer->SetLimitedVariable(3, "rotX", fitResult[3], 1e-4,
                                    -0.05, 0.05);
      minimizer->SetLimitedVariable(4, "rotY", fitResult[4], 1e-3,
                                    -0.05, 0.05);
      minimizer->SetLimitedVariable(5, "rotZ", fitResult[5], 1e-3,
                                    -0.05, 0.05);

      minimizer->Minimize();
      const double* stageResult = minimizer->X();
      bool finiteResult = std::isfinite(minimizer->MinValue());
      for (UInt_t i = 0; i < nPar; ++i)
         finiteResult = finiteResult && std::isfinite(stageResult[i]);
      if (!finiteResult) {
         detector->SetAlignment(0, 0, 0, 0, 0, 0);
         std::cerr << "[DUT Alignment] Fit failed for DUT " << detID
                   << "; keeping configured geometry\n";
         delete minimizer;
         return;
      }
      std::copy(stageResult, stageResult + nPar, fitResult);
      finalMinimum = minimizer->MinValue();
      delete minimizer;
   }

   // 应用结果
   double dx = fitResult[0];
   double dy = fitResult[1];
   double dz = fitResult[2];
   double rotX = fitResult[3];
   double rotY = fitResult[4];
   double rotZ = fitResult[5];

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

   std::cout << "[DUT " << detID << "] Final chi2: " << finalMinimum << std::endl;
}

// ========== DUT对齐私有方法 ==========

LocalHit CalcuDutResidual(std::shared_ptr<Detector> detector,
                          const std::vector<Cluster>& clusters,
                          const TVector3& predL, double& residualX,
                          double& residualY, bool useCentroid) {

   double predX = predL.X();
   double predY = predL.Y();

   // 从探测器配置读取参数，消除魔法数字
   const auto* config = detector->GetPlanarConfig();
   const auto [typeX, typeY] =
       config ? PlanarAxisTypes(*config)
              : std::pair<int, int>{-1, -1};

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
      const auto& cluster = clusters[bestClusterXIndex];
      bestPosX = (useCentroid ? cluster.centroid : cluster.pos) * pitchX;
      residualX = bestPosX - predX;
   } else {
      bestPosX = DUTAnalysisConfig::kInvalidValue;
      residualX = DUTAnalysisConfig::kInvalidValue;
   }

   if (bestClusterYIndex != -1) {
      double pitchY = config->readoutPlanePitch.at(typeY);
      const auto& cluster = clusters[bestClusterYIndex];
      bestPosY = (useCentroid ? cluster.centroid : cluster.pos) * pitchY;
      residualY = bestPosY - predY;
   } else {
      bestPosY = DUTAnalysisConfig::kInvalidValue;
      residualY = DUTAnalysisConfig::kInvalidValue;
   }

   localHit.localPos.SetXYZ(bestPosX, bestPosY, 0);
   localHit.clusterIndices = {bestClusterXIndex, bestClusterYIndex};

   return localHit;
}

double DUTCCResidualSlopeY(const double* par,
                           const std::vector<Event>& events,
                           std::shared_ptr<Detector> detector, int detID,
                           const DUTAlignmentRegion& region,
                           double matchingSigma) {
   detector->SetAlignment(par[0], par[1], par[2], par[3], par[4], par[5]);

   double sumPredY = 0.0;
   double sumResidualY = 0.0;
   double sumPredY2 = 0.0;
   double sumPredYResidualY = 0.0;
   size_t count = 0;
   for (const auto& event : events) {
      const auto frame = event.detectorFramesMap.find(detID);
      if (frame == event.detectorFramesMap.end() ||
          frame->second->Clusters().empty())
         continue;
      const TVector3 predictedGlobal = detector->CalcHitFromTrack(event.track);
      const TVector3 predictedLocal = detector->GlobalToLocal(predictedGlobal);
      const double predX = predictedLocal.X();
      const double predY = predictedLocal.Y();
      if (!std::isfinite(predX) || !std::isfinite(predY) ||
          predX < region.xMin || predX > region.xMax ||
          predY < region.yMin || predY > region.yMax)
         continue;

      const DUTAlignmentMatch match = CalcuDutAlignmentMatch(
          detector, frame->second->Clusters(), predictedLocal,
          matchingSigma);
      const double residualX = match.residualX;
      const double residualY = match.residualY;
      if (!std::isfinite(residualX) || !std::isfinite(residualY) ||
          residualX == DUTAnalysisConfig::kInvalidValue ||
          residualY == DUTAnalysisConfig::kInvalidValue ||
          std::abs(residualX) > 1.0 || std::abs(residualY) > 3.0)
         continue;

      sumPredY += predY;
      sumResidualY += residualY;
      sumPredY2 += predY * predY;
      sumPredYResidualY += predY * residualY;
      ++count;
   }

   if (count < 100) return std::numeric_limits<double>::quiet_NaN();
   const double denominator =
       count * sumPredY2 - sumPredY * sumPredY;
   if (!(denominator > 1e-12) || !std::isfinite(denominator))
      return std::numeric_limits<double>::quiet_NaN();
   return (count * sumPredYResidualY - sumPredY * sumResidualY) /
          denominator;
}

double DUTChi2Objective(const double* par, const std::vector<Event>& events,
                        std::shared_ptr<Detector> detector, int detID,
                        const DUTAlignmentRegion& region,
                        double matchingSigma,
                        double minMatchDominance) {

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
   double sumPredY = 0.0;
   double sumResidualY = 0.0;
   double sumPredY2 = 0.0;
   double sumPredYResidualY = 0.0;

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

      if (!std::isfinite(predX) || !std::isfinite(predY)) continue;
      if (predX < region.xMin || predX > region.xMax ||
          predY < region.yMin || predY > region.yMax)
         continue;
      // Geometry alignment must not absorb uTPC timing distortions.  Match
      // against the charge centroid and reconsider all candidates whenever
      // the trial geometry changes.
      const DUTAlignmentMatch match = CalcuDutAlignmentMatch(
          detector, clusters, predL, matchingSigma);
      residualX = match.residualX;
      residualY = match.residualY;

      if (!std::isfinite(residualX) || !std::isfinite(residualY) ||
          residualX == DUTAnalysisConfig::kInvalidValue ||
          residualY == DUTAnalysisConfig::kInvalidValue ||
          match.dominanceX < minMatchDominance ||
          match.dominanceY < minMatchDominance)
         continue;

      double res = residualX * residualX + residualY * residualY;
      if (!std::isfinite(res) || std::abs(residualX) > 1.0 ||
          std::abs(residualY) > 3.0)
         continue;
      chi2 += res;
      sumPredY += predY;
      sumResidualY += residualY;
      sumPredY2 += predY * predY;
      sumPredYResidualY += predY * residualY;
      nEvents++;
   }

   if (nEvents < 100) return 1e9;
   const double slopeDenominator =
       nEvents * sumPredY2 - sumPredY * sumPredY;
   if (!(slopeDenominator > 1e-12) ||
       !std::isfinite(slopeDenominator))
      return 1e9;
   const double residualSlopeY =
       (nEvents * sumPredYResidualY - sumPredY * sumResidualY) /
       slopeDenominator;
   const double meanResidualY = sumResidualY / nEvents;

   // Keep the broad-residual chi2 for translations and cluster matching,
   // while preventing it from pulling rotX away from the slope-zero seed.
   constexpr double kResidualSlopePenalty = 1e4;
   constexpr double kResidualMeanPenalty = 100.0;
   return chi2 / nEvents +
          kResidualSlopePenalty * residualSlopeY * residualSlopeY +
          kResidualMeanPenalty * meanResidualY * meanResidualY;
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
   result.rms = rms;
   result.efficiency = double(N) / totalEvents;
   result.eff3S = eff3S;
   result.eff5S = eff5S;
   result.eff1mm = eff1mm;

   return result;
}

REGISTER_SCRIPT("DUTAnalysis", DUTAnalysisScript);
