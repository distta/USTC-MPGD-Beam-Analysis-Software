// Plot DUT time resolution and bad-pad-masked efficiency across runs.
//
// Run from the repository root:
//   root -l -b -q tools/PlotLHCb75umRunScan.C
//
// The first two table columns keep a 1050 V difference, corresponding to
// |E_drift| = 3500 V/cm. The horizontal axis is |V_top| ("High voltage").
// Outputs are written below <dataDirectory>/result.

#include <TCanvas.h>
#include <TEfficiency.h>
#include <TFile.h>
#include <TGaxis.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TPad.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

// ======================== USER CONFIGURATION ========================
const std::string dataDirectory = "data/LHCb";
const std::string detectorName = "75 #mum #muRGroove Pad readout";
const std::string gas = "Ar/CF4/ISO=88:10:2";
const int configuredDUTID = 10;
const double driftFieldVPerCm = 600.0;
const std::string configuredOutputDirectory = "auto";

struct ScanSetting {
   int voltage1;
   int voltage2;
   int runID;
};

// Only the explicitly listed settings are plotted.
const std::array<ScanSetting, 20> kScanSettings = {{
    {-800, -600, 149},
    {-800, -610, 151},
    {-800, -620, 153},
    {-800, -630, 155},
    {-800, -640, 157},
    {-800, -650, 159},
    {-800, -660, 161},
    {-800, -670, 163},
    {-800, -680, 165},
    {-800, -690, 167},
    {-800, -700, 169},
    {-800, -710, 171},
    {-800, -720, 173},
    {-800, -730, 175},
    {-800, -740, 177},
    {-800, -750, 179},
    {-800, -760, 181},
    {-800, -770, 183},
    {-800, -780, 185},
}};

// One-based (column, row) pad bins. These four bins remain near zero
// efficiency across the voltage scan and are therefore excluded explicitly.
// Edit this set if the channel mask changes.
const std::set<std::pair<int, int>> kBadPads = {
    {5, 9},
    {5, 5},
    {7, 4},
    {9, 3},
    {9, 7},
};
// ====================== END USER CONFIGURATION ======================

struct EfficiencyResult {
   bool valid{false};
   long long total{0};
   long long passed{0};
   long long excluded{0};
   double value{std::numeric_limits<double>::quiet_NaN()};
   double errorLow{std::numeric_limits<double>::quiet_NaN()};
   double errorHigh{std::numeric_limits<double>::quiet_NaN()};
};

struct TimingResult {
   bool valid{false};
   long long entries{0};
   bool timeWalkApplied{false};
   double resolution{std::numeric_limits<double>::quiet_NaN()};
   double error{std::numeric_limits<double>::quiet_NaN()};
};

std::string RunFile(const std::string& baseDirectory, int run,
                    const char* fileName) {
   return baseDirectory + "/result/" + std::to_string(run) + "/" + fileName;
}

bool FileExists(const std::string& path) {
   return gSystem->AccessPathName(path.c_str()) == kFALSE;
}

EfficiencyResult ReadEfficiency(const std::string& path, int dutID) {
   EfficiencyResult result;
   if (!FileExists(path)) return result;

   TFile file(path.c_str(), "READ");
   if (file.IsZombie()) return result;
   auto* tree = dynamic_cast<TTree*>(file.Get("PadDUTTree"));
   auto* map = dynamic_cast<TH2D*>(
       file.Get(("DUT_" + std::to_string(dutID) +
                 "/Efficiency/hEfficiency2D")
                    .c_str()));
   if (!tree || !map) return result;
   for (const char* branch : {"dutID", "predX", "predY", "hitFlag"}) {
      if (!tree->GetBranch(branch)) return result;
   }

   Int_t eventDUT = -1;
   Int_t hitFlag = 0;
   Double_t predictedX = 0.0;
   Double_t predictedY = 0.0;
   tree->SetBranchAddress("dutID", &eventDUT);
   tree->SetBranchAddress("predX", &predictedX);
   tree->SetBranchAddress("predY", &predictedY);
   tree->SetBranchAddress("hitFlag", &hitFlag);
   const TAxis* xAxis = map->GetXaxis();
   const TAxis* yAxis = map->GetYaxis();

   for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry) {
      tree->GetEntry(entry);
      if (eventDUT != dutID || !std::isfinite(predictedX) ||
          !std::isfinite(predictedY)) {
         continue;
      }
      const int xBin = xAxis->FindFixBin(predictedX);
      const int yBin = yAxis->FindFixBin(predictedY);
      if (xBin < 1 || xBin > xAxis->GetNbins() ||
          yBin < 1 || yBin > yAxis->GetNbins()) {
         continue;
      }
      if (kBadPads.count({xBin, yBin}) != 0) {
         ++result.excluded;
         continue;
      }
      ++result.total;
      if (hitFlag != 0) ++result.passed;
   }
   if (result.total == 0) return result;
   result.value =
       static_cast<double>(result.passed) / static_cast<double>(result.total);
   const double lower = TEfficiency::ClopperPearson(
       result.total, result.passed, 0.683, false);
   const double upper = TEfficiency::ClopperPearson(
       result.total, result.passed, 0.683, true);
   result.errorLow = result.value - lower;
   result.errorHigh = upper - result.value;
   result.valid = true;
   return result;
}

TimingResult ReadTiming(const std::string& path, int dutID) {
   TimingResult result;
   if (!FileExists(path)) return result;

   TFile file(path.c_str(), "READ");
   if (file.IsZombie()) return result;
   auto* results = dynamic_cast<TTree*>(file.Get("TimingResults"));
   if (results && results->GetBranch("analysis") &&
       results->GetBranch("object") &&
       results->GetBranch("detectorA") &&
       results->GetBranch("entries") &&
       results->GetBranch("resolutionNs")) {
      std::string* analysis = nullptr;
      std::string* object = nullptr;
      Int_t detectorA = -1;
      Long64_t entries = 0;
      Double_t resolution =
          std::numeric_limits<double>::quiet_NaN();
      Double_t error = 0.0;
      results->SetBranchAddress("analysis", &analysis);
      results->SetBranchAddress("object", &object);
      results->SetBranchAddress("detectorA", &detectorA);
      results->SetBranchAddress("entries", &entries);
      results->SetBranchAddress("resolutionNs", &resolution);
      if (results->GetBranch("sigmaErrorNs"))
         results->SetBranchAddress("sigmaErrorNs", &error);
      for (Long64_t entry = 0; entry < results->GetEntries(); ++entry) {
         results->GetEntry(entry);
         if (!analysis || !object ||
             *analysis != "without_external_t0" ||
             *object != "dut" || detectorA != dutID ||
             !std::isfinite(resolution)) {
            continue;
         }
         result.valid = true;
         result.entries = entries;
         result.resolution = resolution;
         result.error = std::isfinite(error) ? error : 0.0;
         return result;
      }
   }

   // Backward-compatible fallback for results produced before TimingResults.
   const std::string treePath =
       "DUTTimeResolution/DUT_" + std::to_string(dutID) + "/DUTResolution";
   auto* tree = dynamic_cast<TTree*>(file.Get(treePath.c_str()));
   if (!tree || tree->GetEntries() < 1 || !tree->GetBranch("valid") ||
       !tree->GetBranch("dutResolutionNs")) {
      return result;
   }

   Int_t valid = 0;
   Int_t timeWalkApplied = 0;
   Long64_t entries = 0;
   Double_t resolution = std::numeric_limits<double>::quiet_NaN();
   Double_t error = 0.0;
   tree->SetBranchAddress("valid", &valid);
   tree->SetBranchAddress("dutResolutionNs", &resolution);
   if (tree->GetBranch("dutResolutionErrorNs"))
      tree->SetBranchAddress("dutResolutionErrorNs", &error);
   if (tree->GetBranch("entries"))
      tree->SetBranchAddress("entries", &entries);
   if (tree->GetBranch("timeWalkApplied"))
      tree->SetBranchAddress("timeWalkApplied", &timeWalkApplied);
   tree->GetEntry(0);
   result.valid = valid != 0 && std::isfinite(resolution);
   result.entries = entries;
   result.timeWalkApplied = timeWalkApplied != 0;
   result.resolution = resolution;
   result.error = std::isfinite(error) ? error : 0.0;
   return result;
}

void StyleEfficiencyGraph(TGraphAsymmErrors& graph, int series) {
   const std::array<int, 2> colors = {kBlue + 1, kRed + 1};
   const std::array<int, 2> markers = {20, 21};
   graph.SetLineColor(colors[series]);
   graph.SetMarkerColor(colors[series]);
   graph.SetMarkerStyle(markers[series]);
   graph.SetMarkerSize(1.1);
   graph.SetLineWidth(2);
}

void StyleTimingGraph(TGraphErrors& graph, int series) {
   const std::array<int, 2> colors = {kBlue + 1, kRed + 1};
   const std::array<int, 2> markers = {20, 21};
   graph.SetLineColor(colors[series]);
   graph.SetMarkerColor(colors[series]);
   graph.SetMarkerStyle(markers[series]);
   graph.SetMarkerSize(1.1);
   graph.SetLineWidth(2);
}

}  // namespace

void PlotLHCb75umRunScan(
    const char* baseDirectoryArgument = "",
    const char* outputDirectoryArgument = "",
    int dutIDArgument = -1) {
   const std::string baseDirectory =
       baseDirectoryArgument && baseDirectoryArgument[0] != '\0'
           ? baseDirectoryArgument
           : dataDirectory;
   const int dutID =
       dutIDArgument >= 0 ? dutIDArgument : configuredDUTID;
   gStyle->SetOptStat(0);
   gStyle->SetTitleBorderSize(0);

   std::vector<ScanSetting> settings(kScanSettings.begin(),
                                     kScanSettings.end());
   std::sort(settings.begin(), settings.end(),
             [](const ScanSetting& first, const ScanSetting& second) {
                return std::abs(first.voltage2) <
                       std::abs(second.voltage2);
             });

   TGraphAsymmErrors efficiencyCurve;
   TGraphErrors timingCurve;
   efficiencyCurve.SetName("efficiency");
   timingCurve.SetName("timeResolution");

   std::cout << "Bad pad mask (one-based column,row):";
   for (const auto& [column, padRow] : kBadPads)
      std::cout << " (" << column << ',' << padRow << ')';
   std::cout << "\n\n";

   for (const ScanSetting& setting : settings) {
      const double voltage = std::abs(setting.voltage2);
      const int run = setting.runID;
      const EfficiencyResult efficiency = ReadEfficiency(
          RunFile(baseDirectory, run, "PadDUTInfo.root"), dutID);
      const TimingResult timing = ReadTiming(
          RunFile(baseDirectory, run, "TimeResolution.root"), dutID);

      if (efficiency.valid && timing.valid) {
         const int point = efficiencyCurve.GetN();
         efficiencyCurve.SetPoint(
             point, voltage, 100.0 * efficiency.value);
         efficiencyCurve.SetPointError(
             point, 0.0, 0.0, 100.0 * efficiency.errorLow,
             100.0 * efficiency.errorHigh);
         timingCurve.SetPoint(point, voltage, timing.resolution);
         timingCurve.SetPointError(point, 0.0, timing.error);
      }

      std::cout << "run " << std::setw(3) << run
                << "  V1=" << std::setw(5) << setting.voltage1
                << " V  V2=" << std::setw(5) << setting.voltage2 << " V";
      if (efficiency.valid) {
         std::cout << "  efficiency=" << std::fixed
                   << std::setprecision(3)
                   << 100.0 * efficiency.value << "%"
                   << " (" << efficiency.passed << '/'
                   << efficiency.total << ", excluded "
                   << efficiency.excluded << ')';
      } else {
         std::cout << "  efficiency=missing";
      }
      if (timing.valid) {
         std::cout << "  time=" << std::fixed
                   << std::setprecision(3) << timing.resolution
                   << " +/- " << timing.error << " ns";
         if (timing.timeWalkApplied) std::cout << " [time-walk corrected]";
      } else {
         std::cout << "  time=missing";
      }
      std::cout << '\n';
   }
   StyleEfficiencyGraph(efficiencyCurve, 0);
   StyleTimingGraph(timingCurve, 1);
   efficiencyCurve.SetLineColor(kBlack);
   efficiencyCurve.SetMarkerColor(kBlack);

   double minimumVoltage = std::numeric_limits<double>::infinity();
   double maximumVoltage = -std::numeric_limits<double>::infinity();
   double minimumEfficiency = std::numeric_limits<double>::infinity();
   double maximumTiming = 0.0;
   double minimumTiming = std::numeric_limits<double>::infinity();
   for (int point = 0; point < efficiencyCurve.GetN(); ++point) {
      double voltage = 0.0, efficiency = 0.0;
      double timingVoltage = 0.0, timing = 0.0;
      efficiencyCurve.GetPoint(point, voltage, efficiency);
      timingCurve.GetPoint(point, timingVoltage, timing);
      minimumVoltage = std::min(minimumVoltage, voltage);
      maximumVoltage = std::max(maximumVoltage, voltage);
      minimumEfficiency = std::min(minimumEfficiency, efficiency);
      minimumTiming = std::min(minimumTiming, timing);
      maximumTiming = std::max(maximumTiming, timing);
   }
   if (!std::isfinite(minimumVoltage) || !std::isfinite(maximumVoltage))
      throw std::runtime_error("No valid run results were found");
   const double voltageSpan =
       std::max(1.0, maximumVoltage - minimumVoltage);
   const double voltageMinimum = minimumVoltage - 0.08 * voltageSpan;
   const double voltageMaximum = maximumVoltage + 0.08 * voltageSpan;
   const double timingSpan = std::max(0.2, maximumTiming - minimumTiming);
   const double timingMinimum =
       std::max(0.0, minimumTiming - 0.20 * timingSpan);
   const double timingMaximum = maximumTiming + 0.30 * timingSpan;
   const double efficiencyMinimum =
       std::max(0.0, minimumEfficiency - 3.0);
   const std::string gasSuffix = gas.empty() ? "" : " - " + gas;

   TCanvas canvas("cLHCb75umRunScan", "LHCb 75 um run scan", 1100, 750);
   canvas.SetLeftMargin(0.12);
   canvas.SetRightMargin(0.14);
   canvas.SetTopMargin(0.12);
   canvas.SetBottomMargin(0.13);
   canvas.SetGridx();
   canvas.SetGridy();
   efficiencyCurve.SetTitle(
       Form("%s%s, |E_{drift}| = %.0f V/cm;High voltage [V];"
            "Efficiency [%%]",
            detectorName.c_str(), gasSuffix.c_str(), driftFieldVPerCm));
   efficiencyCurve.GetXaxis()->SetLimits(voltageMinimum, voltageMaximum);
   efficiencyCurve.SetMinimum(efficiencyMinimum);
   efficiencyCurve.SetMaximum(100.3);
   constexpr int axisFont = 42;
   constexpr double axisLabelSize = 0.038;
   constexpr double axisTitleSize = 0.043;
   efficiencyCurve.GetYaxis()->SetLabelFont(axisFont);
   efficiencyCurve.GetYaxis()->SetTitleFont(axisFont);
   efficiencyCurve.GetYaxis()->SetLabelSize(axisLabelSize);
   efficiencyCurve.GetYaxis()->SetTitleSize(axisTitleSize);
   efficiencyCurve.GetYaxis()->SetAxisColor(kBlack);
   efficiencyCurve.GetYaxis()->SetLabelColor(kBlack);
   efficiencyCurve.GetYaxis()->SetTitleColor(kBlack);
   efficiencyCurve.Draw("APL");

   const double efficiencyMaximum = 100.3;
   TGraphErrors scaledTiming(timingCurve.GetN());
   scaledTiming.SetName("timeResolutionScaled");
   StyleTimingGraph(scaledTiming, 1);
   for (int point = 0; point < timingCurve.GetN(); ++point) {
      double voltage = 0.0, timing = 0.0;
      timingCurve.GetPoint(point, voltage, timing);
      const double scaled =
          efficiencyMinimum +
          (timing - timingMinimum) /
              (timingMaximum - timingMinimum) *
              (efficiencyMaximum - efficiencyMinimum);
      const double scaledError =
          timingCurve.GetErrorY(point) /
          (timingMaximum - timingMinimum) *
          (efficiencyMaximum - efficiencyMinimum);
      scaledTiming.SetPoint(point, voltage, scaled);
      scaledTiming.SetPointError(point, 0.0, scaledError);
   }
   scaledTiming.Draw("PL SAME");

   TGaxis timingAxis(
       voltageMaximum, efficiencyMinimum,
       voltageMaximum, efficiencyMaximum,
       timingMinimum, timingMaximum, 510, "+L");
   timingAxis.SetName("timeResolutionAxis");
   timingAxis.SetTitle("Time resolution [ns]");
   timingAxis.SetLineColor(kRed + 1);
   timingAxis.SetLabelColor(kRed + 1);
   timingAxis.SetTitleColor(kRed + 1);
   timingAxis.SetLabelFont(axisFont);
   timingAxis.SetTitleFont(axisFont);
   timingAxis.SetLabelSize(axisLabelSize);
   timingAxis.SetTitleSize(axisTitleSize);
   timingAxis.SetLineWidth(1);
   timingAxis.SetTitleOffset(1.15);
   timingAxis.Draw();

   TLegend legend(0.16, 0.72, 0.43, 0.87);
   legend.SetBorderSize(0);
   legend.SetFillStyle(0);
   legend.AddEntry(&efficiencyCurve, "Efficiency", "lp");
   legend.AddEntry(&scaledTiming, "Time resolution", "lp");
   legend.Draw();

   canvas.Modified();
   canvas.Update();

   const std::string resultDirectory = baseDirectory + "/result";
   const std::string outputDirectory =
       outputDirectoryArgument && outputDirectoryArgument[0] != '\0'
           ? outputDirectoryArgument
           : (configuredOutputDirectory == "auto"
                  ? resultDirectory
                  : configuredOutputDirectory);
   gSystem->mkdir(outputDirectory.c_str(), kTRUE);
   const std::string prefix = outputDirectory + "/LHCb75umRunScan";
   const std::string rootPath = prefix + ".root";
   TFile output(rootPath.c_str(), "RECREATE");
   if (output.IsZombie())
      throw std::runtime_error("Cannot create " + rootPath);
   canvas.Write();
   output.Close();

   std::cout << "\n[PlotLHCb75umRunScan] Saved:\n  "
             << rootPath << '\n';
}
