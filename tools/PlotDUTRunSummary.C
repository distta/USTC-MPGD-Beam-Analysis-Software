// Compare DUT resolution and efficiency across runs.
//
// Edit the USER CONFIGURATION section below, then run:
// root -l -b -q tools/PlotDUTRunSummary.C

// Outputs are written below <dataDirectory>/result when that directory exists,
// otherwise directly below <dataDirectory>:
//   DUTRunSummary_<detector>_<gas>.root  (summary and fake-efficiency canvases)

#include <TAxis.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TPad.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct RunVoltage {
   int runID = 0;
   double voltage = 0.0;
};

// ======================== USER CONFIGURATION ========================
// DUTInfo.root is searched in:
//   <dataDirectory>/result/<runID>/DUTInfo.root
// and then in:
//   <dataDirectory>/<runID>/DUTInfo.root
const std::string dataDirectory = "data/Beam202607/large-area_sector4";
const std::string detectorName = "Large area #muRGroove";
const std::string gas = "Ar/ISO 98/2";
const int dutID = 3;

// Axis ranges. Keep auto=true for automatic ranges; set it to false and edit
// the corresponding minimum/maximum values to use a fixed range.
const bool autoVoltageRange = true;
const double voltageRangeMinimum = 295.0;
const double voltageRangeMaximum = 325.0;
const bool autoResolutionRange = true;
const double resolutionRangeMinimum = 0.0;    // um
const double resolutionRangeMaximum = 250.0;  // um
const bool autoEfficiencyRange = true;
const double efficiencyRangeMinimum = 0.0;    // percent
const double efficiencyRangeMaximum = 100.3;  // percent
const bool autoFakeEfficiencyRange = true;
const double fakeEfficiencyRangeMinimum = 0.0;  // percent
const double fakeEfficiencyRangeMaximum = 1.0;  // percent

// Use "auto" to write into <dataDirectory>/result when it exists.
const std::string configuredOutputDirectory = "auto";

// {run ID, voltage [V]}. Replace the example voltages with measured values.
const std::vector<RunVoltage> runVoltages = {
    {69, 300.0},
    {71, 310},
    {73, 320},
    {75, 330},
    {77, 340},
    {79, 350},
    {81, 360},
    {83, 370},
    {85, 380}};
// ====================== END USER CONFIGURATION ======================

struct ResolutionResult {
   double sigmaUm = 0.0;
   double errorUm = 0.0;
   bool valid = false;
};

struct RunResult {
   int runID = 0;
   double voltage = 0.0;
   ResolutionResult resolutionX;
   ResolutionResult resolutionY;
   double efficiencyX = 0.0;
   double efficiencyY = 0.0;
   double efficiency2D = 0.0;
   double efficiencyErrorLowX = 0.0;
   double efficiencyErrorHighX = 0.0;
   double efficiencyErrorLowY = 0.0;
   double efficiencyErrorHighY = 0.0;
   double efficiencyErrorLow2D = 0.0;
   double efficiencyErrorHigh2D = 0.0;
   double fakeEfficiencyX = 0.0;
   double fakeEfficiencyY = 0.0;
   double fakeEfficiency2D = 0.0;
   double fakeEfficiencyErrorLowX = 0.0;
   double fakeEfficiencyErrorHighX = 0.0;
   double fakeEfficiencyErrorLowY = 0.0;
   double fakeEfficiencyErrorHighY = 0.0;
   double fakeEfficiencyErrorLow2D = 0.0;
   double fakeEfficiencyErrorHigh2D = 0.0;
   double margin = 0.0;
};

std::string JoinPath(const std::string& left, const std::string& right) {
   if (left.empty()) return right;
   return left.back() == '/' ? left + right : left + "/" + right;
}

bool Exists(const std::string& path) {
   return gSystem->AccessPathName(path.c_str()) == kFALSE;
}

std::string ResolveInputFile(const std::string& dataDirectory, int runID) {
   const std::string run = std::to_string(runID);
   const std::string standard = JoinPath(
       dataDirectory, "result/" + run + "/DUTInfo.root");
   if (Exists(standard)) return standard;
   const std::string direct = JoinPath(dataDirectory, run + "/DUTInfo.root");
   if (Exists(direct)) return direct;
   return standard;
}

std::string SanitizeName(const std::string& value) {
   std::string result;
   for (const unsigned char character : value) {
      if (std::isalnum(character))
         result.push_back(character);
      else if (!result.empty() && result.back() != '_')
         result.push_back('_');
   }
   while (!result.empty() && result.back() == '_') result.pop_back();
   return result.empty() ? "gas" : result;
}

ResolutionResult ReadCoreResolution(TFile& input, const char* coordinate,
                                    int requestedDUT) {
   const std::string canvasName =
       Form("hRes%s_DUT%dDouble", coordinate, requestedDUT);
   auto* canvas = dynamic_cast<TCanvas*>(input.Get(canvasName.c_str()));
   if (!canvas) {
      throw std::runtime_error("Missing double-Gaussian canvas '" +
                               canvasName + "'");
   }
   auto* function = dynamic_cast<TF1*>(canvas->FindObject("fDouble"));
   if (!function || function->GetNpar() < 6) {
      throw std::runtime_error("Missing fDouble in canvas '" + canvasName + "'");
   }

   const double sigmaFirst = std::abs(function->GetParameter(2));
   const double sigmaSecond = std::abs(function->GetParameter(5));
   const int coreParameter = sigmaFirst <= sigmaSecond ? 2 : 5;
   const double sigma = std::abs(function->GetParameter(coreParameter));
   const double error = function->GetParError(coreParameter);
   if (!(sigma > 0.0) || !std::isfinite(sigma) ||
       !(error >= 0.0) || !std::isfinite(error)) {
      throw std::runtime_error("Invalid core Gaussian width in canvas '" +
                               canvasName + "'");
   }
   return {1000.0 * sigma, 1000.0 * error, true};
}

struct EfficiencyPoint {
   double value = 0.0;
   double errorLow = 0.0;
   double errorHigh = 0.0;
};

struct EfficiencySummary {
   EfficiencyPoint x;
   EfficiencyPoint y;
   EfficiencyPoint twoD;
   double margin = 0.0;
};

EfficiencySummary ReadEfficiencySummary(TFile& input,
                                        const std::string& treePath) {
   auto* tree = dynamic_cast<TTree*>(input.Get(treePath.c_str()));
   if (!tree || tree->GetEntries() != 1) {
      throw std::runtime_error(
          "Expected exactly one entry in summary tree '" + treePath + "'");
   }
   const std::vector<std::string> requiredBranches = {
       "parameterMm", "efficiencyX", "errorLowX",
       "errorHighX", "efficiencyY", "errorLowY", "errorHighY",
       "efficiency2D", "errorLow2D", "errorHigh2D"};
   for (const auto& branch : requiredBranches) {
      if (!tree->GetBranch(branch.c_str())) {
         throw std::runtime_error("Missing branch '" + branch +
                                  "' in tree '" + treePath + "'");
      }
   }

   Double_t margin = 0.0;
   Double_t efficiencyX = 0.0, errorLowX = 0.0, errorHighX = 0.0;
   Double_t efficiencyY = 0.0, errorLowY = 0.0, errorHighY = 0.0;
   Double_t efficiency2D = 0.0, errorLow2D = 0.0, errorHigh2D = 0.0;
   tree->SetBranchAddress("parameterMm", &margin);
   tree->SetBranchAddress("efficiencyX", &efficiencyX);
   tree->SetBranchAddress("errorLowX", &errorLowX);
   tree->SetBranchAddress("errorHighX", &errorHighX);
   tree->SetBranchAddress("efficiencyY", &efficiencyY);
   tree->SetBranchAddress("errorLowY", &errorLowY);
   tree->SetBranchAddress("errorHighY", &errorHighY);
   tree->SetBranchAddress("efficiency2D", &efficiency2D);
   tree->SetBranchAddress("errorLow2D", &errorLow2D);
   tree->SetBranchAddress("errorHigh2D", &errorHigh2D);
   tree->GetEntry(0);

   return {{100.0 * efficiencyX, 100.0 * errorLowX, 100.0 * errorHighX},
           {100.0 * efficiencyY, 100.0 * errorLowY, 100.0 * errorHighY},
           {100.0 * efficiency2D, 100.0 * errorLow2D, 100.0 * errorHigh2D},
           margin};
}

RunResult ReadRun(const std::string& fileName, const RunVoltage& runVoltage,
                  int requestedDUT) {
   TFile input(fileName.c_str(), "READ");
   if (input.IsZombie()) {
      throw std::runtime_error("Cannot open '" + fileName + "'");
   }
   const std::string dutPath = "DUT_" + std::to_string(requestedDUT);
   const EfficiencySummary efficiency = ReadEfficiencySummary(
       input, dutPath + "/EfficiencySummary");
   const EfficiencySummary fake = ReadEfficiencySummary(
       input, dutPath + "/Fake/FakeSummary");
   if (std::abs(efficiency.margin - fake.margin) > 1e-9) {
      throw std::runtime_error("Efficiency and fake summaries use different margins in '" +
                               fileName + "'");
   }

   RunResult result;
   result.runID = runVoltage.runID;
   result.voltage = runVoltage.voltage;
   result.resolutionX = ReadCoreResolution(input, "X", requestedDUT);
   result.resolutionY = ReadCoreResolution(input, "Y", requestedDUT);
   result.efficiencyX = efficiency.x.value;
   result.efficiencyY = efficiency.y.value;
   result.efficiency2D = efficiency.twoD.value;
   result.efficiencyErrorLowX = efficiency.x.errorLow;
   result.efficiencyErrorHighX = efficiency.x.errorHigh;
   result.efficiencyErrorLowY = efficiency.y.errorLow;
   result.efficiencyErrorHighY = efficiency.y.errorHigh;
   result.efficiencyErrorLow2D = efficiency.twoD.errorLow;
   result.efficiencyErrorHigh2D = efficiency.twoD.errorHigh;
   result.fakeEfficiencyX = fake.x.value;
   result.fakeEfficiencyY = fake.y.value;
   result.fakeEfficiency2D = fake.twoD.value;
   result.fakeEfficiencyErrorLowX = fake.x.errorLow;
   result.fakeEfficiencyErrorHighX = fake.x.errorHigh;
   result.fakeEfficiencyErrorLowY = fake.y.errorLow;
   result.fakeEfficiencyErrorHighY = fake.y.errorHigh;
   result.fakeEfficiencyErrorLow2D = fake.twoD.errorLow;
   result.fakeEfficiencyErrorHigh2D = fake.twoD.errorHigh;
   result.margin = efficiency.margin;
   return result;
}

void StyleGraph(TGraph& graph, int color, int marker) {
   graph.SetLineColor(color);
   graph.SetMarkerColor(color);
   graph.SetMarkerStyle(marker);
   graph.SetMarkerSize(1.1);
   graph.SetLineWidth(2);
   graph.SetLineStyle(1);
}

}  // namespace

void PlotDUTRunSummary() {
   if (dataDirectory.empty()) {
      throw std::runtime_error("dataDirectory must not be empty");
   }
   if (detectorName.empty()) {
      throw std::runtime_error("detectorName must not be empty");
   }
   if (gas.empty()) throw std::runtime_error("gas must not be empty");
   if (dutID < 0) throw std::runtime_error("dutID must be non-negative");
   if (runVoltages.empty()) {
      throw std::runtime_error("runVoltages must not be empty");
   }

   const std::string& baseDirectory = dataDirectory;
   const std::string& gasLabel = gas;
   std::vector<RunVoltage> sortedRunVoltages = runVoltages;
   std::sort(sortedRunVoltages.begin(), sortedRunVoltages.end(),
             [](const auto& left, const auto& right) {
                if (left.voltage != right.voltage)
                   return left.voltage < right.voltage;
                return left.runID < right.runID;
             });
   std::vector<RunResult> results;
   for (const auto& item : sortedRunVoltages) {
      const std::string inputFile = ResolveInputFile(baseDirectory, item.runID);
      results.push_back(ReadRun(inputFile, item, dutID));
   }
   if (results.empty()) throw std::runtime_error("No DUT runs were configured");
   const double commonMargin = results.front().margin;
   for (const auto& result : results) {
      if (std::abs(result.margin - commonMargin) > 1e-9) {
         throw std::runtime_error(
             "All runs must use the same cluster-envelope margin");
      }
   }

   std::printf("\n%-8s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n",
               "Run", "Voltage[V]", "Xres[um]", "Yres[um]",
               "Xeff[%]", "Yeff[%]", "2Deff[%]",
               "Xfake[%]", "Yfake[%]", "2Dfake[%]");
   for (const auto& result : results) {
      std::printf("%-8d %10.2f %10.2f %10.2f %10.3f %10.3f %10.3f %10.5f %10.5f %10.5f\n",
                  result.runID, result.voltage,
                  result.resolutionX.sigmaUm, result.resolutionY.sigmaUm,
                  result.efficiencyX, result.efficiencyY,
                  result.efficiency2D, result.fakeEfficiencyX,
                  result.fakeEfficiencyY, result.fakeEfficiency2D);
   }

   const int points = static_cast<int>(results.size());
   TGraphErrors resolutionX(points), resolutionY(points);
   TGraphAsymmErrors efficiencyX(points), efficiencyY(points), efficiency2D(points);
   TGraphAsymmErrors fakeEfficiencyX(points), fakeEfficiencyY(points),
       fakeEfficiency2D(points);
   resolutionX.SetName("resolutionX");
   resolutionY.SetName("resolutionY");
   efficiencyX.SetName("efficiencyX");
   efficiencyY.SetName("efficiencyY");
   efficiency2D.SetName("efficiency2D");
   fakeEfficiencyX.SetName("fakeEfficiencyX");
   fakeEfficiencyY.SetName("fakeEfficiencyY");
   fakeEfficiency2D.SetName("fakeEfficiency2D");

   double minimumVoltage = std::numeric_limits<double>::infinity();
   double maximumVoltage = -std::numeric_limits<double>::infinity();
   double minimumResolution = std::numeric_limits<double>::infinity();
   double maximumResolution = 0.0;
   double minimumEfficiency = 100.0;
   double maximumFakeEfficiency = 0.0;
   for (int index = 0; index < points; ++index) {
      const auto& result = results[index];
      resolutionX.SetPoint(index, result.voltage, result.resolutionX.sigmaUm);
      resolutionX.SetPointError(index, 0.0, result.resolutionX.errorUm);
      resolutionY.SetPoint(index, result.voltage, result.resolutionY.sigmaUm);
      resolutionY.SetPointError(index, 0.0, result.resolutionY.errorUm);
      efficiencyX.SetPoint(index, result.voltage, result.efficiencyX);
      efficiencyX.SetPointError(index, 0.0, 0.0,
                                result.efficiencyErrorLowX,
                                result.efficiencyErrorHighX);
      efficiencyY.SetPoint(index, result.voltage, result.efficiencyY);
      efficiencyY.SetPointError(index, 0.0, 0.0,
                                result.efficiencyErrorLowY,
                                result.efficiencyErrorHighY);
      efficiency2D.SetPoint(index, result.voltage, result.efficiency2D);
      efficiency2D.SetPointError(index, 0.0, 0.0,
                                 result.efficiencyErrorLow2D,
                                 result.efficiencyErrorHigh2D);
      fakeEfficiencyX.SetPoint(index, result.voltage,
                               result.fakeEfficiencyX);
      fakeEfficiencyX.SetPointError(index, 0.0, 0.0,
                                    result.fakeEfficiencyErrorLowX,
                                    result.fakeEfficiencyErrorHighX);
      fakeEfficiencyY.SetPoint(index, result.voltage,
                               result.fakeEfficiencyY);
      fakeEfficiencyY.SetPointError(index, 0.0, 0.0,
                                    result.fakeEfficiencyErrorLowY,
                                    result.fakeEfficiencyErrorHighY);
      fakeEfficiency2D.SetPoint(index, result.voltage,
                                result.fakeEfficiency2D);
      fakeEfficiency2D.SetPointError(index, 0.0, 0.0,
                                     result.fakeEfficiencyErrorLow2D,
                                     result.fakeEfficiencyErrorHigh2D);
      minimumVoltage = std::min(minimumVoltage, result.voltage);
      maximumVoltage = std::max(maximumVoltage, result.voltage);
      if (result.resolutionX.valid) {
         minimumResolution = std::min(minimumResolution, result.resolutionX.sigmaUm);
         maximumResolution = std::max(maximumResolution, result.resolutionX.sigmaUm);
      }
      if (result.resolutionY.valid) {
         minimumResolution = std::min(minimumResolution, result.resolutionY.sigmaUm);
         maximumResolution = std::max(maximumResolution, result.resolutionY.sigmaUm);
      }
      minimumEfficiency = std::min(
          minimumEfficiency,
          std::min(result.efficiency2D,
                   std::min(result.efficiencyX, result.efficiencyY)));
      maximumFakeEfficiency = std::max(
          maximumFakeEfficiency,
          std::max(result.fakeEfficiency2D,
                   std::max(result.fakeEfficiencyX,
                            result.fakeEfficiencyY)));
   }

   StyleGraph(resolutionX, kBlue + 1, 20);
   StyleGraph(resolutionY, kRed + 1, 21);
   StyleGraph(efficiencyX, kGreen + 2, 22);
   StyleGraph(efficiencyY, kOrange + 7, 33);
   StyleGraph(efficiency2D, kBlack, 23);
   StyleGraph(fakeEfficiencyX, kGreen + 2, 22);
   StyleGraph(fakeEfficiencyY, kOrange + 7, 33);
   StyleGraph(fakeEfficiency2D, kBlack, 23);

   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);
   gStyle->SetEndErrorSize(3);
   TCanvas canvas("cDUTRunSummary", "DUT run summary", 1100, 800);

   const double voltageSpan = maximumVoltage > minimumVoltage
                                  ? maximumVoltage - minimumVoltage
                                  : std::max(1.0, std::abs(minimumVoltage) * 0.05);
   const double voltagePadding = 0.08 * voltageSpan;

   if (!std::isfinite(minimumResolution) || !(maximumResolution > 0.0)) {
      throw std::runtime_error("No valid X/Y resolution fits were obtained");
   }
   const double resolutionSpan =
       std::max(5.0, maximumResolution - minimumResolution);
   const double automaticResolutionMinimum =
       std::max(0.0, minimumResolution - 0.20 * resolutionSpan);
   const double automaticResolutionMaximum =
       maximumResolution + 0.30 * resolutionSpan;
   const double automaticEfficiencyMinimum =
       std::max(0.0, minimumEfficiency - 2.0);
   const double automaticEfficiencyMaximum = 100.3;
   const double automaticFakeEfficiencyMinimum = 0.0;
   const double automaticFakeEfficiencyMaximum =
       std::max(1e-4, 1.25 * maximumFakeEfficiency);
   const double voltageMinimum = autoVoltageRange
                                     ? minimumVoltage - voltagePadding
                                     : voltageRangeMinimum;
   const double voltageMaximum = autoVoltageRange
                                     ? maximumVoltage + voltagePadding
                                     : voltageRangeMaximum;
   const double resolutionMinimum = autoResolutionRange
                                        ? automaticResolutionMinimum
                                        : resolutionRangeMinimum;
   const double resolutionMaximum = autoResolutionRange
                                        ? automaticResolutionMaximum
                                        : resolutionRangeMaximum;
   const double efficiencyMinimum = autoEfficiencyRange
                                        ? automaticEfficiencyMinimum
                                        : efficiencyRangeMinimum;
   const double efficiencyMaximum = autoEfficiencyRange
                                        ? automaticEfficiencyMaximum
                                        : efficiencyRangeMaximum;
   const double fakeEfficiencyMinimum = autoFakeEfficiencyRange
                                            ? automaticFakeEfficiencyMinimum
                                            : fakeEfficiencyRangeMinimum;
   const double fakeEfficiencyMaximum = autoFakeEfficiencyRange
                                            ? automaticFakeEfficiencyMaximum
                                            : fakeEfficiencyRangeMaximum;
   if (!(voltageMaximum > voltageMinimum) ||
       !(resolutionMaximum > resolutionMinimum) ||
       !(efficiencyMaximum > efficiencyMinimum) ||
       !(fakeEfficiencyMaximum > fakeEfficiencyMinimum)) {
      throw std::runtime_error("Each configured axis maximum must exceed its minimum");
   }

   TPad topPad("resolutionPad", "resolution", 0.0, 0.49, 1.0, 1.0);
   TPad bottomPad("efficiencyPad", "efficiency", 0.0, 0.0, 1.0, 0.51);
   topPad.SetLeftMargin(0.12);
   topPad.SetRightMargin(0.04);
   topPad.SetTopMargin(0.14);
   topPad.SetBottomMargin(0.02);
   bottomPad.SetLeftMargin(0.12);
   bottomPad.SetRightMargin(0.04);
   bottomPad.SetTopMargin(0.02);
   bottomPad.SetBottomMargin(0.16);
   topPad.SetGridx();
   topPad.SetGridy();
   bottomPad.SetGridx();
   bottomPad.SetGridy();
   topPad.Draw();
   bottomPad.Draw();

   topPad.cd();
   TMultiGraph resolutionGraph;
   resolutionGraph.SetTitle(
       Form("%s - %s;;Resolution [#mum]",
            detectorName.c_str(), gasLabel.c_str()));
   resolutionGraph.Add(&resolutionX, "LP");
   resolutionGraph.Add(&resolutionY, "LP");
   resolutionGraph.Draw("A");
   resolutionGraph.GetXaxis()->SetLimits(voltageMinimum, voltageMaximum);
   resolutionGraph.SetMinimum(resolutionMinimum);
   resolutionGraph.SetMaximum(resolutionMaximum);
   resolutionGraph.GetXaxis()->SetLabelSize(0.0);
   resolutionGraph.GetXaxis()->SetTickLength(0.0);
   resolutionGraph.GetYaxis()->SetTitleSize(0.065);
   resolutionGraph.GetYaxis()->SetLabelSize(0.055);
   resolutionGraph.GetYaxis()->SetTitleOffset(0.82);
   TLegend resolutionLegend(0.68, 0.68, 0.93, 0.86);
   resolutionLegend.SetBorderSize(0);
   resolutionLegend.SetFillStyle(0);
   resolutionLegend.AddEntry(&resolutionX, "X resolution", "lp");
   resolutionLegend.AddEntry(&resolutionY, "Y resolution", "lp");
   resolutionLegend.Draw();

   bottomPad.cd();
   TMultiGraph efficiencyGraph;
   efficiencyGraph.SetTitle(";Voltage [V];Efficiency [%]");
   efficiencyGraph.Add(&efficiencyX, "LP");
   efficiencyGraph.Add(&efficiencyY, "LP");
   efficiencyGraph.Add(&efficiency2D, "LP");
   efficiencyGraph.Draw("A");
   efficiencyGraph.GetXaxis()->SetLimits(voltageMinimum, voltageMaximum);
   efficiencyGraph.SetMinimum(efficiencyMinimum);
   efficiencyGraph.SetMaximum(efficiencyMaximum);
   efficiencyGraph.GetXaxis()->SetTitleSize(0.065);
   efficiencyGraph.GetXaxis()->SetLabelSize(0.055);
   efficiencyGraph.GetYaxis()->SetTitleSize(0.065);
   efficiencyGraph.GetYaxis()->SetLabelSize(0.055);
   efficiencyGraph.GetYaxis()->SetTitleOffset(0.82);
   TLegend efficiencyLegend(0.68, 0.66, 0.93, 0.89);
   efficiencyLegend.SetBorderSize(0);
   efficiencyLegend.SetFillStyle(0);
   efficiencyLegend.AddEntry(&efficiencyX, "X efficiency", "lp");
   efficiencyLegend.AddEntry(&efficiencyY, "Y efficiency", "lp");
   efficiencyLegend.AddEntry(&efficiency2D, "2D efficiency", "lp");
   efficiencyLegend.Draw();

   canvas.cd();
   canvas.Modified();
   canvas.Update();

   TCanvas fakeEfficiencyCanvas(
       "cDUTFakeEfficiency", "DUT fake-efficiency summary", 1100, 600);
   fakeEfficiencyCanvas.SetLeftMargin(0.12);
   fakeEfficiencyCanvas.SetRightMargin(0.04);
   fakeEfficiencyCanvas.SetTopMargin(0.12);
   fakeEfficiencyCanvas.SetBottomMargin(0.14);
   fakeEfficiencyCanvas.SetGridx();
   fakeEfficiencyCanvas.SetGridy();
   TMultiGraph fakeEfficiencyGraph;
   fakeEfficiencyGraph.SetTitle(
       Form("%s - %s;Voltage [V];Fake efficiency [%%]",
            detectorName.c_str(), gasLabel.c_str()));
   fakeEfficiencyGraph.Add(&fakeEfficiencyX, "LP");
   fakeEfficiencyGraph.Add(&fakeEfficiencyY, "LP");
   fakeEfficiencyGraph.Add(&fakeEfficiency2D, "LP");
   fakeEfficiencyGraph.Draw("A");
   fakeEfficiencyGraph.GetXaxis()->SetLimits(voltageMinimum, voltageMaximum);
   fakeEfficiencyGraph.SetMinimum(fakeEfficiencyMinimum);
   fakeEfficiencyGraph.SetMaximum(fakeEfficiencyMaximum);
   fakeEfficiencyGraph.GetXaxis()->SetTitleSize(0.05);
   fakeEfficiencyGraph.GetXaxis()->SetLabelSize(0.045);
   fakeEfficiencyGraph.GetYaxis()->SetTitleSize(0.05);
   fakeEfficiencyGraph.GetYaxis()->SetLabelSize(0.045);
   fakeEfficiencyGraph.GetYaxis()->SetTitleOffset(1.05);
   TLegend fakeEfficiencyLegend(0.68, 0.66, 0.93, 0.89);
   fakeEfficiencyLegend.SetBorderSize(0);
   fakeEfficiencyLegend.SetFillStyle(0);
   fakeEfficiencyLegend.AddEntry(&fakeEfficiencyX, "X fake efficiency", "lp");
   fakeEfficiencyLegend.AddEntry(&fakeEfficiencyY, "Y fake efficiency", "lp");
   fakeEfficiencyLegend.AddEntry(&fakeEfficiency2D, "2D fake efficiency", "lp");
   fakeEfficiencyLegend.Draw();

   fakeEfficiencyCanvas.Modified();
   fakeEfficiencyCanvas.Update();

   const std::string resultDirectoryCandidate = JoinPath(baseDirectory, "result");
   const std::string outputDirectory =
       configuredOutputDirectory.empty() || configuredOutputDirectory == "auto"
           ? (Exists(resultDirectoryCandidate) ? resultDirectoryCandidate
                                               : baseDirectory)
           : configuredOutputDirectory;
   gSystem->mkdir(outputDirectory.c_str(), kTRUE);
   const std::string outputBase = JoinPath(
       outputDirectory, "DUTRunSummary_" + SanitizeName(detectorName) + "_" +
                            SanitizeName(gasLabel));
   TFile output((outputBase + ".root").c_str(), "RECREATE");
   canvas.Write();
   fakeEfficiencyCanvas.Write();
   output.Close();

   std::printf("\n[PlotDUTRunSummary] Saved:\n  %s.root\n",
               outputBase.c_str());
}
