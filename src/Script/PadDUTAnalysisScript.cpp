#include "Script/PadDUTAnalysisScript.h"

#include "Detector/DetectorFactory.h"
#include "Detector/PlanarPad.h"
#include "Event/DetectorFrame.h"
#include "Script/Base/RawDataParser.h"
#include "Script/Base/ScriptFactory.h"
#include "Terminal.h"

#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include <TCanvas.h>
#include <TDirectory.h>
#include <TEfficiency.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <random>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

REGISTER_SCRIPT("PadDUTAnalysis", PadDUTAnalysisScript)

namespace {

constexpr double kInvalidValue = -999.0;
constexpr int kAxes = 3;
constexpr int kX = 0;
constexpr int kY = 1;
constexpr int k2D = 2;

struct PadMatch {
    int localHitIndex{-1};
    int clusterIndex{-1};
    TVector3 position{kInvalidValue, kInvalidValue, 0.0};
    double residualX{kInvalidValue};
    double residualY{kInvalidValue};
    double distance2{std::numeric_limits<double>::infinity()};

    bool IsValid() const { return clusterIndex >= 0; }
};

struct PadEnvelope {
    int clusterIndex{-1};
    double xLow{0.0};
    double xHigh{0.0};
    double yLow{0.0};
    double yHigh{0.0};
};

struct EfficiencySample {
    int eventID{0};
    int flatBin{-1};
    double predX{0.0};
    double predY{0.0};
    std::vector<PadEnvelope> envelopes;
};

struct FakeEfficiencySample {
    int flatBin{-1};
    const EfficiencySample* source{nullptr};
    double partnerPredX{0.0};
    double partnerPredY{0.0};
};

std::pair<double, double> DistributionRange(
    const std::vector<double>& values) {
    double minimum = std::numeric_limits<double>::infinity();
    double maximum = -std::numeric_limits<double>::infinity();
    for (double value : values) {
        if (!std::isfinite(value)) continue;
        minimum = std::min(minimum, value);
        maximum = std::max(maximum, value);
    }
    if (!std::isfinite(minimum) || !std::isfinite(maximum))
        return {0.0, 1.0};

    double low = std::min(0.0, minimum);
    double high = std::max(0.0, maximum);
    if (high <= low) return {low - 0.5, high + 0.5};
    const double padding = 0.05 * (high - low);
    if (low < 0.0) low -= padding;
    high += padding;
    return {low, high};
}

void StyleDistribution(TH1D& histogram) {
    histogram.SetDirectory(nullptr);
    histogram.SetStats(false);
    histogram.SetLineColor(kBlue + 1);
    histogram.SetLineWidth(2);
}

void WriteDUTHitProperties(
    TDirectory* detectorDirectory, int dutID,
    const std::vector<double>& clusterCharges,
    const std::vector<int>& clusterSizes,
    const std::vector<double>& hitADCs) {
    if (!detectorDirectory) return;
    auto* directory = detectorDirectory->mkdir("HitProperties");
    if (!directory) return;
    TDirectory::TContext context(directory);

    const auto chargeRange = DistributionRange(clusterCharges);
    TH1D charge(
        "hClusterCharge",
        ("DUT " + std::to_string(dutID) +
         " track-matched cluster charge;Cluster charge [ADC];Entries").c_str(),
        100, chargeRange.first, chargeRange.second);
    StyleDistribution(charge);
    for (double value : clusterCharges) charge.Fill(value);
    TCanvas chargeCanvas(
        "cClusterCharge", "Track-matched DUT cluster charge", 900, 700);
    charge.Draw("HIST");
    chargeCanvas.Write();

    int maximumClusterSize = 1;
    for (int size : clusterSizes)
        maximumClusterSize = std::max(maximumClusterSize, size);
    TH1D size(
        "hClusterSize",
        ("DUT " + std::to_string(dutID) +
         " track-matched cluster size;Cluster size [channels];Entries").c_str(),
        maximumClusterSize, 0.5, maximumClusterSize + 0.5);
    StyleDistribution(size);
    for (int value : clusterSizes) size.Fill(value);
    TCanvas sizeCanvas(
        "cClusterSize", "Track-matched DUT cluster size", 900, 700);
    size.Draw("HIST");
    sizeCanvas.Write();

    const auto adcRange = DistributionRange(hitADCs);
    TH1D adc(
        "hHitADC",
        ("DUT " + std::to_string(dutID) +
         " track-matched hit amplitude;Hit amplitude [ADC];Entries").c_str(),
        100, adcRange.first, adcRange.second);
    StyleDistribution(adc);
    for (double value : hitADCs) adc.Fill(value);
    TCanvas adcCanvas(
        "cHitADC", "Track-matched DUT hit amplitude", 900, 700);
    adc.Draw("HIST");
    adcCanvas.Write();
}

struct EfficiencyGrid {
    std::vector<long long> denominator;
    std::array<std::vector<long long>, kAxes> numerator;

    EfficiencyGrid(int xBins, int yBins)
        : denominator(static_cast<size_t>(xBins * yBins), 0) {
        for (auto& values : numerator) values.assign(denominator.size(), 0);
    }

    void Fill(int bin, const std::array<bool, kAxes>& matched) {
        if (bin < 0 || static_cast<size_t>(bin) >= denominator.size()) return;
        ++denominator[static_cast<size_t>(bin)];
        for (int axis = 0; axis < kAxes; ++axis) {
            if (matched[axis]) ++numerator[axis][static_cast<size_t>(bin)];
        }
    }
};

struct EfficiencyEstimate {
    long long total{0};
    long long passed{0};
    double value{0.0};
    double errorLow{0.0};
    double errorHigh{0.0};
};

struct ResidualResult {
    long long entries{0};
    double mean{0.0};
    double rms{0.0};
    double sigma68{0.0};
};

struct PadAlignmentQAPoint {
    double predX{0.0};
    double predY{0.0};
    double resX{0.0};
    double resY{0.0};
};

struct AlignmentMetrics {
    long long candidates{0};
    long long matches{0};
    double meanX{0.0};
    double meanY{0.0};
    double rmsX{0.0};
    double rmsY{0.0};
    std::vector<double> residualsX;
    std::vector<double> residualsY;
};

struct AlignmentSample {
    size_t eventIndex{0};
    TVector3 hitPosition;
};

struct AlignmentResult {
    bool attempted{false};
    bool applied{false};
    int passes{0};
    std::string status{"disabled"};
    double margin{0.0};
    double huberDelta{0.0};
    double objectiveBefore{0.0};
    double objectiveAfter{0.0};
    TVector3 correctionPosition;
    TVector3 correctionRotation;
    TVector3 finalPosition;
    TVector3 finalRotation;
    AlignmentMetrics before;
    AlignmentMetrics after;
};

bool Contains(const std::vector<int>& values, int value) {
    return std::find(values.begin(), values.end(), value) != values.end();
}

double DistanceToInterval(double coordinate, double low, double high) {
    if (coordinate < low) return low - coordinate;
    if (coordinate > high) return coordinate - high;
    return 0.0;
}

bool MatchesEnvelope(const PadEnvelope& envelope,
                     const TVector3& predicted, double margin) {
    return DistanceToInterval(
               predicted.X(), envelope.xLow, envelope.xHigh) <= margin &&
           DistanceToInterval(
               predicted.Y(), envelope.yLow, envelope.yHigh) <= margin;
}

std::array<bool, kAxes> MatchEnvelopes(
    const EfficiencySample& sample, double margin) {
    bool matchedX = false;
    bool matchedY = false;
    bool matched2D = false;
    for (const auto& envelope : sample.envelopes) {
        const bool x = DistanceToInterval(
                           sample.predX, envelope.xLow, envelope.xHigh) <= margin;
        const bool y = DistanceToInterval(
                           sample.predY, envelope.yLow, envelope.yHigh) <= margin;
        matchedX = matchedX || x;
        matchedY = matchedY || y;
        matched2D = matched2D || MatchesEnvelope(
                                      envelope,
                                      {sample.predX, sample.predY, 0.0},
                                      margin);
    }
    return {matchedX, matchedY, matched2D};
}

std::array<bool, kAxes> MatchFakeEnvelopes(
    const FakeEfficiencySample& sample, double margin) {
    if (!sample.source) return {false, false, false};
    bool matchedX = false;
    bool matchedY = false;
    bool matched2D = false;
    for (const auto& envelope : sample.source->envelopes) {
        const bool sourceMatched =
            DistanceToInterval(sample.source->predX,
                               envelope.xLow, envelope.xHigh) <= margin &&
            DistanceToInterval(sample.source->predY,
                               envelope.yLow, envelope.yHigh) <= margin;
        if (sourceMatched) continue;

        const bool x = DistanceToInterval(
                           sample.partnerPredX,
                           envelope.xLow, envelope.xHigh) <= margin;
        const bool y = DistanceToInterval(
                           sample.partnerPredY,
                           envelope.yLow, envelope.yHigh) <= margin;
        matchedX = matchedX || x;
        matchedY = matchedY || y;
        matched2D = matched2D || (x && y);
    }
    return {matchedX, matchedY, matched2D};
}

EfficiencyEstimate EstimateEfficiency(
    const EfficiencyGrid& grid, int axis, int minimumEntries) {
    EfficiencyEstimate result;
    for (size_t bin = 0; bin < grid.denominator.size(); ++bin) {
        if (grid.denominator[bin] < minimumEntries) continue;
        result.total += grid.denominator[bin];
        result.passed += grid.numerator[axis][bin];
    }
    if (result.total <= 0) return result;
    result.value = static_cast<double>(result.passed) / result.total;
    const double lower = TEfficiency::ClopperPearson(
        result.total, result.passed, 0.683, false);
    const double upper = TEfficiency::ClopperPearson(
        result.total, result.passed, 0.683, true);
    result.errorLow = result.value - lower;
    result.errorHigh = upper - result.value;
    return result;
}

int FindFlatBin(double x, double y,
                double xMin, double xMax, int xBins,
                double yMin, double yMax, int yBins,
                const std::vector<int>& excludedX,
                const std::vector<int>& excludedY) {
    if (x < xMin || x >= xMax || y < yMin || y >= yMax) return -1;
    const int xBin = static_cast<int>((x - xMin) / (xMax - xMin) * xBins) + 1;
    const int yBin = static_cast<int>((y - yMin) / (yMax - yMin) * yBins) + 1;
    if (xBin < 1 || xBin > xBins || yBin < 1 || yBin > yBins ||
        Contains(excludedX, xBin) || Contains(excludedY, yBin)) {
        return -1;
    }
    return (xBin - 1) * yBins + (yBin - 1);
}

std::vector<PadEnvelope> BuildPadEnvelopes(
    const DetectorFrame& frame, const planarPadConfig& config) {
    std::vector<PadEnvelope> envelopes;
    const auto& hits = frame.ChannelHits();
    const auto& clusters = frame.Clusters();
    const double halfSizeX = 0.5 * config.sizeX;
    const double halfSizeY = 0.5 * config.sizeY;

    for (size_t clusterIndex = 0; clusterIndex < clusters.size(); ++clusterIndex) {
        const auto& cluster = clusters[clusterIndex];
        if (cluster.type != config.planeType || cluster.channelHitIndices.empty()) {
            continue;
        }
        PadEnvelope envelope;
        envelope.clusterIndex = static_cast<int>(clusterIndex);
        envelope.xLow = std::numeric_limits<double>::infinity();
        envelope.xHigh = -std::numeric_limits<double>::infinity();
        envelope.yLow = std::numeric_limits<double>::infinity();
        envelope.yHigh = -std::numeric_limits<double>::infinity();
        for (int hitIndex : cluster.channelHitIndices) {
            if (hitIndex < 0 || static_cast<size_t>(hitIndex) >= hits.size()) continue;
            const auto& hit = hits[static_cast<size_t>(hitIndex)];
            if (!hit.HasID1() || hit.id0 < 0 || hit.id0 >= config.columns ||
                hit.id1 < 0 || hit.id1 >= config.rows) {
                continue;
            }
            const double centerX = hit.id0 * config.pitchX;
            const double centerY = hit.id1 * config.pitchY;
            envelope.xLow = std::min(
                envelope.xLow, centerX - halfSizeX);
            envelope.xHigh = std::max(envelope.xHigh, centerX + halfSizeX);
            envelope.yLow = std::min(
                envelope.yLow, centerY - halfSizeY);
            envelope.yHigh = std::max(envelope.yHigh, centerY + halfSizeY);
        }
        if (envelope.xLow <= envelope.xHigh &&
            envelope.yLow <= envelope.yHigh) {
            envelopes.push_back(envelope);
        }
    }
    return envelopes;
}

PadMatch FindMatchedPadHit(const DetectorFrame& frame,
                           const TVector3& predicted,
                           const planarPadConfig& config,
                           double margin) {
    PadMatch best;
    const auto envelopes = BuildPadEnvelopes(frame, config);
    const auto& localHits = frame.LocalHits();
    for (size_t index = 0; index < localHits.size(); ++index) {
        const auto& hit = localHits[index];
        if (hit.clusterIndices.empty()) continue;
        const int clusterIndex = hit.clusterIndices.front();
        const auto envelope = std::find_if(
            envelopes.begin(), envelopes.end(),
            [clusterIndex](const PadEnvelope& candidate) {
                return candidate.clusterIndex == clusterIndex;
            });
        if (envelope == envelopes.end() ||
            !MatchesEnvelope(*envelope, predicted, margin)) {
            continue;
        }

        const double residualX = hit.localPos.X() - predicted.X();
        const double residualY = hit.localPos.Y() - predicted.Y();
        const double distance2 = residualX * residualX + residualY * residualY;
        if (distance2 >= best.distance2) continue;
        best.localHitIndex = static_cast<int>(index);
        best.clusterIndex = clusterIndex;
        best.position = hit.localPos;
        best.residualX = residualX;
        best.residualY = residualY;
        best.distance2 = distance2;
    }
    return best;
}

double Mean(const std::vector<double>& values) {
    if (values.empty()) return 0.0;
    return std::accumulate(values.begin(), values.end(), 0.0) / values.size();
}

ResidualResult WriteResidual(const std::vector<double>& input,
                             TDirectory* directory,
                             const std::string& name) {
    ResidualResult result;
    std::vector<double> values;
    values.reserve(input.size());
    for (double value : input) {
        if (std::isfinite(value) && value != kInvalidValue) values.push_back(value);
    }
    if (values.empty() || !directory) return result;

    result.entries = static_cast<long long>(values.size());
    result.mean = Mean(values);
    double variance = 0.0;
    for (double value : values) {
        const double delta = value - result.mean;
        variance += delta * delta;
    }
    result.rms = std::sqrt(variance / values.size());

    std::sort(values.begin(), values.end());
    const auto quantile = [&](double probability) {
        const size_t index = static_cast<size_t>(
            std::llround(probability * (values.size() - 1)));
        return values[index];
    };
    result.sigma68 = 0.5 * (quantile(0.84) - quantile(0.16));
    double minimum = result.mean - 5.0 * result.rms;
    double maximum = result.mean + 5.0 * result.rms;
    if (!(maximum > minimum) || !std::isfinite(minimum) ||
        !std::isfinite(maximum)) {
        minimum = result.mean - 1.0;
        maximum = result.mean + 1.0;
    }

    TDirectory::TContext context(directory);
    TH1D histogram(name.c_str(),
                   (name + ";Residual [mm];Events").c_str(),
                   100, minimum, maximum);
    for (double value : values) histogram.Fill(value);
    histogram.SetStats(false);
    TCanvas canvas(("c" + name.substr(1)).c_str(),
                   (name + " residual distribution").c_str(), 900, 700);
    histogram.Draw("HIST");
    canvas.Write();
    return result;
}

double RobustResidualRange(const std::vector<PadAlignmentQAPoint>& points,
                           bool useX) {
    std::vector<double> absoluteResiduals;
    absoluteResiduals.reserve(points.size());
    for (const auto& point : points) {
        absoluteResiduals.push_back(
            std::abs(useX ? point.resX : point.resY));
    }
    if (absoluteResiduals.empty()) return 1.0;
    const size_t index = std::min(
        absoluteResiduals.size() - 1,
        static_cast<size_t>(0.99 * absoluteResiduals.size()));
    std::nth_element(absoluteResiduals.begin(),
                     absoluteResiduals.begin() + index,
                     absoluteResiduals.end());
    return std::clamp(1.2 * absoluteResiduals[index], 1.0, 100.0);
}

void WriteResidualMaps(TDirectory* directory,
                       const std::vector<PadAlignmentQAPoint>& points,
                       double predXMin, double predXMax,
                       double predYMin, double predYMax) {
    if (!directory || points.empty()) return;
    TDirectory::TContext context(directory);
    const double residualRangeX = RobustResidualRange(points, true);
    const double residualRangeY = RobustResidualRange(points, false);

    TH2D residualXVsHitX(
        "hResidualXVsHitX",
        "X residual vs track hit;Track-predicted Hit X [mm];Residual X [mm]",
        160, predXMin, predXMax, 200, -residualRangeX, residualRangeX);
    TH2D residualYVsHitY(
        "hResidualYVsHitY",
        "Y residual vs track hit;Track-predicted Hit Y [mm];Residual Y [mm]",
        160, predYMin, predYMax, 200, -residualRangeY, residualRangeY);
    TH2D residualXVsHitY(
        "hResidualXVsHitY",
        "X residual vs orthogonal track hit;Track-predicted Hit Y [mm];Residual X [mm]",
        160, predYMin, predYMax, 200, -residualRangeX, residualRangeX);
    TH2D residualYVsHitX(
        "hResidualYVsHitX",
        "Y residual vs orthogonal track hit;Track-predicted Hit X [mm];Residual Y [mm]",
        160, predXMin, predXMax, 200, -residualRangeY, residualRangeY);
    TH2D residual2D(
        "hResidual2D",
        "Two-dimensional residual correlation;Residual X [mm];Residual Y [mm]",
        200, -residualRangeX, residualRangeX,
        200, -residualRangeY, residualRangeY);

    for (const auto& point : points) {
        residualXVsHitX.Fill(point.predX, point.resX);
        residualYVsHitY.Fill(point.predY, point.resY);
        residualXVsHitY.Fill(point.predY, point.resX);
        residualYVsHitX.Fill(point.predX, point.resY);
        residual2D.Fill(point.resX, point.resY);
    }
    for (TH2D* histogram : {
             &residualXVsHitX, &residualYVsHitY,
             &residualXVsHitY, &residualYVsHitX, &residual2D}) {
        histogram->SetStats(false);
        histogram->SetOption("COLZ TEXT");
        histogram->Write();
    }
}

EfficiencyGrid BuildEfficiencyGrid(
    const std::vector<EfficiencySample>& samples,
    int xBins, int yBins, double margin) {
    EfficiencyGrid grid(xBins, yBins);
    for (const auto& sample : samples) {
        grid.Fill(sample.flatBin, MatchEnvelopes(sample, margin));
    }
    return grid;
}

EfficiencyGrid BuildFakeEfficiencyGrid(
    const std::vector<FakeEfficiencySample>& samples,
    int xBins, int yBins, double margin) {
    EfficiencyGrid grid(xBins, yBins);
    for (const auto& sample : samples) {
        grid.Fill(sample.flatBin, MatchFakeEnvelopes(sample, margin));
    }
    return grid;
}

TDirectory* GetOrCreateDirectory(TDirectory* parent, const std::string& name) {
    if (auto* directory = parent->GetDirectory(name.c_str())) return directory;
    return parent->mkdir(name.c_str());
}

void WriteEfficiencyMaps(TDirectory* directory,
                         const EfficiencyGrid& grid,
                         int xBins, double xMin, double xMax,
                         int yBins, double yMin, double yMax,
                         int minimumEntries,
                         const std::string& prefix) {
    if (!directory) return;
    TDirectory::TContext context(directory);
    std::array<std::unique_ptr<TH2D>, kAxes> maps;
    const std::array<const char*, kAxes> labels = {"X", "Y", "2D"};
    const bool isFake = prefix.find("Fake") != std::string::npos;
    for (int axis = 0; axis < kAxes; ++axis) {
        maps[axis] = std::make_unique<TH2D>(
            (prefix + "Efficiency" + labels[axis]).c_str(),
            (std::string(isFake ? "Pad fake efficiency " : "Pad efficiency ") +
             labels[axis] +
             ";Local X [mm];Local Y [mm];Efficiency").c_str(),
            xBins, xMin, xMax, yBins, yMin, yMax);
    }
    for (int xBin = 1; xBin <= xBins; ++xBin) {
        for (int yBin = 1; yBin <= yBins; ++yBin) {
            const int flatBin = (xBin - 1) * yBins + (yBin - 1);
            const long long total = grid.denominator[static_cast<size_t>(flatBin)];
            for (int axis = 0; axis < kAxes; ++axis) {
                const double value = total >= minimumEntries
                                         ? static_cast<double>(grid.numerator[axis][static_cast<size_t>(flatBin)]) / total
                                         : -1.0;
                maps[axis]->SetBinContent(xBin, yBin, value);
            }
        }
    }

    for (int axis = 0; axis < kAxes; ++axis) {
        auto& map = maps[axis];
        map->SetMinimum(0.0);
        map->SetMaximum(isFake ? 0.05 : 1.0);
        map->SetStats(false);
        map->SetMarkerSize(1.0);
        map->SetOption("COLZ TEXT");
        map->Write();
    }
}

void WriteEfficiencyProjections(TDirectory* directory,
                                const EfficiencyGrid& grid,
                                int xBins, double xMin, double xMax,
                                int yBins, double yMin, double yMax,
                                int minimumEntries,
                                const std::string& prefix) {
    if (!directory) return;
    TDirectory::TContext context(directory);
    const bool isFake = !prefix.empty();
    const std::string quantity = isFake ? "fake efficiency" : "efficiency";
    const auto overall = EstimateEfficiency(grid, k2D, minimumEntries);
    for (int projectionAxis : {kX, kY}) {
        const bool projectX = projectionAxis == kX;
        const int bins = projectX ? xBins : yBins;
        const int orthogonalBins = projectX ? yBins : xBins;
        const double minimum = projectX ? xMin : yMin;
        const double maximum = projectX ? xMax : yMax;
        const double binWidth = (maximum - minimum) / bins;
        const std::string label = projectX ? "X" : "Y";

        TGraphAsymmErrors graph;
        graph.SetName(("g" + prefix + "Efficiency2DVs" + label).c_str());
        graph.SetTitle(("2D " + quantity + " vs local " + label +
                        ";Local " + label + " [mm];2D " + quantity).c_str());
        graph.SetMarkerStyle(20);
        graph.SetMarkerSize(0.8);
        graph.SetLineWidth(2);
        int point = 0;
        for (int projectedBin = 0; projectedBin < bins; ++projectedBin) {
            long long total = 0;
            long long passed = 0;
            for (int orthogonalBin = 0; orthogonalBin < orthogonalBins;
                 ++orthogonalBin) {
                const int xBin = projectX ? projectedBin : orthogonalBin;
                const int yBin = projectX ? orthogonalBin : projectedBin;
                const int flatBin = xBin * yBins + yBin;
                if (grid.denominator[static_cast<size_t>(flatBin)] <
                    minimumEntries) {
                    continue;
                }
                total += grid.denominator[static_cast<size_t>(flatBin)];
                passed += grid.numerator[k2D][static_cast<size_t>(flatBin)];
            }
            if (total == 0) continue;
            const double efficiency = static_cast<double>(passed) / total;
            const double lower = TEfficiency::ClopperPearson(
                total, passed, 0.683, false);
            const double upper = TEfficiency::ClopperPearson(
                total, passed, 0.683, true);
            graph.SetPoint(point, minimum + (projectedBin + 0.5) * binWidth,
                           efficiency);
            graph.SetPointError(point, 0.5 * binWidth, 0.5 * binWidth,
                                efficiency - lower, upper - efficiency);
            ++point;
        }

        TCanvas canvas(("c" + prefix + "Efficiency2DVs" + label).c_str(),
                       ("2D " + quantity + " projection on local " +
                        label).c_str(),
                       900, 700);
        graph.SetMinimum(0.0);
        graph.SetMaximum(isFake ? 0.05 : 1.02);
        graph.Draw("APL");
        TLine average(minimum, overall.value, maximum, overall.value);
        average.SetLineColor(kRed + 1);
        average.SetLineStyle(2);
        average.SetLineWidth(2);
        average.Draw();
        TLegend legend(0.58, 0.18, 0.88, 0.31);
        legend.SetBorderSize(0);
        legend.AddEntry(&graph, ("2D " + quantity).c_str(), "lp");
        legend.AddEntry(&average,
                        Form("Average: %.2f%%", 100.0 * overall.value), "l");
        legend.Draw();
        canvas.Write();
    }
}

template <typename GridBuilder>
void WriteMarginScan(TDirectory* directory, GridBuilder&& buildGrid,
                     int minimumEntries,
                     double minimum, double maximum, double step,
                     const std::string& prefix) {
    if (!directory || step <= 0.0 || maximum < minimum) return;
    TDirectory::TContext context(directory);
    const bool isFake = !prefix.empty();
    const std::string quantity = isFake ? "fake efficiency" : "efficiency";
    std::array<TGraph, kAxes> graphs;
    const std::array<const char*, kAxes> labels = {"X", "Y", "2D"};
    for (int axis = 0; axis < kAxes; ++axis) {
        graphs[axis].SetName(
            ("g" + prefix + "EfficiencyVsMargin" + labels[axis]).c_str());
        graphs[axis].SetTitle(
            (std::string("Pad ") + labels[axis] + " " + quantity +
             " vs margin;Margin [mm];" + quantity).c_str());
    }
    int point = 0;
    for (double margin = minimum; margin <= maximum + 0.5 * step;
         margin += step, ++point) {
        const auto grid = buildGrid(margin);
        for (int axis = 0; axis < kAxes; ++axis) {
            graphs[axis].SetPoint(
                point, margin,
                EstimateEfficiency(grid, axis, minimumEntries).value);
        }
    }

    TCanvas canvas(("c" + prefix + "EfficiencyVsMargin").c_str(),
                   ("Pad " + quantity + " vs margin").c_str(),
                   900, 700);
    const std::array<int, kAxes> colors = {kBlue + 1, kRed + 1, kBlack};
    const std::array<int, kAxes> markers = {20, 21, 22};
    for (int axis = 0; axis < kAxes; ++axis) {
        graphs[axis].SetLineColor(colors[axis]);
        graphs[axis].SetMarkerColor(colors[axis]);
        graphs[axis].SetMarkerStyle(markers[axis]);
        graphs[axis].SetLineWidth(2);
    }
    graphs[kX].SetMinimum(0.0);
    graphs[kX].SetMaximum(isFake ? 0.05 : 1.02);
    graphs[kX].Draw("APL");
    graphs[kY].Draw("PL SAME");
    graphs[k2D].Draw("PL SAME");
    TLegend legend(0.65, 0.18, 0.88, 0.36);
    legend.SetBorderSize(0);
    legend.AddEntry(&graphs[kX], "X", "lp");
    legend.AddEntry(&graphs[kY], "Y", "lp");
    legend.AddEntry(&graphs[k2D], "2D", "lp");
    legend.Draw();
    canvas.Write();
}

double HuberLoss(double residual, double delta) {
    const double absolute = std::abs(residual);
    if (absolute <= delta) return 0.5 * residual * residual;
    return delta * (absolute - 0.5 * delta);
}

bool FitsParameter(const PadDUTAlignmentConfig& config,
                   const std::string& name) {
    return std::find(config.parameters.begin(), config.parameters.end(), name) !=
           config.parameters.end();
}

AlignmentMetrics MeasureAlignment(
    const std::vector<Event>& events,
    const std::shared_ptr<Detector>& detector,
    const planarPadConfig& padConfig,
    double margin) {
    AlignmentMetrics metrics;
    for (const auto& event : events) {
        const auto frameIt = event.detectorFramesMap.find(detector->GetID());
        if (frameIt == event.detectorFramesMap.end() ||
            frameIt->second->LocalHits().empty()) {
            continue;
        }
        ++metrics.candidates;
        const TVector3 predicted = detector->GlobalToLocal(
            detector->CalcHitFromTrack(event.track));
        const PadMatch match = FindMatchedPadHit(
            *frameIt->second, predicted, padConfig, margin);
        if (!match.IsValid()) continue;
        ++metrics.matches;
        metrics.residualsX.push_back(match.residualX);
        metrics.residualsY.push_back(match.residualY);
        metrics.meanX += match.residualX;
        metrics.meanY += match.residualY;
    }
    if (metrics.matches == 0) return metrics;
    metrics.meanX /= metrics.matches;
    metrics.meanY /= metrics.matches;
    for (size_t index = 0; index < metrics.residualsX.size(); ++index) {
        const double deltaX = metrics.residualsX[index] - metrics.meanX;
        const double deltaY = metrics.residualsY[index] - metrics.meanY;
        metrics.rmsX += deltaX * deltaX;
        metrics.rmsY += deltaY * deltaY;
    }
    metrics.rmsX = std::sqrt(metrics.rmsX / metrics.matches);
    metrics.rmsY = std::sqrt(metrics.rmsY / metrics.matches);
    return metrics;
}

std::vector<AlignmentSample> SelectAlignmentEvents(
    const std::vector<Event>& events,
    const std::shared_ptr<Detector>& detector,
    const planarPadConfig& padConfig,
    double margin) {
    std::vector<AlignmentSample> selected;
    selected.reserve(events.size());
    for (size_t index = 0; index < events.size(); ++index) {
        const auto frameIt = events[index].detectorFramesMap.find(
            detector->GetID());
        if (frameIt == events[index].detectorFramesMap.end() ||
            frameIt->second->LocalHits().empty()) {
            continue;
        }
        const TVector3 predicted = detector->GlobalToLocal(
            detector->CalcHitFromTrack(events[index].track));
        const PadMatch match = FindMatchedPadHit(
            *frameIt->second, predicted, padConfig, margin);
        if (match.IsValid()) selected.push_back({index, match.position});
    }
    return selected;
}

double AlignmentObjective(
    const double* parameters,
    const std::vector<Event>& events,
    const std::vector<AlignmentSample>& selectedEvents,
    const std::shared_ptr<Detector>& detector,
    const TVector3& baseAlignmentPosition,
    const TVector3& baseAlignmentRotation,
    double huberDelta) {
    detector->SetAlignment(
        baseAlignmentPosition.X() + parameters[0],
        baseAlignmentPosition.Y() + parameters[1],
        baseAlignmentPosition.Z() + parameters[2],
        baseAlignmentRotation.X() + parameters[3],
        baseAlignmentRotation.Y() + parameters[4],
        baseAlignmentRotation.Z() + parameters[5]);
    double objective = 0.0;
    long long count = 0;
    for (const auto& sample : selectedEvents) {
        const auto& event = events[sample.eventIndex];
        const TVector3 predicted = detector->GlobalToLocal(
            detector->CalcHitFromTrack(event.track));
        const double residualX = sample.hitPosition.X() - predicted.X();
        const double residualY = sample.hitPosition.Y() - predicted.Y();
        objective += HuberLoss(residualX, huberDelta) +
                     HuberLoss(residualY, huberDelta);
        ++count;
    }
    return count > 0 ? objective / count : 1e9;
}

AlignmentResult RunAlignmentPass(
    const std::vector<Event>& events,
    const std::shared_ptr<Detector>& detector,
    const PadDUTAlignmentConfig& requestedConfig,
    double margin) {
    AlignmentResult result;
    result.attempted = true;
    const auto* pad = detector->GetPlanarPadConfig();
    if (!pad) {
        result.status = "not a planar_pad detector";
        return result;
    }

    PadDUTAlignmentConfig config = requestedConfig;
    result.margin = margin;
    result.huberDelta = config.huberDelta > 0.0
                            ? config.huberDelta
                            : 0.5 * std::min(pad->pitchX, pad->pitchY);

    const TVector3 baseAlignmentPosition = detector->GetAlignPos();
    const TVector3 baseAlignmentRotation = detector->GetAlignRot();
    result.before = MeasureAlignment(
        events, detector, *pad, margin);
    const auto selectedEvents = SelectAlignmentEvents(
        events, detector, *pad, margin);
    if (selectedEvents.size() < static_cast<size_t>(config.minMatches)) {
        result.status = "insufficient initial matches";
        result.after = result.before;
        result.finalPosition = detector->GetPos();
        result.finalRotation = detector->GetRot();
        std::cerr << "[PadDUTAnalysis] DUT " << detector->GetID()
                  << " alignment skipped: " << selectedEvents.size()
                  << " initial matches, need at least " << config.minMatches
                  << '\n';
        return result;
    }

    const std::array<const char*, 6> names = {
        "dx", "dy", "dz", "rotX", "rotY", "rotZ"};
    const bool hasFreeParameter = std::any_of(
        names.begin(), names.end(), [&](const char* name) {
            return FitsParameter(config, name);
        });
    if (!hasFreeParameter) {
        result.status = "no free alignment parameters";
        result.after = result.before;
        result.finalPosition = detector->GetPos();
        result.finalRotation = detector->GetRot();
        return result;
    }

    auto minimizer = std::unique_ptr<ROOT::Math::Minimizer>(
        ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
    if (!minimizer) {
        result.status = "cannot create minimizer";
        result.after = result.before;
        return result;
    }
    auto objective = [&](const double* parameters) {
        return AlignmentObjective(
            parameters, events, selectedEvents, detector,
            baseAlignmentPosition, baseAlignmentRotation,
            result.huberDelta);
    };
    ROOT::Math::Functor function(objective, 6);
    minimizer->SetFunction(function);
    minimizer->SetTolerance(config.tolerance);
    minimizer->SetMaxFunctionCalls(config.maxFunctionCalls);
    minimizer->SetPrintLevel(0);
    for (int parameter = 0; parameter < 6; ++parameter) {
        if (!FitsParameter(config, names[parameter])) {
            minimizer->SetFixedVariable(parameter, names[parameter], 0.0);
            continue;
        }
        const bool translation = parameter < 3;
        const double limit = translation ? config.maxTranslation
                                         : config.maxRotation;
        const double step = translation ? 0.01 : 0.0001;
        minimizer->SetLimitedVariable(
            parameter, names[parameter], 0.0, step, -limit, limit);
    }

    const std::array<double, 6> zero{};
    result.objectiveBefore = objective(zero.data());
    const bool converged = minimizer->Minimize();
    const double* fitted = minimizer->X();
    if (!fitted) {
        detector->SetAlignment(
            baseAlignmentPosition.X(), baseAlignmentPosition.Y(),
            baseAlignmentPosition.Z(), baseAlignmentRotation.X(),
            baseAlignmentRotation.Y(), baseAlignmentRotation.Z());
        result.status = "minimizer returned no parameters";
        result.after = result.before;
        return result;
    }

    std::array<double, 6> correction{};
    std::copy(fitted, fitted + correction.size(), correction.begin());
    result.objectiveAfter = objective(correction.data());
    result.after = MeasureAlignment(
        events, detector, *pad, margin);
    const bool improved = std::isfinite(result.objectiveAfter) &&
                          result.objectiveAfter <= result.objectiveBefore;
    const bool enoughMatches =
        result.after.matches >= static_cast<long long>(config.minMatches);
    if (!converged || !improved || !enoughMatches) {
        detector->SetAlignment(
            baseAlignmentPosition.X(), baseAlignmentPosition.Y(),
            baseAlignmentPosition.Z(), baseAlignmentRotation.X(),
            baseAlignmentRotation.Y(), baseAlignmentRotation.Z());
        result.applied = false;
        result.status = !converged
                            ? "minimizer did not converge"
                            : !improved ? "objective did not improve"
                                        : "insufficient matches after fit";
        result.after = result.before;
        result.objectiveAfter = result.objectiveBefore;
    } else {
        result.applied = true;
        result.status = "applied";
        result.correctionPosition.SetXYZ(
            correction[0], correction[1], correction[2]);
        result.correctionRotation.SetXYZ(
            correction[3], correction[4], correction[5]);
    }
    result.finalPosition = detector->GetPos();
    result.finalRotation = detector->GetRot();

    return result;
}

AlignmentResult RunAlignment(
    const std::vector<Event>& events,
    const std::shared_ptr<Detector>& detector,
    const PadDUTAlignmentConfig& config,
    double margin) {
    AlignmentResult combined;
    combined.attempted = true;
    const TVector3 initialAlignmentPosition = detector->GetAlignPos();
    const TVector3 initialAlignmentRotation = detector->GetAlignRot();
    const auto* pad = detector->GetPlanarPadConfig();
    if (!pad) {
        combined.status = "not a planar_pad detector";
        return combined;
    }
    combined.margin = margin;
    combined.huberDelta = config.huberDelta > 0.0
                              ? config.huberDelta
                              : 0.5 * std::min(pad->pitchX, pad->pitchY);
    combined.before = MeasureAlignment(
        events, detector, *pad, margin);
    const auto diagnosticEvents = SelectAlignmentEvents(
        events, detector, *pad, margin);
    const std::array<double, 6> zero{};
    combined.objectiveBefore = diagnosticEvents.empty()
                                   ? 0.0
                                   : AlignmentObjective(
                                         zero.data(), events, diagnosticEvents,
                                         detector, initialAlignmentPosition,
                                         initialAlignmentRotation,
                                         combined.huberDelta);
    bool appliedAnyPass = false;

    for (int pass = 0; pass < config.maxPasses; ++pass) {
        AlignmentResult current = RunAlignmentPass(
            events, detector, config, margin);
        if (!current.applied) {
            combined.status = appliedAnyPass
                                  ? "applied; stopped: " + current.status
                                  : current.status;
            break;
        }

        appliedAnyPass = true;
        ++combined.passes;
        combined.after = current.after;
        combined.status = "applied";
        std::cout << std::fixed << std::setprecision(5)
                  << "[PadDUTAnalysis] DUT " << detector->GetID()
                  << " alignment pass " << combined.passes
                  << ": matches " << current.before.matches << " -> "
                  << current.after.matches
                  << ", mean residual [mm] (" << current.before.meanX
                  << ", " << current.before.meanY << ") -> ("
                  << current.after.meanX << ", " << current.after.meanY
                  << "), correction=[" << current.correctionPosition.X()
                  << ", " << current.correctionPosition.Y() << ", "
                  << current.correctionRotation.Z() << "]\n";

        const bool translationConverged =
            current.correctionPosition.Mag() <=
            config.convergenceTranslation;
        const bool rotationConverged =
            current.correctionRotation.Mag() <=
            config.convergenceRotation;
        if (translationConverged && rotationConverged) {
            combined.status = "applied; converged";
            break;
        }
    }

    combined.applied = appliedAnyPass;
    if (!appliedAnyPass) combined.after = combined.before;
    const TVector3 finalAlignmentPosition = detector->GetAlignPos();
    const TVector3 finalAlignmentRotation = detector->GetAlignRot();
    combined.correctionPosition =
        finalAlignmentPosition - initialAlignmentPosition;
    combined.correctionRotation =
        finalAlignmentRotation - initialAlignmentRotation;
    const std::array<double, 6> totalCorrection = {
        combined.correctionPosition.X(),
        combined.correctionPosition.Y(),
        combined.correctionPosition.Z(),
        combined.correctionRotation.X(),
        combined.correctionRotation.Y(),
        combined.correctionRotation.Z()};
    combined.objectiveAfter = diagnosticEvents.empty()
                                  ? 0.0
                                  : AlignmentObjective(
                                        totalCorrection.data(), events,
                                        diagnosticEvents, detector,
                                        initialAlignmentPosition,
                                        initialAlignmentRotation,
                                        combined.huberDelta);
    combined.finalPosition = detector->GetPos();
    combined.finalRotation = detector->GetRot();

    std::cout << std::fixed << std::setprecision(5)
              << "[PadDUTAnalysis] DUT " << detector->GetID()
              << " alignment " << combined.status
              << " after " << combined.passes << " pass(es)"
              << ": matches " << combined.before.matches << " -> "
              << combined.after.matches
              << ", mean residual [mm] (" << combined.before.meanX << ", "
              << combined.before.meanY << ") -> (" << combined.after.meanX
              << ", " << combined.after.meanY << ")"
              << ", RMS [mm] (" << combined.before.rmsX << ", "
              << combined.before.rmsY << ") -> (" << combined.after.rmsX
              << ", " << combined.after.rmsY << ")\n"
              << "[PadDUTAnalysis] DUT " << detector->GetID()
              << " total correction: position=["
              << combined.correctionPosition.X() << ", "
              << combined.correctionPosition.Y() << ", "
              << combined.correctionPosition.Z() << "], rotation=["
              << combined.correctionRotation.X() << ", "
              << combined.correctionRotation.Y() << ", "
              << combined.correctionRotation.Z() << "]\n"
              << "[PadDUTAnalysis] DUT " << detector->GetID()
              << " final geometry: position=[" << combined.finalPosition.X()
              << ", " << combined.finalPosition.Y() << ", "
              << combined.finalPosition.Z() << "], rotation=["
              << combined.finalRotation.X() << ", "
              << combined.finalRotation.Y() << ", "
              << combined.finalRotation.Z() << "]\n";
    return combined;
}

void WriteAlignmentDiagnostics(TDirectory* detectorDirectory,
                               int detectorID,
                               const AlignmentResult& result) {
    if (!detectorDirectory || !result.attempted) return;
    auto* directory = GetOrCreateDirectory(detectorDirectory, "Alignment");
    TDirectory::TContext context(directory);

    const auto residualRange = [](const std::vector<double>& before,
                                  const std::vector<double>& after) {
        double maximum = 0.0;
        for (double value : before) maximum = std::max(maximum, std::abs(value));
        for (double value : after) maximum = std::max(maximum, std::abs(value));
        return std::max(1.0, 1.05 * maximum);
    };
    const double rangeX = residualRange(
        result.before.residualsX, result.after.residualsX);
    const double rangeY = residualRange(
        result.before.residualsY, result.after.residualsY);

    TH1D residualXBefore("hResidualXBefore",
                         "X residual before alignment;Residual X [mm];Events",
                         120, -rangeX, rangeX);
    TH1D residualXAfter("hResidualXAfter",
                        "X residual after alignment;Residual X [mm];Events",
                        120, -rangeX, rangeX);
    TH1D residualYBefore("hResidualYBefore",
                         "Y residual before alignment;Residual Y [mm];Events",
                         120, -rangeY, rangeY);
    TH1D residualYAfter("hResidualYAfter",
                        "Y residual after alignment;Residual Y [mm];Events",
                        120, -rangeY, rangeY);
    for (double value : result.before.residualsX) residualXBefore.Fill(value);
    for (double value : result.after.residualsX) residualXAfter.Fill(value);
    for (double value : result.before.residualsY) residualYBefore.Fill(value);
    for (double value : result.after.residualsY) residualYAfter.Fill(value);
    residualXBefore.Write();
    residualXAfter.Write();
    residualYBefore.Write();
    residualYAfter.Write();

    const auto writeComparison = [](TH1D& before, TH1D& after,
                                    const char* canvasName,
                                    const char* canvasTitle) {
        before.SetStats(false);
        after.SetStats(false);
        before.SetLineColor(kBlue + 1);
        after.SetLineColor(kRed + 1);
        before.SetLineWidth(2);
        after.SetLineWidth(2);
        TCanvas canvas(canvasName, canvasTitle, 900, 700);
        const double maximum = std::max(before.GetMaximum(), after.GetMaximum());
        before.SetMaximum(1.15 * maximum);
        before.Draw("HIST");
        after.Draw("HIST SAME");
        TLegend legend(0.66, 0.74, 0.88, 0.88);
        legend.SetBorderSize(0);
        legend.AddEntry(&before, "Before", "l");
        legend.AddEntry(&after, "After", "l");
        legend.Draw();
        canvas.Write();
    };
    writeComparison(residualXBefore, residualXAfter,
                    "cResidualXBeforeAfter", "X residual before/after alignment");
    writeComparison(residualYBefore, residualYAfter,
                    "cResidualYBeforeAfter", "Y residual before/after alignment");

    Int_t dutID = detectorID;
    Bool_t applied = result.applied;
    Int_t passes = result.passes;
    Long64_t matchesBefore = result.before.matches;
    Long64_t matchesAfter = result.after.matches;
    Double_t margin = result.margin;
    Double_t objectiveBefore = result.objectiveBefore;
    Double_t objectiveAfter = result.objectiveAfter;
    Double_t meanXBefore = result.before.meanX;
    Double_t meanXAfter = result.after.meanX;
    Double_t meanYBefore = result.before.meanY;
    Double_t meanYAfter = result.after.meanY;
    Double_t rmsXBefore = result.before.rmsX;
    Double_t rmsXAfter = result.after.rmsX;
    Double_t rmsYBefore = result.before.rmsY;
    Double_t rmsYAfter = result.after.rmsY;
    Double_t dx = result.correctionPosition.X();
    Double_t dy = result.correctionPosition.Y();
    Double_t dz = result.correctionPosition.Z();
    Double_t rotX = result.correctionRotation.X();
    Double_t rotY = result.correctionRotation.Y();
    Double_t rotZ = result.correctionRotation.Z();
    TTree summary("AlignmentSummary", "Pad DUT alignment summary");
    summary.Branch("dutID", &dutID);
    summary.Branch("applied", &applied);
    summary.Branch("passes", &passes);
    summary.Branch("matchesBefore", &matchesBefore);
    summary.Branch("matchesAfter", &matchesAfter);
    summary.Branch("margin", &margin);
    summary.Branch("objectiveBefore", &objectiveBefore);
    summary.Branch("objectiveAfter", &objectiveAfter);
    summary.Branch("meanXBefore", &meanXBefore);
    summary.Branch("meanXAfter", &meanXAfter);
    summary.Branch("meanYBefore", &meanYBefore);
    summary.Branch("meanYAfter", &meanYAfter);
    summary.Branch("rmsXBefore", &rmsXBefore);
    summary.Branch("rmsXAfter", &rmsXAfter);
    summary.Branch("rmsYBefore", &rmsYBefore);
    summary.Branch("rmsYAfter", &rmsYAfter);
    summary.Branch("dx", &dx);
    summary.Branch("dy", &dy);
    summary.Branch("dz", &dz);
    summary.Branch("rotX", &rotX);
    summary.Branch("rotY", &rotY);
    summary.Branch("rotZ", &rotZ);
    summary.Fill();
    summary.Write();
}

}  // namespace

void PadDUTAnalysisScript::LoadConfig(const json& config) {
    m_runAlignment = config.value("runAlignment", m_runAlignment);
    m_progressInterval = config.value("progressInterval", m_progressInterval);
    m_maxEvents = config.value("maxEvents", m_maxEvents);

    const json alignment = config.value("alignment", json::object());
    m_alignment.parameters = alignment.value(
        "parameters", m_alignment.parameters);
    const std::array<std::string, 6> validParameters = {
        "dx", "dy", "dz", "rotX", "rotY", "rotZ"};
    m_alignment.parameters.erase(
        std::remove_if(
            m_alignment.parameters.begin(), m_alignment.parameters.end(),
            [&](const std::string& name) {
                return std::find(validParameters.begin(), validParameters.end(),
                                 name) == validParameters.end();
            }),
        m_alignment.parameters.end());
    std::sort(m_alignment.parameters.begin(), m_alignment.parameters.end());
    m_alignment.parameters.erase(
        std::unique(m_alignment.parameters.begin(),
                    m_alignment.parameters.end()),
        m_alignment.parameters.end());
    m_alignment.minMatches = std::max(
        1, alignment.value("minMatches", m_alignment.minMatches));
    m_alignment.huberDelta = alignment.value(
        "huberDelta", m_alignment.huberDelta);
    m_alignment.maxTranslation = std::max(
        1e-6, alignment.value(
                  "maxTranslation", m_alignment.maxTranslation));
    m_alignment.maxRotation = std::max(
        1e-8, alignment.value("maxRotation", m_alignment.maxRotation));
    m_alignment.tolerance = std::max(
        1e-9, alignment.value("tolerance", m_alignment.tolerance));
    m_alignment.maxFunctionCalls = std::max(
        1, alignment.value(
               "maxFunctionCalls", m_alignment.maxFunctionCalls));
    m_alignment.maxPasses = std::max(
        1, alignment.value("maxPasses", m_alignment.maxPasses));
    m_alignment.convergenceTranslation = std::max(
        0.0, alignment.value(
                 "convergenceTranslation",
                 m_alignment.convergenceTranslation));
    m_alignment.convergenceRotation = std::max(
        0.0, alignment.value(
                 "convergenceRotation", m_alignment.convergenceRotation));

    const json efficiency = config.value("efficiencyMap", json::object());
    m_effXBins = std::max(1, efficiency.value("xBins", m_effXBins));
    m_effYBins = std::max(1, efficiency.value("yBins", m_effYBins));
    m_effMinEntriesPerBin = std::max(
        1, efficiency.value("minEntriesPerBin", m_effMinEntriesPerBin));
    m_effXMin = efficiency.value("xMin", m_effXMin);
    m_effXMax = efficiency.value("xMax", m_effXMax);
    m_effYMin = efficiency.value("yMin", m_effYMin);
    m_effYMax = efficiency.value("yMax", m_effYMax);
    if (m_effXMax < m_effXMin) std::swap(m_effXMin, m_effXMax);
    if (m_effYMax < m_effYMin) std::swap(m_effYMin, m_effYMax);
    if (m_effXMax == m_effXMin) m_effXMax = m_effXMin + 1.0;
    if (m_effYMax == m_effYMin) m_effYMax = m_effYMin + 1.0;
    m_effExcludedXBins = efficiency.value(
        "excludeXBins", std::vector<int>{});
    m_effExcludedYBins = efficiency.value(
        "excludeYBins", std::vector<int>{});
    const auto cleanBins = [](std::vector<int>& bins, int maximum) {
        bins.erase(std::remove_if(bins.begin(), bins.end(),
                                  [maximum](int bin) {
                                      return bin < 1 || bin > maximum;
                                  }),
                   bins.end());
        std::sort(bins.begin(), bins.end());
        bins.erase(std::unique(bins.begin(), bins.end()), bins.end());
    };
    cleanBins(m_effExcludedXBins, m_effXBins);
    cleanBins(m_effExcludedYBins, m_effYBins);

    m_margin = efficiency.value("margin", m_margin);
    const json scan = efficiency.value("marginScan", json::object());
    m_marginScanMin = scan.value("min", m_marginScanMin);
    m_marginScanMax = scan.value("max", m_marginScanMax);
    m_marginScanStep = scan.value("step", m_marginScanStep);
    const json fake = efficiency.value("fake", json::object());
    m_enableFakeEfficiency = fake.value("enabled", m_enableFakeEfficiency);
    m_fakeSeed = fake.value("seed", m_fakeSeed);
    m_fakePartnersPerEvent = std::max(
        1, fake.value("partnersPerEvent", m_fakePartnersPerEvent));
}

void PadDUTAnalysisScript::Print() const {
    if (!Terminal::Verbose()) return;
    std::cout << "PadDUTAnalysisScript Configuration:\n"
              << "  Run Alignment: " << (m_runAlignment ? "Yes" : "No") << '\n'
              << "  Alignment Parameters:";
    if (m_alignment.parameters.empty()) {
        std::cout << " none";
    } else {
        for (const auto& parameter : m_alignment.parameters) {
            std::cout << ' ' << parameter;
        }
    }
    std::cout << '\n'
              << "  Alignment Min Matches: " << m_alignment.minMatches << '\n'
              << "  Alignment Max Passes: " << m_alignment.maxPasses << '\n'
              << "  Cluster Match: track inside cluster envelope + margin\n"
              << "  Efficiency X: [" << m_effXMin << ", " << m_effXMax
              << "] mm, bins=" << m_effXBins << '\n'
              << "  Efficiency Y: [" << m_effYMin << ", " << m_effYMax
              << "] mm, bins=" << m_effYBins << '\n'
              << "  Pad Envelope Margin: " << m_margin << " mm\n"
              << "  Fake Efficiency: "
              << (m_enableFakeEfficiency ? "enabled" : "disabled") << '\n'
              << "  Max Events: "
              << (m_maxEvents > 0 ? std::to_string(m_maxEvents) : "All")
              << '\n';
}

bool PadDUTAnalysisScript::Execute() {
    auto parser = GetParser();
    if (!parser) {
        std::cerr << "[PadDUTAnalysis] Parser is not set\n";
        return false;
    }
    Print();

    auto& factory = DetectorFactory::GetInstance();
    std::vector<std::shared_ptr<Detector>> padDUTs;
    for (const auto& detector :
         factory.GetDetectorsByRole(Detector::Role::DUT)) {
        if (detector->GetPlanarPadConfig()) {
            padDUTs.push_back(detector);
        } else {
            std::cerr << "[PadDUTAnalysis] WARNING: skipping non-pad DUT "
                      << detector->GetID() << '\n';
        }
    }
    if (padDUTs.empty()) {
        std::cerr << "[PadDUTAnalysis] No planar_pad DUT is configured\n";
        return false;
    }

    const std::string trackPath = GetOutputDir() + "TrackInfo.root";
    TFile trackFile(trackPath.c_str(), "READ");
    auto* trackTree = dynamic_cast<TTree*>(trackFile.Get("Tracks"));
    if (trackFile.IsZombie() || !trackTree) {
        std::cerr << "[PadDUTAnalysis] Cannot read Tracks from "
                  << trackPath << '\n';
        return false;
    }

    Int_t eventID = 0;
    Track* track = nullptr;
    double eventTime = 0.0;
    trackTree->SetBranchStatus("*", false);
    trackTree->SetBranchStatus("eventID", true);
    trackTree->SetBranchAddress("eventID", &eventID);
    std::unordered_map<int, int> trackMultiplicity;
    for (Long64_t entry = 0; entry < trackTree->GetEntries(); ++entry) {
        trackTree->GetEntry(entry);
        ++trackMultiplicity[eventID];
    }
    trackTree->SetBranchStatus("*", true);
    trackTree->SetBranchAddress("track", &track);
    if (trackTree->GetBranch("t0")) {
        trackTree->SetBranchAddress("t0", &eventTime);
    }

    const auto& trackers = factory.GetDetectorsByRole(Detector::Role::Tracker);
    std::unordered_set<int> strictSingleHitTrackerEvents;
    std::vector<Event> events;
    const Long64_t maximumEvents = m_maxEvents > 0
                                       ? m_maxEvents
                                       : std::numeric_limits<Long64_t>::max();
    Long64_t processed = 0;
    for (Long64_t entry = 0;
         entry < trackTree->GetEntries() && processed < maximumEvents;
         ++entry) {
        trackTree->GetEntry(entry);
        if (!track || trackMultiplicity[eventID] != 1) continue;
        auto rawByDetector = parser->LoadEvent(eventID);
        if (rawByDetector.empty()) continue;

        if (m_enableFakeEfficiency) {
            bool strict = true;
            for (const auto& trackerDetector : trackers) {
                DetectorFrame frame(*trackerDetector);
                const auto raw = rawByDetector.find(trackerDetector->GetID());
                if (raw != rawByDetector.end()) frame.SetRawData(raw->second);
                frame.Process();
                if (frame.LocalHits().size() != 1) {
                    strict = false;
                    break;
                }
            }
            if (strict) strictSingleHitTrackerEvents.insert(eventID);
        }

        Event event{.eventID = eventID, .track = *track, .t0 = eventTime};
        for (const auto& detector : padDUTs) {
            auto frame = std::make_shared<DetectorFrame>(*detector);
            const auto raw = rawByDetector.find(detector->GetID());
            if (raw != rawByDetector.end()) frame->SetRawData(raw->second);
            frame->Process(eventTime);
            event.detectorFramesMap[detector->GetID()] = std::move(frame);
        }
        events.push_back(std::move(event));
        ++processed;
        if (Terminal::Interactive() && m_progressInterval > 0 &&
            processed % (5 * m_progressInterval) == 0) {
            std::cout << "\r      processing "
                      << Terminal::Count(processed)
                      << " single-track events" << std::flush;
        }
    }
    Terminal::ClearProgress();

    std::unordered_map<int, AlignmentResult> alignmentResults;
    if (m_runAlignment) {
        for (const auto& detector : padDUTs) {
            alignmentResults.emplace(
                detector->GetID(),
                RunAlignment(events, detector, m_alignment, m_margin));
        }
    }

    const std::string outputPath = GetOutputDir() + "PadDUTInfo.root";
    TFile outputFile(outputPath.c_str(), "RECREATE");
    if (outputFile.IsZombie()) {
        std::cerr << "[PadDUTAnalysis] Cannot create " << outputPath << '\n';
        return false;
    }

    Int_t dutID = 0;
    Double_t predX = kInvalidValue, predY = kInvalidValue;
    Double_t hitX = kInvalidValue, hitY = kInvalidValue;
    Double_t resX = kInvalidValue, resY = kInvalidValue;
    Int_t hitFlag = 0;
    Int_t clusterIndex = -1;
    std::vector<ChannelHit> channelHits;
    std::vector<Cluster> clusters;
    Cluster selectedCluster{};
    std::vector<ChannelHit> selectedChannelHits;

    TTree outputTree("PadDUTTree", "Two-dimensional pad DUT data");
    outputTree.Branch("eventID", &eventID);
    outputTree.Branch("dutID", &dutID);
    outputTree.Branch("predX", &predX);
    outputTree.Branch("predY", &predY);
    outputTree.Branch("hitX", &hitX);
    outputTree.Branch("hitY", &hitY);
    outputTree.Branch("resX", &resX);
    outputTree.Branch("resY", &resY);
    outputTree.Branch("hitFlag", &hitFlag);
    outputTree.Branch("clusterIndex", &clusterIndex);
    outputTree.Branch("channelHits", &channelHits);
    outputTree.Branch("clusters", &clusters);
    outputTree.Branch("selectedCluster", &selectedCluster);
    outputTree.Branch("selectedChannelHits", &selectedChannelHits);

    for (const auto& detector : padDUTs) {
        dutID = detector->GetID();
        const auto& padConfig = *detector->GetPlanarPadConfig();
        std::vector<double> residualsX;
        std::vector<double> residualsY;
        std::vector<PadAlignmentQAPoint> alignmentQAPoints;
        std::vector<EfficiencySample> samples;
        std::vector<double> matchedClusterCharges;
        std::vector<int> matchedClusterSizes;
        std::vector<double> matchedHitADCs;
        TH2D reconstructedHitMap(
            "hReconstructedHitMap",
            "All reconstructed DUT hits;Reconstructed X [mm];Reconstructed Y [mm]",
            m_effXBins, m_effXMin, m_effXMax,
            m_effYBins, m_effYMin, m_effYMax);
        TH2D matchedHitMap(
            "hMatchedHitMap",
            "Track-matched reconstructed DUT hits;Reconstructed X [mm];Reconstructed Y [mm]",
            m_effXBins, m_effXMin, m_effXMax,
            m_effYBins, m_effYMin, m_effYMax);
        for (TH1* histogram : {
                 static_cast<TH1*>(&reconstructedHitMap),
                 static_cast<TH1*>(&matchedHitMap)}) {
            histogram->SetDirectory(nullptr);
            histogram->SetStats(false);
            histogram->SetOption("COLZ TEXT");
        }
        for (const auto& event : events) {
            eventID = event.eventID;
            predX = predY = hitX = hitY = resX = resY = kInvalidValue;
            hitFlag = 0;
            clusterIndex = -1;
            channelHits.clear();
            clusters.clear();
            selectedCluster = Cluster{};
            selectedChannelHits.clear();

            const TVector3 predicted = detector->GlobalToLocal(
                detector->CalcHitFromTrack(event.track));
            predX = predicted.X();
            predY = predicted.Y();
            const auto frameIt = event.detectorFramesMap.find(dutID);
            if (frameIt != event.detectorFramesMap.end()) {
                const auto& frame = *frameIt->second;
                channelHits = frame.ChannelHits();
                clusters = frame.Clusters();
                for (const auto& localHit : frame.LocalHits()) {
                    if (std::isfinite(localHit.localPos.X()) &&
                        std::isfinite(localHit.localPos.Y())) {
                        reconstructedHitMap.Fill(
                            localHit.localPos.X(), localHit.localPos.Y());
                    }
                }
                const PadMatch match = FindMatchedPadHit(
                    frame, predicted, padConfig,
                    m_margin);
                if (match.IsValid()) {
                    hitFlag = 3;
                    clusterIndex = match.clusterIndex;
                    hitX = match.position.X();
                    hitY = match.position.Y();
                    resX = match.residualX;
                    resY = match.residualY;
                    selectedCluster = clusters[static_cast<size_t>(clusterIndex)];
                    if (std::isfinite(selectedCluster.charge))
                        matchedClusterCharges.push_back(selectedCluster.charge);
                    matchedClusterSizes.push_back(selectedCluster.size);
                    matchedHitMap.Fill(hitX, hitY);
                    for (int index : selectedCluster.channelHitIndices) {
                        if (index >= 0 &&
                            static_cast<size_t>(index) < channelHits.size()) {
                            const auto& selectedHit =
                                channelHits[static_cast<size_t>(index)];
                            selectedChannelHits.push_back(selectedHit);
                            if (selectedHit.isValid &&
                                std::isfinite(selectedHit.amp))
                                matchedHitADCs.push_back(selectedHit.amp);
                        }
                    }
                    residualsX.push_back(resX);
                    residualsY.push_back(resY);
                    alignmentQAPoints.push_back(
                        {predX, predY, resX, resY});
                }

                const int flatBin = FindFlatBin(
                    predX, predY,
                    m_effXMin, m_effXMax, m_effXBins,
                    m_effYMin, m_effYMax, m_effYBins,
                    m_effExcludedXBins, m_effExcludedYBins);
                if (flatBin >= 0) {
                    EfficiencySample sample;
                    sample.eventID = eventID;
                    sample.flatBin = flatBin;
                    sample.predX = predX;
                    sample.predY = predY;
                    sample.envelopes = BuildPadEnvelopes(
                        frame, padConfig);
                    samples.push_back(std::move(sample));
                }
            }
            outputTree.Fill();
        }

        auto* detectorDirectory = GetOrCreateDirectory(
            &outputFile, "DUT_" + std::to_string(dutID));
        WriteDUTHitProperties(
            detectorDirectory, dutID, matchedClusterCharges,
            matchedClusterSizes, matchedHitADCs);
        if (Terminal::Verbose()) {
            Terminal::Detail(
                "DUT" + std::to_string(dutID) + " hit properties · " +
                Terminal::Count(matchedClusterSizes.size()) +
                " clusters · " + Terminal::Count(matchedHitADCs.size()) +
                " channel hits");
        }
        const auto alignment = alignmentResults.find(dutID);
        if (alignment != alignmentResults.end()) {
            WriteAlignmentDiagnostics(
                detectorDirectory, dutID, alignment->second);
        }
        auto* residualDirectory = GetOrCreateDirectory(
            detectorDirectory, "Residuals");
        const ResidualResult resultX = WriteResidual(
            residualsX, residualDirectory, "hResidualX");
        const ResidualResult resultY = WriteResidual(
            residualsY, residualDirectory, "hResidualY");
        WriteResidualMaps(
            residualDirectory, alignmentQAPoints,
            m_effXMin, m_effXMax, m_effYMin, m_effYMax);

        auto* hitMapDirectory = GetOrCreateDirectory(
            detectorDirectory, "HitMaps");
        {
            TDirectory::TContext context(hitMapDirectory);
            reconstructedHitMap.Write();
            matchedHitMap.Write();
        }

        auto* efficiencyDirectory = GetOrCreateDirectory(
            detectorDirectory, "Efficiency");
        const EfficiencyGrid efficiency = BuildEfficiencyGrid(
            samples, m_effXBins, m_effYBins, m_margin);
        WriteEfficiencyMaps(
            efficiencyDirectory, efficiency,
            m_effXBins, m_effXMin, m_effXMax,
            m_effYBins, m_effYMin, m_effYMax,
            m_effMinEntriesPerBin, "h");
        WriteEfficiencyProjections(
            efficiencyDirectory, efficiency,
            m_effXBins, m_effXMin, m_effXMax,
            m_effYBins, m_effYMin, m_effYMax,
            m_effMinEntriesPerBin, "");
        WriteMarginScan(
            efficiencyDirectory,
            [&](double margin) {
                return BuildEfficiencyGrid(
                    samples, m_effXBins, m_effYBins, margin);
            },
            m_effMinEntriesPerBin,
            m_marginScanMin, m_marginScanMax, m_marginScanStep, "");

        EfficiencyEstimate fake2D;
        if (m_enableFakeEfficiency) {
            std::vector<const EfficiencySample*> strictSamples;
            for (const auto& sample : samples) {
                if (strictSingleHitTrackerEvents.count(sample.eventID)) {
                    strictSamples.push_back(&sample);
                }
            }
            std::vector<FakeEfficiencySample> fakeSamples;
            if (strictSamples.size() > 1) {
                std::mt19937 random(m_fakeSeed);
                std::uniform_int_distribution<size_t> distribution(
                    0, strictSamples.size() - 2);
                const size_t partners = std::min(
                    static_cast<size_t>(m_fakePartnersPerEvent),
                    strictSamples.size() - 1);
                for (size_t sourceIndex = 0;
                     sourceIndex < strictSamples.size(); ++sourceIndex) {
                    const auto& source = *strictSamples[sourceIndex];
                    std::set<size_t> selectedPartners;
                    while (selectedPartners.size() < partners) {
                        size_t partner = distribution(random);
                        if (partner >= sourceIndex) ++partner;
                        selectedPartners.insert(partner);
                    }
                    for (size_t partnerIndex : selectedPartners) {
                        const auto& partner = *strictSamples[partnerIndex];
                        fakeSamples.push_back(
                            {partner.flatBin,
                             &source,
                             partner.predX, partner.predY});
                    }
                }
            }
            const EfficiencyGrid fakeGrid = BuildFakeEfficiencyGrid(
                fakeSamples, m_effXBins, m_effYBins, m_margin);
            auto* fakeDirectory = GetOrCreateDirectory(
                detectorDirectory, "Fake");
            WriteEfficiencyMaps(
                fakeDirectory, fakeGrid,
                m_effXBins, m_effXMin, m_effXMax,
                m_effYBins, m_effYMin, m_effYMax,
                m_effMinEntriesPerBin, "hFake");
            WriteEfficiencyProjections(
                fakeDirectory, fakeGrid,
                m_effXBins, m_effXMin, m_effXMax,
                m_effYBins, m_effYMin, m_effYMax,
                m_effMinEntriesPerBin, "Fake");
            fake2D = EstimateEfficiency(
                fakeGrid, k2D, m_effMinEntriesPerBin);
        }

        const auto efficiencyX = EstimateEfficiency(
            efficiency, kX, m_effMinEntriesPerBin);
        const auto efficiencyY = EstimateEfficiency(
            efficiency, kY, m_effMinEntriesPerBin);
        const auto efficiency2D = EstimateEfficiency(
            efficiency, k2D, m_effMinEntriesPerBin);
        std::ostringstream summary;
        summary << "DUT" << dutID << " · "
                << Terminal::Count(matchedClusterSizes.size()) << '/'
                << Terminal::Count(efficiency2D.total) << " matched · "
                << std::fixed << std::setprecision(2)
                << "efficiency " << 100.0 * efficiency2D.value << "% · "
                << std::setprecision(3)
                << "fake " << 100.0 * fake2D.value << "% · "
                << std::setprecision(2)
                << "σ68 " << resultX.sigma68 << " × "
                << resultY.sigma68 << " mm";
        Terminal::Detail(summary.str());
        if (Terminal::Verbose()) {
            std::ostringstream axes;
            axes << std::fixed << std::setprecision(2)
                 << "DUT" << dutID << " efficiency X/Y "
                 << 100.0 * efficiencyX.value << "% / "
                 << 100.0 * efficiencyY.value << '%';
            Terminal::Detail(Terminal::Muted(axes.str()));
        }
    }

    outputFile.cd();
    outputTree.Write();
    outputFile.Close();
    trackFile.Close();
    return true;
}
