#include "Script/DUTEfficiencyAnalysis.h"
#include "Event/DetectorFrame.h"

#include <TCanvas.h>
#include <TEfficiency.h>
#include <TGraphAsymmErrors.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TMath.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iomanip>
#include <limits>
#include <random>
#include <sstream>
#include <string>

namespace DUTEfficiency {
namespace {

constexpr int kAxes = 3;
constexpr int kX = 0;
constexpr int kY = 1;
constexpr int k2D = 2;

struct Interval {
    double low = 0.0;
    double high = 0.0;
};

struct Sample {
    int eventID = 0;
    double predX = 0.0;
    double predY = 0.0;
    int flatBin = -1;
    std::vector<Interval> intervalsX;
    std::vector<Interval> intervalsY;
};

struct Estimate {
    double value = 0.0;
    double errorLow = 0.0;
    double errorHigh = 0.0;
    long long total = 0;
    long long passed = 0;
    int bins = 0;
};

struct Nonuniformity {
    double value = 0.0;
    double mean = 0.0;
    double standardDeviation = 0.0;
    int bins = 0;
};

struct Grid {
    std::vector<long long> denominator;
    std::array<std::vector<long long>, kAxes> numerator;

    Grid(int nx, int ny)
        : denominator(nx * ny, 0) {
        for (auto& values : numerator) values.assign(nx * ny, 0);
    }

    void Fill(int bin, const std::array<bool, kAxes>& matched) {
        if (bin < 0 || static_cast<size_t>(bin) >= denominator.size()) return;
        ++denominator[bin];
        for (int axis = 0; axis < kAxes; ++axis)
            if (matched[axis]) ++numerator[axis][bin];
    }
};

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

bool Contains(const std::vector<int>& bins, int value) {
    return std::find(bins.begin(), bins.end(), value) != bins.end();
}

int FindFlatBin(double x, double y, const Config& config) {
    if (x < config.xMin || x >= config.xMax ||
        y < config.yMin || y >= config.yMax) return -1;
    const int bx = static_cast<int>((x - config.xMin) /
                                    (config.xMax - config.xMin) * config.xBins) + 1;
    const int by = static_cast<int>((y - config.yMin) /
                                    (config.yMax - config.yMin) * config.yBins) + 1;
    if (bx < 1 || bx > config.xBins || by < 1 || by > config.yBins ||
        Contains(config.excludedXBins, bx) || Contains(config.excludedYBins, by)) return -1;
    return (bx - 1) * config.yBins + (by - 1);
}

Estimate Calculate(const Grid& grid, int axis, int minimumEntries) {
    Estimate result;
    for (size_t bin = 0; bin < grid.denominator.size(); ++bin) {
        const long long total = grid.denominator[bin];
        if (total < minimumEntries) continue;
        const long long passed = grid.numerator[axis][bin];
        result.total += total;
        result.passed += passed;
        ++result.bins;
    }
    if (result.total == 0) return result;
    result.value = static_cast<double>(result.passed) / result.total;
    const double lower = TEfficiency::ClopperPearson(
        result.total, result.passed, 0.683, false);
    const double upper = TEfficiency::ClopperPearson(
        result.total, result.passed, 0.683, true);
    result.errorLow = result.value - lower;
    result.errorHigh = upper - result.value;
    return result;
}

Nonuniformity CalculateNonuniformity(const Grid& grid, int axis,
                                     int minimumEntries) {
    Nonuniformity result;
    for (size_t bin = 0; bin < grid.denominator.size(); ++bin) {
        const long long total = grid.denominator[bin];
        if (total < minimumEntries) continue;
        result.mean += static_cast<double>(grid.numerator[axis][bin]) / total;
        ++result.bins;
    }
    if (result.bins == 0) return result;
    result.mean /= result.bins;

    double squaredDeviation = 0.0;
    for (size_t bin = 0; bin < grid.denominator.size(); ++bin) {
        const long long total = grid.denominator[bin];
        if (total < minimumEntries) continue;
        const double efficiency =
            static_cast<double>(grid.numerator[axis][bin]) / total;
        const double deviation = efficiency - result.mean;
        squaredDeviation += deviation * deviation;
    }
    result.standardDeviation = std::sqrt(squaredDeviation / result.bins);
    if (result.mean > 0.0) {
        result.value = result.standardDeviation / result.mean;
    }
    return result;
}

double RequiredExtension(double coordinate, const std::vector<Interval>& intervals) {
    double required = std::numeric_limits<double>::infinity();
    for (const auto& interval : intervals) {
        const double distance = coordinate < interval.low
                                    ? interval.low - coordinate
                                    : coordinate > interval.high
                                          ? coordinate - interval.high
                                          : 0.0;
        required = std::min(required, distance);
    }
    return required;
}

std::vector<Interval> RemoveMatchedClusters(
    const std::vector<Interval>& intervals, double coordinate, double margin) {
    std::vector<Interval> remaining;
    remaining.reserve(intervals.size());
    for (const auto& interval : intervals) {
        const double extension = coordinate < interval.low
                                     ? interval.low - coordinate
                                     : coordinate > interval.high
                                           ? coordinate - interval.high
                                           : 0.0;
        if (extension > margin) remaining.push_back(interval);
    }
    return remaining;
}

std::array<bool, kAxes> EnvelopeMatch(const Sample& sample, double margin) {
    const bool x = RequiredExtension(sample.predX, sample.intervalsX) <=
                   margin;
    const bool y = RequiredExtension(sample.predY, sample.intervalsY) <=
                   margin;
    return {x, y, x && y};
}

std::vector<Interval> BuildIntervals(const DetectorFrame& frame, int type,
                                     double pitch) {
    std::vector<Interval> result;
    const auto& hits = frame.ChannelHits();
    for (const auto& cluster : frame.Clusters()) {
        if (cluster.type != type || cluster.channelHitIndices.empty()) continue;
        int minimum = std::numeric_limits<int>::max();
        int maximum = std::numeric_limits<int>::min();
        for (int index : cluster.channelHitIndices) {
            if (index < 0 || static_cast<size_t>(index) >= hits.size()) continue;
            minimum = std::min(minimum, hits[index].id0);
            maximum = std::max(maximum, hits[index].id0);
        }
        if (minimum <= maximum) {
            result.push_back({minimum * pitch, maximum * pitch});
        }
    }
    return result;
}

TDirectory* GetOrCreate(TDirectory* parent, const std::string& name) {
    if (auto* existing = parent->GetDirectory(name.c_str())) return existing;
    return parent->mkdir(name.c_str());
}

void StyleGraph(TGraphAsymmErrors& graph, int color, int marker);

void WriteEventCount(TDirectory* directory, const Grid& grid,
                     const Config& config) {
    TDirectory::TContext context(directory);
    TH2D eventCount("hEventCount", "DUT event count;Local X [mm];Local Y [mm];Events",
                    config.xBins, config.xMin, config.xMax,
                    config.yBins, config.yMin, config.yMax);
    for (int bx = 1; bx <= config.xBins; ++bx) {
        for (int by = 1; by <= config.yBins; ++by) {
            const int flat = (bx - 1) * config.yBins + (by - 1);
            eventCount.SetBinContent(bx, by, grid.denominator[flat]);
        }
    }
    eventCount.SetStats(false);
    eventCount.SetMarkerSize(1.1);
    TCanvas canvas("cEventCount", "DUT event count", 900, 700);
    canvas.SetRightMargin(0.15);
    eventCount.Draw("COLZ TEXT");
    canvas.Write();
}

std::string PercentLabel(const std::string& label, double value) {
    std::ostringstream stream;
    stream << label << ": " << std::fixed << std::setprecision(2)
           << 100.0 * value << "%";
    return stream.str();
}

void WriteMaps(TDirectory* directory, const Grid& grid,
               const Config& config) {
    TDirectory::TContext context(directory);
    std::array<std::unique_ptr<TH2D>, kAxes> efficiency;
    const std::array<const char*, kAxes> labels = {"X", "Y", "2D"};
    for (int axis = 0; axis < kAxes; ++axis) {
        efficiency[axis] = std::make_unique<TH2D>(
            ("hEff" + std::string(labels[axis])).c_str(),
            (std::string(labels[axis]) +
             " efficiency map;Local X [mm];Local Y [mm];Efficiency").c_str(),
            config.xBins, config.xMin, config.xMax,
            config.yBins, config.yMin, config.yMax);
    }
    for (int bx = 1; bx <= config.xBins; ++bx) {
        for (int by = 1; by <= config.yBins; ++by) {
            const int flat = (bx - 1) * config.yBins + (by - 1);
            const long long total = grid.denominator[flat];
            const bool valid = total >= config.minEntriesPerBin;
            for (int axis = 0; axis < kAxes; ++axis) {
                efficiency[axis]->SetBinContent(
                    bx, by, valid ? static_cast<double>(grid.numerator[axis][flat]) / total : -1.0);
            }
        }
    }
    for (int axis = 0; axis < kAxes; ++axis) {
        efficiency[axis]->SetMinimum(0.0);
        efficiency[axis]->SetMaximum(1.0);
        efficiency[axis]->SetStats(false);
        efficiency[axis]->SetMarkerSize(1.1);
    }

    for (int axis = 0; axis < kAxes; ++axis) {
        const std::string canvasName = "cEffMap" + std::string(labels[axis]);
        TCanvas canvas(canvasName.c_str(),
                       (std::string(labels[axis]) + " efficiency map").c_str(),
                       900, 700);
        canvas.SetRightMargin(0.15);
        canvas.SetTopMargin(0.12);
        efficiency[axis]->Draw("COLZ TEXT");
        const auto nonuniformity = CalculateNonuniformity(
            grid, axis, config.minEntriesPerBin);
        TLatex annotation;
        annotation.SetNDC();
        annotation.SetTextFont(42);
        annotation.SetTextSize(0.035);
        annotation.DrawLatex(
            0.12, 0.93,
            Form("Non-uniformity #sigma(#epsilon_{i})/<#epsilon_{i}> = %.2f%%",
                 100.0 * nonuniformity.value));
        canvas.Write();
    }
}

void WriteFakeMaps(TDirectory* directory, const Grid& grid,
                   const Config& config) {
    TDirectory::TContext context(directory);
    const std::array<const char*, kAxes> labels = {"X", "Y", "2D"};
    for (int axis = 0; axis < kAxes; ++axis) {
        TH2D map(("hFakeEfficiencyMap" + std::string(labels[axis])).c_str(),
                 ("Fake efficiency " + std::string(labels[axis]) +
                  ";Local X [mm];Local Y [mm];Fake efficiency [%]").c_str(),
                 config.xBins, config.xMin, config.xMax,
                 config.yBins, config.yMin, config.yMax);
        double maximum = 0.0;
        for (int bx = 1; bx <= config.xBins; ++bx) {
            for (int by = 1; by <= config.yBins; ++by) {
                const int flat = (bx - 1) * config.yBins + (by - 1);
                const long long total = grid.denominator[flat];
                if (total < config.minEntriesPerBin) {
                    map.SetBinContent(bx, by, -1.0);
                    continue;
                }
                const double efficiencyPercent =
                    100.0 * grid.numerator[axis][flat] / total;
                map.SetBinContent(bx, by, efficiencyPercent);
                maximum = std::max(maximum, efficiencyPercent);
            }
        }
        map.SetMinimum(0.0);
        if (maximum > 0.0) map.SetMaximum(1.1 * maximum);
        map.SetStats(false);
        map.Write();
    }
}

void WriteProjections(TDirectory* directory, const Grid& grid,
                      const Config& config) {
    TDirectory::TContext context(directory);
    const Estimate overall = Calculate(grid, k2D, config.minEntriesPerBin);

    for (int projectionAxis : {kX, kY}) {
        const bool projectX = projectionAxis == kX;
        const int bins = projectX ? config.xBins : config.yBins;
        const double minimum = projectX ? config.xMin : config.yMin;
        const double maximum = projectX ? config.xMax : config.yMax;
        const double width = (maximum - minimum) / bins;
        const std::string label = projectX ? "X" : "Y";

        TGraphAsymmErrors graph;
        graph.SetName(("gEff2DVs" + label).c_str());
        StyleGraph(graph, kBlue + 1, 20);
        int point = 0;
        for (int projectedBin = 0; projectedBin < bins; ++projectedBin) {
            long long total = 0;
            long long passed = 0;
            const int orthogonalBins = projectX ? config.yBins : config.xBins;
            for (int orthogonalBin = 0; orthogonalBin < orthogonalBins;
                 ++orthogonalBin) {
                const int bx = projectX ? projectedBin : orthogonalBin;
                const int by = projectX ? orthogonalBin : projectedBin;
                const int flat = bx * config.yBins + by;
                if (grid.denominator[flat] < config.minEntriesPerBin) continue;
                total += grid.denominator[flat];
                passed += grid.numerator[k2D][flat];
            }
            if (total == 0) continue;
            const double efficiency = static_cast<double>(passed) / total;
            const double lower = TEfficiency::ClopperPearson(
                total, passed, 0.683, false);
            const double upper = TEfficiency::ClopperPearson(
                total, passed, 0.683, true);
            graph.SetPoint(point, minimum + (projectedBin + 0.5) * width,
                           efficiency);
            graph.SetPointError(point, 0.5 * width, 0.5 * width,
                                efficiency - lower, upper - efficiency);
            ++point;
        }

        const std::string canvasName = "cEff2DProjection" + label;
        TCanvas canvas(canvasName.c_str(),
                       ("2D efficiency projection on local " + label).c_str(),
                       900, 700);
        graph.SetTitle(
            ("2D efficiency vs local " + label + ";Local " + label +
             " [mm];2D efficiency").c_str());
        graph.SetMinimum(0.0);
        graph.SetMaximum(1.02);
        graph.Draw("APL");

        TLine overallLine(minimum, overall.value, maximum, overall.value);
        overallLine.SetLineColor(kRed + 1);
        overallLine.SetLineStyle(2);
        overallLine.SetLineWidth(2);
        overallLine.Draw();

        TLegend legend(0.58, 0.18, 0.88, 0.31);
        legend.SetBorderSize(0);
        legend.AddEntry(&graph, "2D efficiency", "lp");
        const std::string averageLabel = PercentLabel(
            "Average efficiency", overall.value);
        legend.AddEntry(&overallLine, averageLabel.c_str(), "l");
        legend.Draw();
        canvas.Write();
    }
}

void StyleGraph(TGraphAsymmErrors& graph, int color, int marker) {
    graph.SetLineColor(color);
    graph.SetMarkerColor(color);
    graph.SetMarkerStyle(marker);
    graph.SetMarkerSize(0.8);
    graph.SetLineWidth(2);
}

void WriteScan(TDirectory* directory, const std::string& xTitle,
               double configuredParameter,
               const std::vector<double>& parameters,
               const std::vector<std::array<Estimate, kAxes>>& values) {
    TDirectory::TContext context(directory);
    const std::array<const char*, kAxes> labels = {"X", "Y", "2D"};
    const std::array<int, kAxes> colors = {kBlue + 1, kRed + 1, kBlack};
    const std::array<int, kAxes> markers = {20, 21, 22};
    std::array<TGraphAsymmErrors, kAxes> graphs;
    for (int axis = 0; axis < kAxes; ++axis) {
        graphs[axis].SetName(("gEfficiency" + std::string(labels[axis])).c_str());
        StyleGraph(graphs[axis], colors[axis], markers[axis]);
        for (size_t point = 0; point < parameters.size(); ++point) {
            const auto& estimate = values[point][axis];
            graphs[axis].SetPoint(point, parameters[point], estimate.value);
            graphs[axis].SetPointError(
                point, 0, 0, estimate.errorLow, estimate.errorHigh);
        }
    }

    TCanvas canvas("cEffVsMargin", "Efficiency vs margin", 900, 700);
    graphs[kX].SetTitle(
        (std::string("Efficiency vs margin;") + xTitle + ";Efficiency").c_str());
    graphs[kX].SetMinimum(0.5);
    graphs[kX].SetMaximum(1.02);
    graphs[kX].Draw("APL");
    graphs[kY].Draw("PL SAME");
    graphs[k2D].Draw("PL SAME");
    TLine line(configuredParameter, 0.0, configuredParameter, 1.02);
    line.SetLineColor(kGray + 2);
    line.SetLineStyle(2);
    line.SetLineWidth(2);
    line.Draw();
    TLegend legend(0.62, 0.18, 0.88, 0.38);
    legend.SetBorderSize(0);
    size_t configuredPoint = 0;
    if (!parameters.empty()) {
        configuredPoint = std::distance(
            parameters.begin(),
            std::min_element(parameters.begin(), parameters.end(),
                [&](double left, double right) {
                    return std::abs(left - configuredParameter) <
                           std::abs(right - configuredParameter);
                }));
    }
    for (int axis = 0; axis < kAxes; ++axis) {
        const double value = values.empty()
                                 ? 0.0
                                 : values[configuredPoint][axis].value;
        const std::string resultLabel = PercentLabel(
            std::string(labels[axis]) + " efficiency", value);
        legend.AddEntry(&graphs[axis], resultLabel.c_str(), "lp");
    }
    std::ostringstream marginLabel;
    marginLabel << "Margin: " << std::fixed << std::setprecision(2)
                << configuredParameter << " mm";
    legend.AddEntry(&line, marginLabel.str().c_str(), "l");
    legend.Draw();
    canvas.Write();
}

Grid BuildGrid(const std::vector<Sample>& samples, const Config& config,
               const std::function<std::array<bool, kAxes>(const Sample&)>& matcher) {
    Grid grid(config.xBins, config.yBins);
    for (const auto& sample : samples) grid.Fill(sample.flatBin, matcher(sample));
    return grid;
}

void FillSummary(TTree& tree, int dutID, double parameter,
                 const Grid& grid, const Config& config) {
    static Int_t outDutID;
    static Double_t outParameter, efficiencyX, efficiencyY, efficiency2D;
    static Double_t errorLowX, errorHighX, errorLowY, errorHighY, errorLow2D, errorHigh2D;
    static Double_t nonuniformityX, nonuniformityY, nonuniformity2D;
    static Long64_t totalEvents;
    static Int_t validBins;
    if (!tree.GetNbranches()) {
        tree.Branch("dutID", &outDutID);
        tree.Branch("parameterMm", &outParameter);
        tree.Branch("totalEvents", &totalEvents);
        tree.Branch("validBins", &validBins);
        tree.Branch("efficiencyX", &efficiencyX);
        tree.Branch("errorLowX", &errorLowX);
        tree.Branch("errorHighX", &errorHighX);
        tree.Branch("efficiencyY", &efficiencyY);
        tree.Branch("errorLowY", &errorLowY);
        tree.Branch("errorHighY", &errorHighY);
        tree.Branch("efficiency2D", &efficiency2D);
        tree.Branch("errorLow2D", &errorLow2D);
        tree.Branch("errorHigh2D", &errorHigh2D);
        tree.Branch("nonuniformityX", &nonuniformityX);
        tree.Branch("nonuniformityY", &nonuniformityY);
        tree.Branch("nonuniformity2D", &nonuniformity2D);
    }
    const auto x = Calculate(grid, kX, config.minEntriesPerBin);
    const auto y = Calculate(grid, kY, config.minEntriesPerBin);
    const auto xy = Calculate(grid, k2D, config.minEntriesPerBin);
    const auto nonuniformX = CalculateNonuniformity(
        grid, kX, config.minEntriesPerBin);
    const auto nonuniformY = CalculateNonuniformity(
        grid, kY, config.minEntriesPerBin);
    const auto nonuniform2D = CalculateNonuniformity(
        grid, k2D, config.minEntriesPerBin);
    outDutID = dutID;
    outParameter = parameter;
    totalEvents = xy.total;
    validBins = xy.bins;
    efficiencyX = x.value; errorLowX = x.errorLow; errorHighX = x.errorHigh;
    efficiencyY = y.value; errorLowY = y.errorLow; errorHighY = y.errorHigh;
    efficiency2D = xy.value; errorLow2D = xy.errorLow; errorHigh2D = xy.errorHigh;
    nonuniformityX = nonuniformX.value;
    nonuniformityY = nonuniformY.value;
    nonuniformity2D = nonuniform2D.value;
    tree.Fill();
}

}  // namespace

Result Analyze(const std::vector<Event>& events,
               const std::unordered_set<int>& strictSingleHitTrackerEvents,
               const std::shared_ptr<Detector>& detector,
               const Config& config,
               TDirectory* outputDirectory) {
    Result result;
    if (!detector || !outputDirectory || config.xBins <= 0 || config.yBins <= 0) return result;
    const auto* planar = detector->GetPlanarConfig();
    if (!planar) return result;
    const auto [typeX, typeY] = PlanarAxisTypes(*planar);
    if (typeX < 0 || typeY < 0 ||
        !planar->readoutPlanePitch.count(typeX) ||
        !planar->readoutPlanePitch.count(typeY)) return result;
    const double pitchX = planar->readoutPlanePitch.at(typeX);
    const double pitchY = planar->readoutPlanePitch.at(typeY);

    std::vector<Sample> samples;
    samples.reserve(events.size());
    for (const auto& event : events) {
        const auto frameIt = event.detectorFramesMap.find(detector->GetID());
        if (frameIt == event.detectorFramesMap.end()) continue;
        const auto& frame = *frameIt->second;
        const TVector3 predicted = detector->GlobalToLocal(detector->CalcHitFromTrack(event.track));
        const int flatBin = FindFlatBin(predicted.X(), predicted.Y(), config);
        if (flatBin < 0) continue;
        Sample sample;
        sample.eventID = event.eventID;
        sample.predX = predicted.X();
        sample.predY = predicted.Y();
        sample.flatBin = flatBin;
        sample.intervalsX = BuildIntervals(frame, typeX, pitchX);
        sample.intervalsY = BuildIntervals(frame, typeY, pitchY);
        samples.push_back(std::move(sample));
    }

    auto* dutDirectory = GetOrCreate(outputDirectory, "DUT_" + std::to_string(detector->GetID()));
    auto* efficiencyDirectory = GetOrCreate(dutDirectory, "Efficiency");
    auto* fakeDirectory = GetOrCreate(dutDirectory, "Fake");

    const Grid efficiency = BuildGrid(samples, config, [&](const Sample& sample) {
        return EnvelopeMatch(sample, config.margin);
    });
    WriteEventCount(efficiencyDirectory, efficiency, config);
    WriteMaps(efficiencyDirectory, efficiency, config);
    WriteProjections(efficiencyDirectory, efficiency, config);

    TDirectory::TContext dutContext(dutDirectory);
    TTree summary("EfficiencySummary", "Configured DUT efficiency results");
    FillSummary(summary, detector->GetID(), config.margin, efficiency, config);
    summary.Write();

    const auto eventWeighted = Calculate(
        efficiency, k2D, config.minEntriesPerBin);
    result.eligibleEvents = eventWeighted.total;
    result.matchedEvents = eventWeighted.passed;
    result.eventWeightedX = Calculate(
        efficiency, kX, config.minEntriesPerBin).value;
    result.eventWeightedY = Calculate(
        efficiency, kY, config.minEntriesPerBin).value;
    result.eventWeighted2D = eventWeighted.value;
    result.nonuniformityX = CalculateNonuniformity(
        efficiency, kX, config.minEntriesPerBin).value;
    result.nonuniformityY = CalculateNonuniformity(
        efficiency, kY, config.minEntriesPerBin).value;
    result.nonuniformity2D = CalculateNonuniformity(
        efficiency, k2D, config.minEntriesPerBin).value;

    auto scan = [&](TDirectory* directory, const std::string& xTitle,
                    double minimum, double maximum,
                    double step, const auto& matcher) {
        std::vector<double> parameters;
        std::vector<std::array<Estimate, kAxes>> values;
        if (step <= 0.0) return;
        for (double parameter = minimum; parameter <= maximum + 0.5 * step; parameter += step) {
            const Grid grid = BuildGrid(samples, config, [&](const Sample& sample) {
                return matcher(sample, parameter);
            });
            std::array<Estimate, kAxes> estimates{};
            for (int axis = 0; axis < kAxes; ++axis)
                estimates[axis] = Calculate(
                    grid, axis, config.minEntriesPerBin);
            parameters.push_back(parameter);
            values.push_back(estimates);
        }
        WriteScan(directory, xTitle, config.margin,
                  parameters, values);
    };
    scan(efficiencyDirectory, "Margin [mm]", config.envelopeScanMin,
         config.envelopeScanMax, config.envelopeScanStep,
         [](const Sample& sample, double margin) { return EnvelopeMatch(sample, margin); });

    Grid fakeGrid(config.xBins, config.yBins);
    std::vector<const Sample*> strict;
    for (const auto& sample : samples)
        if (strictSingleHitTrackerEvents.count(sample.eventID)) strict.push_back(&sample);
    if (config.enableFakeEfficiency && strict.size() > 1) {
        std::mt19937 random(config.fakeSeed);
        std::uniform_int_distribution<size_t> partnerDistribution(
            0, strict.size() - 2);
        const size_t partners = std::min(
            static_cast<size_t>(std::max(1, config.fakePartnersPerEvent)),
            strict.size() - 1);
        for (size_t sourceIndex = 0; sourceIndex < strict.size(); ++sourceIndex) {
            const Sample& source = *strict[sourceIndex];
            const auto remainingX = RemoveMatchedClusters(
                source.intervalsX, source.predX, config.margin);
            const auto remainingY = RemoveMatchedClusters(
                source.intervalsY, source.predY, config.margin);
            std::unordered_set<size_t> selectedPartners;
            while (selectedPartners.size() < partners) {
                size_t partnerIndex = partnerDistribution(random);
                if (partnerIndex >= sourceIndex) ++partnerIndex;
                selectedPartners.insert(partnerIndex);
            }
            for (size_t partnerIndex : selectedPartners) {
                const Sample& partner = *strict[partnerIndex];
                const bool matchedX = RequiredExtension(
                    partner.predX, remainingX) <= config.margin;
                const bool matchedY = RequiredExtension(
                    partner.predY, remainingY) <= config.margin;
                fakeGrid.Fill(partner.flatBin,
                              {matchedX, matchedY, matchedX && matchedY});
            }
        }
        WriteFakeMaps(fakeDirectory, fakeGrid, config);
        TDirectory::TContext fakeContext(fakeDirectory);
        TTree fakeSummary(
            "FakeSummary",
            "Matched-cluster-removed random-event fake efficiency");
        FillSummary(fakeSummary, detector->GetID(), config.margin,
                    fakeGrid, config);
        fakeSummary.Write();
        result.fake2D = Calculate(
            fakeGrid, k2D, config.minEntriesPerBin).value;
    }

    TDirectory::TContext metadataContext(dutDirectory);
    TTree metadata("AnalysisConfig", "DUT efficiency analysis metadata");
    Config storedConfig = config;
    Int_t strictFakeEvents = strict.size();
    metadata.Branch("xBins", &storedConfig.xBins);
    metadata.Branch("yBins", &storedConfig.yBins);
    metadata.Branch("xMin", &storedConfig.xMin);
    metadata.Branch("xMax", &storedConfig.xMax);
    metadata.Branch("yMin", &storedConfig.yMin);
    metadata.Branch("yMax", &storedConfig.yMax);
    metadata.Branch("excludedXBins", &storedConfig.excludedXBins);
    metadata.Branch("excludedYBins", &storedConfig.excludedYBins);
    metadata.Branch("margin", &storedConfig.margin);
    metadata.Branch("minEntriesPerBin", &storedConfig.minEntriesPerBin);
    metadata.Branch("envelopeScanMin", &storedConfig.envelopeScanMin);
    metadata.Branch("envelopeScanMax", &storedConfig.envelopeScanMax);
    metadata.Branch("envelopeScanStep", &storedConfig.envelopeScanStep);
    metadata.Branch("fakeSeed", &storedConfig.fakeSeed);
    metadata.Branch("fakePartnersPerEvent", &storedConfig.fakePartnersPerEvent);
    metadata.Branch("strictFakeEvents", &strictFakeEvents);
    metadata.Fill();
    metadata.Write();
    return result;
}

}  // namespace DUTEfficiency
