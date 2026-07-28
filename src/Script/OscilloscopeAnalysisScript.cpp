#include "Script/OscilloscopeAnalysisScript.h"

#include "Algorithm/Oscilloscope/OscilloscopeDataProcessor.h"
#include "Script/Base/ScriptFactory.h"

#include <TError.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

using namespace std;

namespace {

constexpr size_t kScintillatorCount =
    OscilloscopeDataProcessor::kChannelCount;

struct GaussianFit {
    double mean = numeric_limits<double>::quiet_NaN();
    double sigma = numeric_limits<double>::quiet_NaN();
    double sigmaError = numeric_limits<double>::quiet_NaN();
    bool valid = false;
};

double Quantile(vector<double> values, double fraction) {
    if (values.empty()) return numeric_limits<double>::quiet_NaN();
    const size_t index = min(
        values.size() - 1,
        static_cast<size_t>(fraction * static_cast<double>(values.size())));
    nth_element(values.begin(), values.begin() + index, values.end());
    return values[index];
}

GaussianFit FitCoreGaussian(TH1D& histogram,
                            const vector<double>& values) {
    GaussianFit result;
    if (values.size() < 20) return result;
    vector<double> sorted = values;
    const double median = Quantile(sorted, 0.5);
    vector<double> deviations;
    deviations.reserve(values.size());
    for (double value : values) deviations.push_back(abs(value - median));
    const double robustSigma = 1.4826 * Quantile(deviations, 0.5);
    if (!isfinite(robustSigma) || robustSigma <= 0.0) return result;

    TF1 gaussian(("fit_" + string(histogram.GetName())).c_str(), "gaus",
                 median - 4.0 * robustSigma,
                 median + 4.0 * robustSigma);
    gaussian.SetParameters(histogram.GetMaximum(), median, robustSigma);
    auto quietFit = [&](const char* options) {
        const int previousLevel = gErrorIgnoreLevel;
        gErrorIgnoreLevel = max(gErrorIgnoreLevel, static_cast<int>(kWarning));
        const int status = histogram.Fit(&gaussian, options);
        gErrorIgnoreLevel = previousLevel;
        return status;
    };
    if (quietFit("QNR") != 0) return result;
    gaussian.SetRange(gaussian.GetParameter(1) -
                          2.5 * abs(gaussian.GetParameter(2)),
                      gaussian.GetParameter(1) +
                          2.5 * abs(gaussian.GetParameter(2)));
    if (quietFit("QR") != 0) return result;
    result.mean = gaussian.GetParameter(1);
    result.sigma = abs(gaussian.GetParameter(2));
    result.sigmaError = gaussian.GetParError(2);
    result.valid = isfinite(result.sigma) && result.sigma > 0.0;
    return result;
}

filesystem::path ResolveInputDirectory(const string& configured,
                                       const string& outputDirectory) {
    filesystem::path output =
        filesystem::path(outputDirectory).lexically_normal();
    if (output.filename().empty()) output = output.parent_path();
    const string runID = output.filename().string();
    const filesystem::path base =
        output.parent_path().parent_path();
    if (configured.empty()) return base / "scope" / runID;

    filesystem::path path(configured);
    if (path.is_absolute() || filesystem::is_directory(path)) return path;
    return base / path;
}

}  // namespace

void OscilloscopeAnalysisScript::LoadConfig(const json& config) {
    m_csvDirectory = config.value("csvDirectory", m_csvDirectory);
    m_outputFile = config.value("outputFile", m_outputFile);
    m_cfdFraction = config.value("cfdFraction", m_cfdFraction);
    m_minPulseAmplitude =
        config.value("minPulseAmplitude", m_minPulseAmplitude);
    m_histogramBins = config.value("histogramBins", m_histogramBins);
    m_maxWaveformFiles =
        config.value("maxWaveformFiles", m_maxWaveformFiles);
    m_eventIDTimeWindowNs =
        config.value("eventIDTimeWindowNs", m_eventIDTimeWindowNs);
    m_eventIDMedianFilterSamples = config.value(
        "eventIDMedianFilterSamples", m_eventIDMedianFilterSamples);
    m_signalTimeWindowNs =
        config.value("signalTimeWindowNs", m_signalTimeWindowNs);
}

void OscilloscopeAnalysisScript::Print() const {
    cout << "OscilloscopeAnalysisScript: CSV=" << m_csvDirectory
         << (m_csvDirectory.empty() ? "scope/<runID> (auto)" : "")
         << ", output=" << m_outputFile
         << ", CFD=" << m_cfdFraction
         << ", minimum pulse=" << m_minPulseAmplitude << " V"
         << ", event-ID window=[" << m_eventIDTimeWindowNs[0] << ", "
         << m_eventIDTimeWindowNs[1] << "] ns"
         << ", event-ID median filter=" << m_eventIDMedianFilterSamples
         << " samples"
         << ", signal window=[" << m_signalTimeWindowNs[0] << ", "
         << m_signalTimeWindowNs[1] << "] ns\n";
}

bool OscilloscopeAnalysisScript::Validate() const {
    return !m_outputFile.empty() &&
           m_cfdFraction > 0.0 && m_cfdFraction < 1.0 &&
           m_minPulseAmplitude > 0.0 && m_histogramBins > 0 &&
           isfinite(m_eventIDTimeWindowNs[0]) &&
           isfinite(m_eventIDTimeWindowNs[1]) &&
           m_eventIDTimeWindowNs[0] < m_eventIDTimeWindowNs[1] &&
           m_eventIDMedianFilterSamples > 0 &&
           m_eventIDMedianFilterSamples % 2 == 1 &&
           isfinite(m_signalTimeWindowNs[0]) &&
           isfinite(m_signalTimeWindowNs[1]) &&
           m_signalTimeWindowNs[0] < m_signalTimeWindowNs[1] &&
           m_signalTimeWindowNs[1] <= 0.0;
}

bool OscilloscopeAnalysisScript::Execute() {
    const filesystem::path inputDirectory =
        ResolveInputDirectory(m_csvDirectory, GetOutputDir());

    OscilloscopeProcessingConfig processingConfig;
    processingConfig.cfdFraction = m_cfdFraction;
    processingConfig.minPulseAmplitude = m_minPulseAmplitude;
    processingConfig.maxWaveformFiles = m_maxWaveformFiles;
    processingConfig.eventIDTimeWindowNs = m_eventIDTimeWindowNs;
    processingConfig.eventIDMedianFilterSamples =
        m_eventIDMedianFilterSamples;
    processingConfig.signalTimeWindowNs = m_signalTimeWindowNs;

    OscilloscopeDataProcessor processor;
    const auto printInput = [&](size_t traceFiles) {
        if (traceFiles == 0) return;
        cout << "\n[OscilloscopeAnalysis] Input\n"
             << "  Directory               : " << inputDirectory << '\n'
             << "  Complete C1-C4 traces   : " << traceFiles << '\n'
             << "  Event-ID search window  : ["
             << m_eventIDTimeWindowNs[0] << ", "
             << m_eventIDTimeWindowNs[1] << "] ns\n"
             << "  Event-ID median filter  : "
             << m_eventIDMedianFilterSamples << " samples\n"
             << "  CFD fraction            : " << m_cfdFraction << '\n'
             << "  Minimum pulse amplitude : "
             << m_minPulseAmplitude << " V\n"
             << "  Signal-time window      : ["
             << m_signalTimeWindowNs[0] << ", "
             << m_signalTimeWindowNs[1]
             << "] ns relative to eventIDTime\n";
    };
    const auto printProgress = [](size_t processed, size_t total) {
        if (processed % 25 == 0 || processed == total)
            cout << "\r[OscilloscopeAnalysis] trace files "
                 << processed << '/' << total << flush;
    };
    const OscilloscopeProcessingResult processing = processor.Process(
        inputDirectory, processingConfig, printInput, printProgress);
    if (processing.traceFiles > 0) cout << '\n';
    if (!processing.success) {
        cerr << "[OscilloscopeAnalysis] " << processing.error << '\n';
        return false;
    }

    const filesystem::path outputPath =
        filesystem::path(m_outputFile).is_absolute()
            ? filesystem::path(m_outputFile)
            : filesystem::path(GetOutputDir()) / m_outputFile;
    unique_ptr<TFile> output(TFile::Open(outputPath.c_str(), "RECREATE"));
    if (!output || output->IsZombie()) {
        cerr << "[OscilloscopeAnalysis] cannot create " << outputPath
             << '\n';
        return false;
    }

    const auto& outputEvents = processing.events;
    const auto& pairDifferences = processing.pairDifferences;

    output->cd();
    array<GaussianFit, 3> pairFits;
    const array<pair<int, int>, 3> pairs = {
        pair<int, int>{1, 2}, {1, 3}, {2, 3}};
    for (size_t pairIndex = 0; pairIndex < pairs.size(); ++pairIndex) {
        const auto& values = pairDifferences[pairIndex];
        if (values.empty()) continue;
        const double median = Quantile(values, 0.5);
        vector<double> deviations;
        deviations.reserve(values.size());
        for (double value : values)
            deviations.push_back(abs(value - median));
        const double robustSigma =
            max(0.05, 1.4826 * Quantile(deviations, 0.5));
        const string name =
            "hDeltaT_C" + to_string(pairs[pairIndex].first) + "_C" +
            to_string(pairs[pairIndex].second);
        const string title =
            "C" + to_string(pairs[pairIndex].first) + " - C" +
            to_string(pairs[pairIndex].second) +
            ";#Deltat [ns];Events";
        TH1D histogram(name.c_str(), title.c_str(), m_histogramBins,
                       median - 8.0 * robustSigma,
                       median + 8.0 * robustSigma);
        for (double value : values) histogram.Fill(value);
        pairFits[pairIndex] = FitCoreGaussian(histogram, values);
        histogram.Write();
    }

    array<double, 3> pairSigma{};
    for (size_t pairIndex = 0; pairIndex < pairFits.size(); ++pairIndex) {
        pairSigma[pairIndex] = pairFits[pairIndex].sigma;
    }
    const array<double, kScintillatorCount> variance = {
        0.5 * (pairSigma[0] * pairSigma[0] +
               pairSigma[1] * pairSigma[1] -
               pairSigma[2] * pairSigma[2]),
        0.5 * (pairSigma[0] * pairSigma[0] +
               pairSigma[2] * pairSigma[2] -
               pairSigma[1] * pairSigma[1]),
        0.5 * (pairSigma[1] * pairSigma[1] +
               pairSigma[2] * pairSigma[2] -
               pairSigma[0] * pairSigma[0])};
    const bool resolutionValid =
        all_of(pairFits.begin(), pairFits.end(),
               [](const GaussianFit& fit) { return fit.valid; }) &&
        all_of(variance.begin(), variance.end(),
               [](double value) {
                   return isfinite(value) && value > 0.0;
               });
    if (!resolutionValid) {
        cerr << "[OscilloscopeAnalysis] cannot derive three positive channel "
                "resolutions; referenceTime cannot be calculated\n";
        output->Close();
        return false;
    }

    array<double, kScintillatorCount> resolution{};
    array<double, kScintillatorCount> weight{};
    double inverseVarianceSum = 0.0;
    for (size_t channel = 0; channel < resolution.size(); ++channel) {
        resolution[channel] = sqrt(variance[channel]);
        inverseVarianceSum += 1.0 / variance[channel];
    }
    for (size_t channel = 0; channel < weight.size(); ++channel)
        weight[channel] =
            (1.0 / variance[channel]) / inverseVarianceSum;
    const double referenceResolution =
        1.0 / sqrt(inverseVarianceSum);

    TTree events("Events", "Oscilloscope event times relative to eventIDTime");
    ULong64_t eventID = 0;
    Double_t eventIDTime = numeric_limits<double>::quiet_NaN();
    Double_t referenceTime = numeric_limits<double>::quiet_NaN();
    array<Double_t, kScintillatorCount> time{};
    array<Double_t, kScintillatorCount> amplitude{};
    events.Branch("eventID", &eventID);
    events.Branch("eventIDTime", &eventIDTime);
    events.Branch("referenceTime", &referenceTime);
    for (size_t channel = 0; channel < time.size(); ++channel) {
        const string prefix = "C" + to_string(channel + 1);
        events.Branch((prefix + "Time").c_str(), &time[channel]);
        events.Branch((prefix + "Amp").c_str(), &amplitude[channel]);
    }

    vector<double> referenceTimes;
    referenceTimes.reserve(outputEvents.size());
    for (const OscilloscopeEvent& event : outputEvents) {
        eventID = event.eventID;
        eventIDTime = event.eventIDTime;
        referenceTime = 0.0;
        for (size_t channel = 0; channel < time.size(); ++channel) {
            time[channel] = event.time[channel];
            amplitude[channel] = event.amplitude[channel];
            referenceTime += weight[channel] * time[channel];
        }
        referenceTimes.push_back(referenceTime);
        events.Fill();
    }
    output->cd();
    events.Write();

    const double referenceMedian = Quantile(referenceTimes, 0.5);
    vector<double> referenceDeviations;
    referenceDeviations.reserve(referenceTimes.size());
    for (double value : referenceTimes)
        referenceDeviations.push_back(abs(value - referenceMedian));
    const double referenceRobustSigma =
        max(0.05, 1.4826 * Quantile(referenceDeviations, 0.5));
    TH1D referenceHistogram(
        "hReferenceTime",
        "Inverse-variance weighted reference time;"
        "Reference time relative to eventIDTime [ns];Events",
        m_histogramBins,
        referenceMedian - 8.0 * referenceRobustSigma,
        referenceMedian + 8.0 * referenceRobustSigma);
    for (double value : referenceTimes)
        referenceHistogram.Fill(value);
    referenceHistogram.Write();
    output->Close();

    cout << "\n[OscilloscopeAnalysis] Event summary\n"
         << "  Processed segments   : " << processing.processedSegments << '\n'
         << "  Decoded event IDs    : " << processing.decodedEventIDs << '\n'
         << "  Unique event IDs     : " << processing.uniqueEventIDs.size()
         << '\n'
         << "  Invalid event IDs    : " << processing.invalidEventIDs << '\n'
         << "  Complete Events rows : " << outputEvents.size() << '\n';

    cout << "\n[OscilloscopeAnalysis] Channel selection\n"
         << "  Times are relative to eventIDTime; selection window is ["
         << m_signalTimeWindowNs[0] << ", " << m_signalTimeWindowNs[1]
         << "] ns.\n"
         << "  " << left << setw(8) << "Channel"
         << right << setw(10) << "Valid"
         << setw(12) << "Low amp"
         << setw(14) << "No crossing"
         << setw(12) << "No window" << '\n';
    for (size_t channel = 0; channel < kScintillatorCount; ++channel) {
        const auto& statistics =
            processing.channelStatistics[channel];
        cout << "  " << left << setw(8)
             << ("C" + to_string(channel + 1))
             << right << setw(10) << statistics.valid
             << setw(12) << statistics.belowAmplitude
             << setw(14) << statistics.noCrossing
             << setw(12) << statistics.noWindow << '\n';
    }

    cout << fixed << setprecision(4)
         << "\n[OscilloscopeAnalysis] Pair Gaussian fits\n"
         << "  " << left << setw(10) << "Pair"
         << right << setw(10) << "Entries"
         << setw(14) << "Mean [ns]"
         << setw(15) << "Sigma [ns]"
         << setw(15) << "Error [ns]" << '\n';
    for (size_t pair = 0; pair < pairs.size(); ++pair) {
        const string label =
            "C" + to_string(pairs[pair].first) + "-C" +
            to_string(pairs[pair].second);
        cout << "  " << left << setw(10) << label
             << right << setw(10) << pairDifferences[pair].size()
             << setw(14) << pairFits[pair].mean
             << setw(15) << pairFits[pair].sigma
             << setw(15) << pairFits[pair].sigmaError << '\n';
    }

    cout << "\n[OscilloscopeAnalysis] Reference time\n"
         << "  Channel resolutions [ns]: C1=" << resolution[0]
         << ", C2=" << resolution[1] << ", C3=" << resolution[2] << '\n'
         << "  Channel weights         : C1=" << weight[0]
         << ", C2=" << weight[1] << ", C3=" << weight[2] << '\n'
         << "  Combined resolution     : " << referenceResolution
         << " ns\n"
         << "  Stored channel time     : C_iTime = raw CFD time - "
            "eventIDTime\n"
         << "  Definition              : referenceTime = "
            "sum_i(w_i * C_iTime)\n"
         << "\n[OscilloscopeAnalysis] Output\n"
         << "  File    : " << outputPath << '\n'
         << "  Objects : Events, hDeltaT_C1_C2, hDeltaT_C1_C3, "
            "hDeltaT_C2_C3, hReferenceTime\n"
         << defaultfloat;
    return !outputEvents.empty();
}

REGISTER_SCRIPT("OscilloscopeAnalysis", OscilloscopeAnalysisScript);
