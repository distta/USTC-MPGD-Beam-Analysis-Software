#include "Script/OscilloscopeAnalysisScript.h"

#include "Script/Base/ScriptFactory.h"

#include <TDirectory.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using namespace std;

namespace {

constexpr int kScintillatorCount = 3;
constexpr int kEventIDBits = 16;
constexpr int kFrameBits = kEventIDBits + 2;
constexpr double kEventBitPeriodSeconds = 25.0e-9;

struct Waveform {
    vector<double> time;
    vector<double> amplitude;
};

struct TimingResult {
    double timeNs = numeric_limits<double>::quiet_NaN();
    double amplitude = numeric_limits<double>::quiet_NaN();
    double baseline = numeric_limits<double>::quiet_NaN();
    double threshold = numeric_limits<double>::quiet_NaN();
    int pulseCandidates = 0;
    bool valid = false;
};

struct EventIDResult {
    uint64_t eventID = 0;
    double frameStartNs = numeric_limits<double>::quiet_NaN();
    double threshold = numeric_limits<double>::quiet_NaN();
    bool valid = false;
};

struct GaussianFit {
    double mean = numeric_limits<double>::quiet_NaN();
    double sigma = numeric_limits<double>::quiet_NaN();
    double sigmaError = numeric_limits<double>::quiet_NaN();
    bool valid = false;
};

bool ParsePair(const string& line, double& first, double& second) {
    const char* begin = line.c_str();
    char* separator = nullptr;
    first = strtod(begin, &separator);
    if (separator == begin || *separator != ',') return false;
    const char* secondBegin = separator + 1;
    char* end = nullptr;
    second = strtod(secondBegin, &end);
    return end != secondBegin;
}

class TraceReader {
   public:
    bool Open(const filesystem::path& path, string& error) {
        m_path = path;
        m_input.open(path);
        if (!m_input) {
            error = "cannot open " + path.string();
            return false;
        }

        string line;
        if (!getline(m_input, line) || !getline(m_input, line)) {
            error = "incomplete header in " + path.string();
            return false;
        }
        replace(line.begin(), line.end(), '\r', ' ');
        vector<string> fields;
        string field;
        stringstream header(line);
        while (getline(header, field, ',')) fields.push_back(field);
        if (fields.size() < 4 || fields[0] != "Segments") {
            error = "invalid segment header in " + path.string();
            return false;
        }
        try {
            m_segmentCount = stoi(fields[1]);
            m_segmentSize = stoi(fields[3]);
        } catch (...) {
            error = "invalid segment counts in " + path.string();
            return false;
        }
        if (m_segmentCount <= 0 || m_segmentSize <= 1) {
            error = "non-positive segment dimensions in " + path.string();
            return false;
        }
        while (getline(m_input, line)) {
            if (line.rfind("Time,Ampl", 0) == 0) return true;
        }
        error = "Time,Ampl header missing in " + path.string();
        return false;
    }

    bool ReadSegment(Waveform& waveform) {
        waveform.time.clear();
        waveform.amplitude.clear();
        waveform.time.reserve(m_segmentSize);
        waveform.amplitude.reserve(m_segmentSize);
        string line;
        double firstTime = 0.0;
        double pitch = 0.0;
        for (int sample = 0; sample < m_segmentSize; ++sample) {
            if (!getline(m_input, line)) return false;
            double time = 0.0;
            double amplitude = 0.0;
            if (sample < 2) {
                if (!ParsePair(line, time, amplitude)) continue;
                if (sample == 0)
                    firstTime = time;
                else
                    pitch = time - firstTime;
            } else {
                time = firstTime + sample * pitch;
                const char* comma = strchr(line.c_str(), ',');
                if (!comma) continue;
                char* end = nullptr;
                amplitude = strtod(comma + 1, &end);
                if (end == comma + 1) continue;
            }
            waveform.time.push_back(time);
            waveform.amplitude.push_back(amplitude);
        }
        return true;
    }

    int SegmentCount() const { return m_segmentCount; }

   private:
    ifstream m_input;
    filesystem::path m_path;
    int m_segmentCount = 0;
    int m_segmentSize = 0;
};

double Quantile(vector<double> values, double fraction) {
    if (values.empty()) return numeric_limits<double>::quiet_NaN();
    const size_t index = min(
        values.size() - 1,
        static_cast<size_t>(fraction * static_cast<double>(values.size())));
    nth_element(values.begin(), values.begin() + index, values.end());
    return values[index];
}

double InterpolatedCrossing(const Waveform& waveform, size_t second,
                            double threshold) {
    if (second == 0 || second >= waveform.time.size() ||
        second >= waveform.amplitude.size())
        return numeric_limits<double>::quiet_NaN();
    const double firstAmplitude = waveform.amplitude[second - 1];
    const double secondAmplitude = waveform.amplitude[second];
    const double difference = secondAmplitude - firstAmplitude;
    if (difference == 0.0) return waveform.time[second];
    const double fraction = (threshold - firstAmplitude) / difference;
    return waveform.time[second - 1] +
           fraction * (waveform.time[second] - waveform.time[second - 1]);
}

EventIDResult DecodeEventID(const Waveform& waveform) {
    EventIDResult result;
    const size_t count =
        min(waveform.time.size(), waveform.amplitude.size());
    if (count < static_cast<size_t>(kFrameBits)) return result;

    vector<double> amplitudes(waveform.amplitude.begin(),
                              waveform.amplitude.begin() + count);
    const double low = Quantile(amplitudes, 0.05);
    const double high = Quantile(amplitudes, 0.95);
    result.threshold = 0.5 * (low + high);
    if (!isfinite(result.threshold) || high - low < 0.2) return result;

    size_t firstFalling = count;
    for (size_t sample = 1; sample < count; ++sample) {
        if (waveform.amplitude[sample - 1] > result.threshold &&
            waveform.amplitude[sample] <= result.threshold) {
            firstFalling = sample;
            break;
        }
    }
    if (firstFalling == count) return result;

    const double frameStart =
        InterpolatedCrossing(waveform, firstFalling, result.threshold);
    array<int, kFrameBits> bits{};
    for (int bit = 0; bit < kFrameBits; ++bit) {
        const double sampleTime =
            frameStart + (bit + 0.5) * kEventBitPeriodSeconds;
        auto found = lower_bound(waveform.time.begin(),
                                 waveform.time.begin() + count, sampleTime);
        if (found == waveform.time.begin() + count) return result;
        size_t sample = static_cast<size_t>(found - waveform.time.begin());
        if (sample > 0 &&
            sampleTime - waveform.time[sample - 1] <=
                waveform.time[sample] - sampleTime)
            --sample;
        bits[bit] = waveform.amplitude[sample] < result.threshold ? 1 : 0;
    }

    // C4 uses a negative-logic serial frame: start(1), 16-bit ID, stop(0).
    if (bits.front() != 1 || bits.back() != 0) return result;
    uint64_t eventID = 0;
    for (int bit = 1; bit <= kEventIDBits; ++bit)
        eventID = (eventID << 1U) | static_cast<uint64_t>(bits[bit]);
    result.eventID = eventID;
    result.frameStartNs = frameStart * 1.0e9;
    result.valid = true;
    return result;
}

TimingResult MeasureNegativePulse(const Waveform& waveform,
                                  double eventIDFrameStartNs,
                                  double cfdFraction,
                                  double minPulseAmplitude) {
    TimingResult result;
    const size_t count =
        min(waveform.time.size(), waveform.amplitude.size());
    if (count < 10 || !isfinite(eventIDFrameStartNs)) return result;
    const double eventIDFrameStart = eventIDFrameStartNs * 1.0e-9;
    const auto searchEndIterator =
        upper_bound(waveform.time.begin(), waveform.time.begin() + count,
                    eventIDFrameStart);
    const size_t searchEnd = static_cast<size_t>(
        searchEndIterator - waveform.time.begin());
    if (searchEnd < 2) return result;
    const size_t baselineSamples = min<size_t>(100, count / 4);
    vector<double> baselineValues(
        waveform.amplitude.begin(),
        waveform.amplitude.begin() + baselineSamples);
    result.baseline = Quantile(baselineValues, 0.5);
    const auto minimum = min_element(waveform.amplitude.begin(),
                                     waveform.amplitude.begin() + searchEnd);
    result.amplitude = result.baseline - *minimum;
    if (!isfinite(result.amplitude) ||
        result.amplitude < minPulseAmplitude)
        return result;
    result.threshold =
        result.baseline - cfdFraction * result.amplitude;
    // Search backward from the C4 frame and select the nearest preceding
    // scintillator pulse. Later pulses and earlier pile-up are ignored.
    size_t selectedSample = 0;
    for (size_t sample = searchEnd - 1; sample > 0; --sample) {
        if (waveform.amplitude[sample - 1] > result.threshold &&
            waveform.amplitude[sample] <= result.threshold) {
            ++result.pulseCandidates;
            if (selectedSample == 0) selectedSample = sample;
        }
    }
    if (selectedSample > 0) {
        result.timeNs =
            InterpolatedCrossing(waveform, selectedSample,
                                 result.threshold) *
            1.0e9;
        result.valid = isfinite(result.timeNs);
    }
    return result;
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

    TF1 gaussian("coreGaussian", "gaus",
                 median - 4.0 * robustSigma,
                 median + 4.0 * robustSigma);
    gaussian.SetParameters(histogram.GetMaximum(), median, robustSigma);
    if (histogram.Fit(&gaussian, "QNR") != 0) return result;
    gaussian.SetRange(gaussian.GetParameter(1) -
                          2.5 * abs(gaussian.GetParameter(2)),
                      gaussian.GetParameter(1) +
                          2.5 * abs(gaussian.GetParameter(2)));
    if (histogram.Fit(&gaussian, "QNR") != 0) return result;
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

map<int, map<int, filesystem::path>> FindTraceFiles(
    const filesystem::path& directory) {
    map<int, map<int, filesystem::path>> files;
    const regex lecroyPattern(
        R"(^C([0-9]+)--[^-]+--([0-9]+)\.csv$)");
    const regex legacyPattern(R"(^C([0-9]+)Trace([0-9]+)\.csv$)");
    for (const auto& item : filesystem::directory_iterator(directory)) {
        if (!item.is_regular_file()) continue;
        smatch match;
        const string name = item.path().filename().string();
        if (regex_match(name, match, lecroyPattern) ||
            regex_match(name, match, legacyPattern)) {
            const int channel = stoi(match[1].str());
            const int traceIndex = stoi(match[2].str());
            files[traceIndex][channel] = item.path();
        }
    }
    return files;
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
}

void OscilloscopeAnalysisScript::Print() const {
    cout << "OscilloscopeAnalysisScript: CSV=" << m_csvDirectory
         << (m_csvDirectory.empty() ? "scope/<runID> (auto)" : "")
         << ", output=" << m_outputFile
         << ", CFD=" << m_cfdFraction
         << ", minimum pulse=" << m_minPulseAmplitude << " V\n";
}

bool OscilloscopeAnalysisScript::Validate() const {
    return !m_outputFile.empty() &&
           m_cfdFraction > 0.0 && m_cfdFraction < 1.0 &&
           m_minPulseAmplitude > 0.0 && m_histogramBins > 0;
}

bool OscilloscopeAnalysisScript::Execute() {
    const filesystem::path inputDirectory =
        ResolveInputDirectory(m_csvDirectory, GetOutputDir());
    if (!filesystem::is_directory(inputDirectory)) {
        cerr << "[OscilloscopeAnalysis] CSV directory does not exist: "
             << inputDirectory << '\n';
        return false;
    }
    auto traceFiles = FindTraceFiles(inputDirectory);
    vector<int> traceIndices;
    for (const auto& [index, channels] : traceFiles) {
        bool complete = true;
        for (int channel = 1; channel <= 4; ++channel)
            complete = complete && channels.count(channel) > 0;
        if (complete) traceIndices.push_back(index);
    }
    if (m_maxWaveformFiles > 0 &&
        traceIndices.size() > static_cast<size_t>(m_maxWaveformFiles))
        traceIndices.resize(static_cast<size_t>(m_maxWaveformFiles));
    if (traceIndices.empty()) {
        cerr << "[OscilloscopeAnalysis] no complete C1-C4 trace sets in "
             << inputDirectory << '\n';
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

    TTree events("OscilloscopeEvents",
                 "Decoded event ID and scintillator CFD times");
    ULong64_t eventID = 0;
    Int_t traceFileIndex = -1, segmentIndex = -1;
    Bool_t eventIDValid = false;
    Double_t eventIDFrameStart = numeric_limits<double>::quiet_NaN();
    Double_t eventIDThreshold = numeric_limits<double>::quiet_NaN();
    array<Bool_t, kScintillatorCount> timeValid{};
    array<Double_t, kScintillatorCount> time{};
    array<Double_t, kScintillatorCount> amplitude{};
    array<Double_t, kScintillatorCount> baseline{};
    array<Double_t, kScintillatorCount> threshold{};
    array<Int_t, kScintillatorCount> pulseCandidates{};
    events.Branch("eventID", &eventID);
    events.Branch("eventIDValid", &eventIDValid);
    events.Branch("traceFileIndex", &traceFileIndex);
    events.Branch("segmentIndex", &segmentIndex);
    events.Branch("eventIDFrameStart", &eventIDFrameStart);
    events.Branch("eventIDThreshold", &eventIDThreshold);
    for (int channel = 0; channel < kScintillatorCount; ++channel) {
        const string prefix = "C" + to_string(channel + 1);
        events.Branch((prefix + "TimeValid").c_str(), &timeValid[channel]);
        events.Branch((prefix + "Time").c_str(), &time[channel]);
        events.Branch((prefix + "Amplitude").c_str(),
                      &amplitude[channel]);
        events.Branch((prefix + "Baseline").c_str(), &baseline[channel]);
        events.Branch((prefix + "CFDThreshold").c_str(),
                      &threshold[channel]);
        events.Branch((prefix + "PulseCandidates").c_str(),
                      &pulseCandidates[channel]);
    }

    array<vector<double>, 3> pairDifferences;
    set<uint64_t> uniqueEventIDs;
    size_t decodedEventIDs = 0;
    size_t invalidEventIDs = 0;
    array<size_t, kScintillatorCount> invalidTimes{};
    array<size_t, kScintillatorCount> multiplePulseEvents{};
    size_t processedFiles = 0;
    for (int index : traceIndices) {
        array<TraceReader, 4> readers;
        string error;
        bool opened = true;
        for (int channel = 0; channel < 4; ++channel)
            opened = opened &&
                     readers[channel].Open(
                         traceFiles[index][channel + 1], error);
        if (!opened) {
            cerr << "[OscilloscopeAnalysis] " << error << '\n';
            return false;
        }
        int segments = readers.front().SegmentCount();
        for (const TraceReader& reader : readers)
            segments = min(segments, reader.SegmentCount());
        array<Waveform, 4> waveforms;
        for (int segment = 0; segment < segments; ++segment) {
            for (size_t channel = 0; channel < readers.size(); ++channel) {
                if (!readers[channel].ReadSegment(waveforms[channel])) {
                    cerr << "[OscilloscopeAnalysis] incomplete segment "
                         << segment << " in trace " << index << '\n';
                    return false;
                }
            }
            const EventIDResult decoded = DecodeEventID(waveforms[3]);
            eventID = decoded.eventID;
            eventIDValid = decoded.valid;
            eventIDFrameStart = decoded.frameStartNs;
            eventIDThreshold = decoded.threshold;
            traceFileIndex = index;
            segmentIndex = segment;
            if (decoded.valid) {
                ++decodedEventIDs;
                uniqueEventIDs.insert(decoded.eventID);
            } else {
                ++invalidEventIDs;
            }

            array<TimingResult, kScintillatorCount> timing;
            for (int channel = 0; channel < kScintillatorCount;
                 ++channel) {
                timing[channel] =
                    MeasureNegativePulse(waveforms[channel],
                                         decoded.frameStartNs,
                                         m_cfdFraction,
                                         m_minPulseAmplitude);
                timeValid[channel] = timing[channel].valid;
                time[channel] = timing[channel].timeNs;
                amplitude[channel] = timing[channel].amplitude;
                baseline[channel] = timing[channel].baseline;
                threshold[channel] = timing[channel].threshold;
                pulseCandidates[channel] =
                    timing[channel].pulseCandidates;
                if (timing[channel].pulseCandidates > 1)
                    ++multiplePulseEvents[channel];
                if (!timing[channel].valid)
                    ++invalidTimes[channel];
            }
            if (timing[0].valid && timing[1].valid)
                pairDifferences[0].push_back(
                    timing[0].timeNs - timing[1].timeNs);
            if (timing[0].valid && timing[2].valid)
                pairDifferences[1].push_back(
                    timing[0].timeNs - timing[2].timeNs);
            if (timing[1].valid && timing[2].valid)
                pairDifferences[2].push_back(
                    timing[1].timeNs - timing[2].timeNs);
            events.Fill();
        }
        ++processedFiles;
        if (processedFiles % 25 == 0 ||
            processedFiles == traceIndices.size())
            cout << "\r[OscilloscopeAnalysis] trace files "
                 << processedFiles << '/' << traceIndices.size() << flush;
    }
    cout << '\n';
    events.Write();

    TDirectory* timingDirectory = output->mkdir("Timing");
    timingDirectory->cd();
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

    output->cd();
    TTree summary("TimingResolution",
                  "Three-scintillator timing resolution summary");
    array<Double_t, 3> pairSigma{}, pairSigmaError{};
    array<Double_t, kScintillatorCount> resolution{}, resolutionError{};
    array<Bool_t, kScintillatorCount> resolutionValid{};
    ULong64_t decodedEvents = decodedEventIDs;
    ULong64_t uniqueEvents = uniqueEventIDs.size();
    ULong64_t invalidEvents = invalidEventIDs;
    for (size_t pairIndex = 0; pairIndex < pairFits.size(); ++pairIndex) {
        pairSigma[pairIndex] = pairFits[pairIndex].sigma;
        pairSigmaError[pairIndex] = pairFits[pairIndex].sigmaError;
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
    for (int channel = 0; channel < kScintillatorCount; ++channel) {
        resolutionValid[channel] =
            all_of(pairFits.begin(), pairFits.end(),
                   [](const GaussianFit& fit) { return fit.valid; }) &&
            variance[channel] > 0.0;
        resolution[channel] =
            resolutionValid[channel]
                ? sqrt(variance[channel])
                : numeric_limits<double>::quiet_NaN();
        if (resolutionValid[channel]) {
            array<size_t, 3> terms{};
            if (channel == 0)
                terms = {0, 1, 2};
            else if (channel == 1)
                terms = {0, 2, 1};
            else
                terms = {1, 2, 0};
            double varianceErrorSquared = 0.0;
            for (size_t term : terms) {
                const double contribution =
                    pairSigma[term] * pairSigmaError[term];
                varianceErrorSquared += contribution * contribution;
            }
            resolutionError[channel] =
                sqrt(varianceErrorSquared) /
                (2.0 * resolution[channel]);
        } else {
            resolutionError[channel] =
                numeric_limits<double>::quiet_NaN();
        }
    }
    summary.Branch("decodedEvents", &decodedEvents);
    summary.Branch("uniqueEventIDs", &uniqueEvents);
    summary.Branch("invalidEventIDs", &invalidEvents);
    summary.Branch("sigmaC1C2", &pairSigma[0]);
    summary.Branch("sigmaC1C3", &pairSigma[1]);
    summary.Branch("sigmaC2C3", &pairSigma[2]);
    summary.Branch("sigmaErrorC1C2", &pairSigmaError[0]);
    summary.Branch("sigmaErrorC1C3", &pairSigmaError[1]);
    summary.Branch("sigmaErrorC2C3", &pairSigmaError[2]);
    summary.Branch("resolutionC1", &resolution[0]);
    summary.Branch("resolutionC2", &resolution[1]);
    summary.Branch("resolutionC3", &resolution[2]);
    summary.Branch("resolutionErrorC1", &resolutionError[0]);
    summary.Branch("resolutionErrorC2", &resolutionError[1]);
    summary.Branch("resolutionErrorC3", &resolutionError[2]);
    summary.Branch("resolutionC1Valid", &resolutionValid[0]);
    summary.Branch("resolutionC2Valid", &resolutionValid[1]);
    summary.Branch("resolutionC3Valid", &resolutionValid[2]);
    summary.Fill();
    summary.Write();
    output->Close();

    cout << "[OscilloscopeAnalysis] event IDs=" << decodedEventIDs
         << ", unique IDs=" << uniqueEventIDs.size()
         << ", invalid IDs=" << invalidEventIDs;
    for (int channel = 0; channel < kScintillatorCount; ++channel)
        cout << ", invalid C" << channel + 1 << " times="
             << invalidTimes[channel]
             << ", multiple C" << channel + 1 << " pulses="
             << multiplePulseEvents[channel];
    cout << '\n';
    cout << "[OscilloscopeAnalysis] pair sigmas: C1-C2="
         << pairSigma[0] << " ns, C1-C3=" << pairSigma[1]
         << " ns, C2-C3=" << pairSigma[2] << " ns\n";
    cout << "[OscilloscopeAnalysis] resolutions: C1="
         << resolution[0] << " ns, C2=" << resolution[1]
         << " ns, C3=" << resolution[2] << " ns\n";
    cout << "[OscilloscopeAnalysis] output=" << outputPath << '\n';
    return decodedEventIDs > 0;
}

REGISTER_SCRIPT("OscilloscopeAnalysis", OscilloscopeAnalysisScript);
