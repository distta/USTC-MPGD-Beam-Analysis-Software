#include "Script/TimeResolutionScript.h"

#include "Algorithm/algorithms/WaveformProcessor.h"
#include "Event/DataModel.h"
#include "Script/Base/RawDataParser.h"
#include "Script/Base/ScriptFactory.h"

#include <TFile.h>
#include <TDirectory.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLine.h>
#include <TNamed.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <future>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace std;

namespace {

constexpr size_t kEstimatorCount = 6;
constexpr int kTimingCacheVersion = 5;
constexpr double kEventBitPeriodSeconds = 25.0e-9;
const array<string, kEstimatorCount> kEstimatorNames = {
    "XMean", "XMin", "YMean", "YMin", "XYMean", "XYMin"};

struct WaveformSegment {
    vector<double> time;
    vector<double> amplitude;
};

struct TraceFile {
    int segmentCount = 0;
    int segmentSize = 0;
    vector<WaveformSegment> segments;
};

struct ReferenceTime {
    uint64_t decodedEventID = 0;
    double triggerTimeNs = 0.0;
    double firstC3TransitionTimeNs = 0.0;
    double truthT0Ns = 0.0;
    int fileIndex = -1;
    int segmentIndex = -1;
};

struct DecodeDiagnostics {
    vector<double> eventSampleTimesNs;
    vector<double> eventSampleAmplitudes;
    string decodedBits;
};

using TrackKey = pair<uint64_t, int>;

struct DetectorTimes {
    array<double, kEstimatorCount> value{};
    DetectorTimes() { value.fill(numeric_limits<double>::quiet_NaN()); }
};

struct FitResult {
    long long entries = 0;
    double mean = numeric_limits<double>::quiet_NaN();
    double sigma = numeric_limits<double>::quiet_NaN();
    double sigmaError = numeric_limits<double>::quiet_NaN();
};

struct TimeSample {
    double value = 0.0;
    double t0 = numeric_limits<double>::quiet_NaN();
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

bool ReadTraceFile(const filesystem::path& path, TraceFile& trace, string& error,
                   double timeMin, double timeMax) {
    ifstream input(path);
    if (!input) {
        error = "cannot open " + path.string();
        return false;
    }

    string line;
    if (!getline(input, line) || !getline(input, line)) {
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
        trace.segmentCount = stoi(fields[1]);
        trace.segmentSize = stoi(fields[3]);
    } catch (...) {
        error = "invalid segment counts in " + path.string();
        return false;
    }
    if (trace.segmentCount <= 0 || trace.segmentSize <= 1) {
        error = "non-positive segment dimensions in " + path.string();
        return false;
    }

    bool foundDataHeader = false;
    while (getline(input, line)) {
        if (line.rfind("Time,Ampl", 0) == 0) {
            foundDataHeader = true;
            break;
        }
    }
    if (!foundDataHeader) {
        error = "Time,Ampl header missing in " + path.string();
        return false;
    }

    trace.segments.assign(trace.segmentCount, {});
    for (auto& segment : trace.segments) {
        segment.time.reserve(trace.segmentSize / 2);
        segment.amplitude.reserve(trace.segmentSize / 2);
    }

    const long long expected = static_cast<long long>(trace.segmentCount) * trace.segmentSize;
    long long sample = 0;
    while (sample < expected && getline(input, line)) {
        const int segment = static_cast<int>(sample / trace.segmentSize);
        const int localSample = static_cast<int>(sample % trace.segmentSize);
        static thread_local vector<double> firstTimes;
        if (firstTimes.size() != static_cast<size_t>(trace.segmentCount * 2))
            firstTimes.assign(static_cast<size_t>(trace.segmentCount * 2), 0.0);

        if (localSample < 2) {
            double time = 0.0, amplitude = 0.0;
            if (ParsePair(line, time, amplitude)) {
                firstTimes[static_cast<size_t>(segment * 2 + localSample)] = time;
                if (time >= timeMin && time <= timeMax) {
                    trace.segments[segment].time.push_back(time);
                    trace.segments[segment].amplitude.push_back(amplitude);
                }
            }
        } else {
            const double first = firstTimes[static_cast<size_t>(segment * 2)];
            const double pitch = firstTimes[static_cast<size_t>(segment * 2 + 1)] - first;
            const double time = first + localSample * pitch;
            if (time >= timeMin && time <= timeMax) {
                const char* comma = strchr(line.c_str(), ',');
                if (comma) {
                    char* end = nullptr;
                    const double amplitude = strtod(comma + 1, &end);
                    if (end != comma + 1) {
                        trace.segments[segment].time.push_back(time);
                        trace.segments[segment].amplitude.push_back(amplitude);
                    }
                }
            }
        }
        ++sample;
    }
    if (sample != expected) {
        error = "expected " + to_string(expected) + " samples but read " +
                to_string(sample) + " from " + path.string();
        return false;
    }
    return true;
}

vector<double> RisingCrossings(const WaveformSegment& waveform, double threshold) {
    vector<double> result;
    const size_t count = min(waveform.time.size(), waveform.amplitude.size());
    for (size_t i = 1; i < count; ++i) {
        const double y0 = waveform.amplitude[i - 1];
        const double y1 = waveform.amplitude[i];
        if (!(y0 < threshold && y1 >= threshold) || y1 == y0) continue;
        const double fraction = (threshold - y0) / (y1 - y0);
        result.push_back(waveform.time[i - 1] +
                         fraction * (waveform.time[i] - waveform.time[i - 1]));
    }
    return result;
}

double FirstLevelTransitionAfter(const WaveformSegment& waveform, double threshold,
                                 double startTime) {
    const size_t count = min(waveform.time.size(), waveform.amplitude.size());
    if (count < 2) return numeric_limits<double>::quiet_NaN();
    auto timeIt = lower_bound(waveform.time.begin(), waveform.time.end(), startTime);
    size_t begin = timeIt == waveform.time.begin() ? 1 : static_cast<size_t>(timeIt - waveform.time.begin());
    for (size_t i = begin; i < count; ++i) {
        const double y0 = waveform.amplitude[i - 1];
        const double y1 = waveform.amplitude[i];
        const bool rising = y0 < threshold && y1 >= threshold;
        const bool falling = y0 >= threshold && y1 < threshold;
        if ((!rising && !falling) || y1 == y0) continue;
        const double fraction = (threshold - y0) / (y1 - y0);
        const double crossing = waveform.time[i - 1] +
                                fraction * (waveform.time[i] - waveform.time[i - 1]);
        if (crossing >= startTime) return crossing;
    }
    return numeric_limits<double>::quiet_NaN();
}

double InterpolateAmplitude(const WaveformSegment& waveform, double time) {
    if (waveform.time.empty() || time < waveform.time.front() || time > waveform.time.back())
        return numeric_limits<double>::quiet_NaN();
    auto upper = lower_bound(waveform.time.begin(), waveform.time.end(), time);
    if (upper == waveform.time.begin()) return waveform.amplitude.front();
    if (upper == waveform.time.end()) return waveform.amplitude.back();
    const size_t i = static_cast<size_t>(upper - waveform.time.begin());
    const double t0 = waveform.time[i - 1], t1 = waveform.time[i];
    if (t1 == t0) return waveform.amplitude[i];
    const double fraction = (time - t0) / (t1 - t0);
    return waveform.amplitude[i - 1] + fraction * (waveform.amplitude[i] - waveform.amplitude[i - 1]);
}

bool DecodeReference(const WaveformSegment& eventCode, const WaveformSegment& trigger,
                     double triggerThreshold, double eventThreshold, int eventBits,
                     uint64_t expectedEventID, ReferenceTime& reference,
                     DecodeDiagnostics* diagnostics = nullptr) {
    const auto triggerCrossings = RisingCrossings(trigger, triggerThreshold);
    if (triggerCrossings.empty()) return false;
    const double triggerTime = triggerCrossings.front();
    const double firstEventTransition =
        FirstLevelTransitionAfter(eventCode, eventThreshold, triggerTime);
    if (!isfinite(firstEventTransition)) return false;

    uint64_t decodedEventID = 0;
    if (diagnostics) {
        diagnostics->eventSampleTimesNs.clear();
        diagnostics->eventSampleAmplitudes.clear();
        diagnostics->decodedBits.clear();
    }
    for (int bit = 0; bit < eventBits; ++bit) {
        const double sampleTime = firstEventTransition + (bit + 0.5) * kEventBitPeriodSeconds;
        const double amplitude = InterpolateAmplitude(eventCode, sampleTime);
        if (!isfinite(amplitude)) return false;
        const bool value = amplitude >= eventThreshold;
        decodedEventID = (decodedEventID << 1U) | static_cast<uint64_t>(value);
        if (diagnostics) {
            diagnostics->eventSampleTimesNs.push_back(sampleTime * 1.0e9);
            diagnostics->eventSampleAmplitudes.push_back(amplitude);
            diagnostics->decodedBits.push_back(value ? '1' : '0');
        }
    }
    // A malformed C3 cell must not move the running counter far away and
    // poison every later segment. Genuine acquisition gaps here are small;
    // larger discontinuities are handled as one invalid segment.
    constexpr uint64_t kMaxCounterJump = 32;
    const uint64_t distance = decodedEventID > expectedEventID
                                  ? decodedEventID - expectedEventID
                                  : expectedEventID - decodedEventID;
    if (distance > kMaxCounterJump) return false;

    reference.decodedEventID = decodedEventID;
    reference.triggerTimeNs = triggerTime * 1.0e9;
    reference.firstC3TransitionTimeNs = firstEventTransition * 1.0e9;
    reference.truthT0Ns = reference.triggerTimeNs - reference.firstC3TransitionTimeNs;
    return true;
}

void WriteWaveformDiagnostic(const WaveformSegment& eventCode, const WaveformSegment& trigger,
                             const ReferenceTime* reference,
                             const DecodeDiagnostics* diagnostics, int fileIndex, int segmentIndex,
                             double eventThreshold, double triggerThreshold, TDirectory& directory) {
    auto makeGraph = [](const WaveformSegment& waveform) {
        const size_t count = min(waveform.time.size(), waveform.amplitude.size());
        vector<double> timeNs(count);
        for (size_t i = 0; i < count; ++i) timeNs[i] = waveform.time[i] * 1.0e9;
        return TGraph(static_cast<int>(count), timeNs.data(), waveform.amplitude.data());
    };
    auto amplitudeRange = [](const WaveformSegment& waveform) {
        if (waveform.amplitude.empty()) return pair<double, double>{-1.0, 1.0};
        const auto [low, high] = minmax_element(waveform.amplitude.begin(), waveform.amplitude.end());
        const double margin = max(0.05, 0.05 * (*high - *low));
        return pair<double, double>{*low - margin, *high + margin};
    };
    auto drawThreshold = [](const WaveformSegment& waveform, double threshold, int color) {
        if (waveform.time.empty()) return;
        TLine line(waveform.time.front() * 1.0e9, threshold,
                   waveform.time.back() * 1.0e9, threshold);
        line.SetLineColor(color);
        line.SetLineStyle(2);
        line.DrawClone();
    };

    ostringstream baseName;
    baseName << "segment_" << setw(2) << setfill('0') << segmentIndex + 1;
    directory.cd();
    TCanvas canvas(baseName.str().c_str(), "C3/C4 waveform decoding", 1400, 1000);
    canvas.Divide(1, 2);

    TGraph eventGraph = makeGraph(eventCode);
    eventGraph.SetLineColor(kMagenta + 1);
    eventGraph.SetTitle("C3 event code;Time relative to oscilloscope trigger [ns];Amplitude [V]");
    canvas.cd(1);
    eventGraph.Draw("AL");
    drawThreshold(eventCode, eventThreshold, kRed + 1);
    if (diagnostics && !diagnostics->eventSampleTimesNs.empty()) {
        TGraph sampleGraph(static_cast<int>(diagnostics->eventSampleTimesNs.size()),
                           diagnostics->eventSampleTimesNs.data(),
                           diagnostics->eventSampleAmplitudes.data());
        sampleGraph.SetMarkerStyle(20);
        sampleGraph.SetMarkerSize(0.8);
        sampleGraph.SetMarkerColor(kBlack);
        sampleGraph.DrawClone("P SAME");
    }

    TGraph triggerGraph = makeGraph(trigger);
    triggerGraph.SetLineColor(kOrange + 7);
    triggerGraph.SetTitle("C4 trigger;Time relative to oscilloscope trigger [ns];Amplitude [V]");
    canvas.cd(2);
    triggerGraph.Draw("AL");
    drawThreshold(trigger, triggerThreshold, kRed + 1);
    if (reference) {
        const auto [low, high] = amplitudeRange(trigger);
        TLine crossing(reference->triggerTimeNs, low, reference->triggerTimeNs, high);
        crossing.SetLineColor(kGreen + 2);
        crossing.SetLineWidth(2);
        crossing.DrawClone();
    }

    canvas.cd(1);
    TLatex label;
    label.SetNDC();
    label.SetTextSize(0.045);
    ostringstream annotation;
    annotation << "Trace " << setw(5) << setfill('0') << fileIndex
               << ", segment " << segmentIndex + 1;
    if (reference) {
        annotation << ", eventID=" << reference->decodedEventID
                   << ", first C3 transition=" << fixed << setprecision(3)
                   << reference->firstC3TransitionTimeNs << " ns"
                   << ", T0=" << fixed << setprecision(3) << reference->truthT0Ns << " ns";
    } else {
        annotation << ", decode FAILED";
    }
    label.DrawLatex(0.10, 0.92, annotation.str().c_str());
    label.SetTextSize(0.035);
    if (diagnostics && !diagnostics->decodedBits.empty())
        label.DrawLatex(0.10, 0.84, ("Decoded bits: " + diagnostics->decodedBits).c_str());
    canvas.Write();
    string metadataText = annotation.str();
    if (diagnostics && !diagnostics->decodedBits.empty())
        metadataText += ", bits=" + diagnostics->decodedBits;
    TNamed metadata((baseName.str() + "_decoding").c_str(), metadataText.c_str());
    metadata.Write();
}

array<double, kEstimatorCount> ExtractTimes(const vector<int>& clusterIndices,
                                            const vector<Cluster>& clusters,
                                            const vector<StripHit>& stripHits,
                                            const vector<RawData>& detectorRawData,
                                            WaveformProcessor& timingFitter) {
    array<double, kEstimatorCount> result;
    result.fill(numeric_limits<double>::quiet_NaN());
    array<vector<double>, 2> planeStripTimes;
    map<int, double> fittedStripTimes;

    for (int clusterIndex : clusterIndices) {
        if (clusterIndex < 0 || clusterIndex >= static_cast<int>(clusters.size())) continue;
        const Cluster& cluster = clusters[clusterIndex];
        if (cluster.type < 0 || cluster.type > 1) continue;

        vector<double> stripTimes;
        for (int stripIndex : cluster.stripHitIndices) {
            if (stripIndex < 0 || stripIndex >= static_cast<int>(stripHits.size())) continue;
            const int rawIndex = stripHits[stripIndex].rawIndices;
            if (rawIndex < 0 || rawIndex >= static_cast<int>(detectorRawData.size())) continue;

            auto cached = fittedStripTimes.find(rawIndex);
            if (cached == fittedStripTimes.end()) {
                const StripHit fitted = timingFitter.ProcessWaveform(detectorRawData[rawIndex]);
                const double fittedTime = fitted.isValid && isfinite(fitted.time)
                                              ? fitted.time
                                              : numeric_limits<double>::quiet_NaN();
                cached = fittedStripTimes.emplace(rawIndex, fittedTime).first;
            }
            if (isfinite(cached->second)) stripTimes.push_back(cached->second);
        }
        if (stripTimes.empty()) continue;
        planeStripTimes[cluster.type].insert(
            planeStripTimes[cluster.type].end(), stripTimes.begin(), stripTimes.end());
    }

    vector<double> combined;
    for (size_t plane = 0; plane < 2; ++plane) {
        const auto& times = planeStripTimes[plane];
        if (times.empty()) continue;
        result[plane * 2] = accumulate(times.begin(), times.end(), 0.0) / times.size();
        result[plane * 2 + 1] = *min_element(times.begin(), times.end());
        combined.insert(combined.end(), times.begin(), times.end());
    }
    if (!combined.empty()) {
        result[4] = accumulate(combined.begin(), combined.end(), 0.0) / combined.size();
        result[5] = *min_element(combined.begin(), combined.end());
    }
    return result;
}

string SafeName(string value) {
    for (char& c : value)
        if (!isalnum(static_cast<unsigned char>(c))) c = '_';
    return value;
}

vector<double> SampleValues(const vector<TimeSample>& samples) {
    vector<double> values;
    values.reserve(samples.size());
    for (const auto& sample : samples) values.push_back(sample.value);
    return values;
}

FitResult WriteAndFit(const vector<double>& values, TDirectory& directory, const string& name,
                      const string& title, int bins) {
    FitResult result;
    result.entries = static_cast<long long>(values.size());
    if (values.size() < 3) return result;

    const double rawMean = accumulate(values.begin(), values.end(), 0.0) / values.size();
    double variance = 0.0;
    for (double value : values) variance += (value - rawMean) * (value - rawMean);
    const double rawSigma = sqrt(variance / (values.size() - 1));
    if (!isfinite(rawSigma) || rawSigma <= 0.0) return result;

    double low = rawMean - 6.0 * rawSigma;
    double high = rawMean + 6.0 * rawSigma;
    if (low == high) high = low + 1.0;
    directory.cd();
    TH1D histogram(name.c_str(), title.c_str(), bins, low, high);
    for (double value : values) histogram.Fill(value);

    TF1 gaussian((name + "_gaus").c_str(), "gaus", rawMean - 2.0 * rawSigma, rawMean + 2.0 * rawSigma);
    gaussian.SetParameters(histogram.GetMaximum(), rawMean, rawSigma);
    const int status = histogram.Fit(&gaussian, "Q");
    if (status == 0 && isfinite(gaussian.GetParameter(2))) {
        result.mean = gaussian.GetParameter(1);
        result.sigma = abs(gaussian.GetParameter(2));
        result.sigmaError = gaussian.GetParError(2);
    } else {
        result.mean = rawMean;
        result.sigma = rawSigma;
    }
    histogram.Write();
    return result;
}

}  // namespace

void TimeResolutionScript::LoadConfig(const json& config) {
    m_csvDirectory = config.value("csvDirectory", m_csvDirectory);
    m_trackFile = config.value("trackFile", m_trackFile);
    m_outputFile = config.value("outputFile", m_outputFile);
    m_timingCacheFile = config.value("timingCacheFile", m_timingCacheFile);
    m_waveformDiagnosticFile =
        config.value("waveformDiagnosticFile", m_waveformDiagnosticFile);
    m_rebuildTimingCache = config.value("rebuildTimingCache", m_rebuildTimingCache);
    m_writeWaveformDiagnostics =
        config.value("writeWaveformDiagnostics", m_writeWaveformDiagnostics);
    m_triggerThreshold = config.value("triggerThreshold", m_triggerThreshold);
    m_eventThreshold = config.value("eventThreshold", m_eventThreshold);
    m_eventBits = config.value("eventBits", m_eventBits);
    m_histogramBins = config.value("histogramBins", m_histogramBins);
    m_maxWaveformFiles = config.value("maxWaveformFiles", m_maxWaveformFiles);
}

void TimeResolutionScript::Print() const {
    cout << "TimeResolutionScript: CSV=" << m_csvDirectory
         << ", trigger threshold=" << m_triggerThreshold << " V, event bits=" << m_eventBits
         << ", rebuild cache=" << (m_rebuildTimingCache ? "yes" : "no") << '\n';
    if (m_writeWaveformDiagnostics)
        cout << "  waveform diagnostics=" << m_waveformDiagnosticFile << '\n';
}

bool TimeResolutionScript::Validate() const {
    return m_eventBits > 0 && m_eventBits <= 63 && m_histogramBins > 0;
}

bool TimeResolutionScript::Execute() {
    const string trackPath = m_trackFile.empty() ? GetOutputDir() + "TrackInfo.root" : m_trackFile;
    const string outputPath = filesystem::path(m_outputFile).is_absolute()
                                  ? m_outputFile
                                  : GetOutputDir() + m_outputFile;

    auto parser = GetParser();
    if (!parser) {
        cerr << "[TimeResolution] parser not set; raw waveforms are required for Fit timing\n";
        return false;
    }

    unique_ptr<TFile> trackFile(TFile::Open(trackPath.c_str(), "READ"));
    if (!trackFile || trackFile->IsZombie()) {
        cerr << "[TimeResolution] cannot open track file " << trackPath << '\n';
        return false;
    }
    auto* validation = dynamic_cast<TTree*>(trackFile->Get("TrackerValidation"));
    if (!validation) {
        cerr << "[TimeResolution] TrackerValidation is missing; rerun TrackAnalysis with saveValidationData=true\n";
        return false;
    }
    for (const char* branch : {"eventID", "rawEventID", "trackIndex", "detID",
                               "clusterIndices", "stripHits", "clusters"}) {
        if (!validation->GetBranch(branch)) {
            cerr << "[TimeResolution] branch " << branch << " is missing; rerun the updated TrackAnalysis\n";
            return false;
        }
    }

    Int_t eventID = 0;
    ULong64_t rawEventID = 0;
    Int_t trackIndex = 0, detID = 0;
    vector<Int_t>* clusterIndices = nullptr;
    vector<StripHit>* stripHits = nullptr;
    vector<Cluster>* clusters = nullptr;
    validation->SetBranchAddress("eventID", &eventID);
    validation->SetBranchAddress("rawEventID", &rawEventID);
    validation->SetBranchAddress("trackIndex", &trackIndex);
    validation->SetBranchAddress("detID", &detID);
    validation->SetBranchAddress("clusterIndices", &clusterIndices);
    validation->SetBranchAddress("stripHits", &stripHits);
    validation->SetBranchAddress("clusters", &clusters);

    map<TrackKey, map<int, DetectorTimes>> trackTimes;
    set<uint64_t> wantedEventIDs;
    WaveformProcessor timingFitter;
    timingFitter.LoadConfig(json{{"mode", "Fit"}});
    int loadedEventID = numeric_limits<int>::min();
    unordered_map<int, vector<RawData>> rawHits;
    size_t missingRawDetectors = 0;
    size_t emptyTimingDetectors = 0;
    for (Long64_t entry = 0; entry < validation->GetEntries(); ++entry) {
        validation->GetEntry(entry);
        if (!clusterIndices || !stripHits || !clusters) continue;
        if (eventID != loadedEventID) {
            rawHits = parser->LoadEvent(eventID);
            loadedEventID = eventID;
        }
        const auto rawDetector = rawHits.find(detID);
        if (rawDetector == rawHits.end()) {
            ++missingRawDetectors;
            continue;
        }
        DetectorTimes times;
        times.value = ExtractTimes(*clusterIndices, *clusters, *stripHits,
                                   rawDetector->second, timingFitter);
        bool hasAnyTime = false;
        for (double value : times.value) {
            if (isfinite(value)) {
                hasAnyTime = true;
                break;
            }
        }
        if (!hasAnyTime) ++emptyTimingDetectors;
        trackTimes[{rawEventID, trackIndex}][detID] = times;
        wantedEventIDs.insert(rawEventID);
    }
    cout << "[TimeResolution] loaded " << trackTimes.size() << " tracks and "
         << wantedEventIDs.size() << " raw event IDs\n";
    if (missingRawDetectors > 0 || emptyTimingDetectors > 0)
        cout << "[TimeResolution] Fit timing diagnostics: missing raw detector entries="
             << missingRawDetectors << ", empty detector times=" << emptyTimingDetectors << '\n';

    vector<pair<int, filesystem::path>> traceFiles;
    const regex filenamePattern(R"(^C3Trace([0-9]+)\.csv$)");
    const bool useExternalT0 = !m_csvDirectory.empty() && filesystem::is_directory(m_csvDirectory);
    if (!m_csvDirectory.empty() && !useExternalT0) {
        cerr << "[TimeResolution] CSV directory does not exist: " << m_csvDirectory
             << "; external T0 outputs will be skipped\n";
    }
    if (useExternalT0) {
        for (const auto& item : filesystem::directory_iterator(m_csvDirectory)) {
            smatch match;
            const string name = item.path().filename().string();
            if (item.is_regular_file() && regex_match(name, match, filenamePattern))
                traceFiles.emplace_back(stoi(match[1].str()), item.path());
        }
    }
    sort(traceFiles.begin(), traceFiles.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
    if (m_maxWaveformFiles > 0 && traceFiles.size() > static_cast<size_t>(m_maxWaveformFiles))
        traceFiles.resize(m_maxWaveformFiles);

    unique_ptr<TFile> waveformDiagnostics;
    if (m_writeWaveformDiagnostics && useExternalT0) {
        const string diagnosticPath = filesystem::path(m_waveformDiagnosticFile).is_absolute()
                                          ? m_waveformDiagnosticFile
                                          : GetOutputDir() + m_waveformDiagnosticFile;
        waveformDiagnostics.reset(TFile::Open(diagnosticPath.c_str(), "RECREATE"));
        if (!waveformDiagnostics || waveformDiagnostics->IsZombie()) {
            cerr << "[TimeResolution] cannot create waveform diagnostic file "
                 << diagnosticPath << '\n';
            return false;
        }
        cout << "[TimeResolution] writing every oscilloscope segment to " << diagnosticPath << '\n';
    } else if (m_writeWaveformDiagnostics) {
        cerr << "[TimeResolution] waveform diagnostics require a valid CSV directory; skipped\n";
    }

    map<uint64_t, ReferenceTime> references;
    vector<ReferenceTime> decodedReferences;
    size_t invalidSegments = 0, duplicateIDs = 0, processedFiles = 0;
    uint64_t expectedEventID = traceFiles.empty() ? 0 : static_cast<uint64_t>((traceFiles.front().first + 1) * 10);
    uint64_t firstDecodedEventID = 0, lastDecodedEventID = 0, missingEventCount = 0;
    size_t decodedSegments = 0, continuousSegments = 0, repeatedSegments = 0;
    size_t forwardGapSegments = 0, backwardSegments = 0;
    vector<string> sequenceAnomalies;
    const string cachePath = filesystem::path(m_timingCacheFile).is_absolute()
                                 ? m_timingCacheFile
                                 : GetOutputDir() + m_timingCacheFile;
    bool cacheLoaded = false;
    if (useExternalT0 && !m_writeWaveformDiagnostics && !m_rebuildTimingCache && filesystem::exists(cachePath)) {
        unique_ptr<TFile> cache(TFile::Open(cachePath.c_str(), "READ"));
        auto* metadata = cache ? dynamic_cast<TTree*>(cache->Get("CacheMetadata")) : nullptr;
        auto* cacheTree = cache ? dynamic_cast<TTree*>(cache->Get("TimingReferences")) : nullptr;
        Int_t cachedVersion = 0, cachedFiles = 0, cachedBits = 0;
        Double_t cachedTriggerThreshold = 0.0, cachedEventThreshold = 0.0;
        if (metadata && cacheTree && metadata->GetBranch("cacheVersion") &&
            cacheTree->GetBranch("firstC3TransitionTimeNs")) {
            metadata->SetBranchAddress("cacheVersion", &cachedVersion);
            metadata->SetBranchAddress("traceFileCount", &cachedFiles);
            metadata->SetBranchAddress("eventBits", &cachedBits);
            metadata->SetBranchAddress("triggerThreshold", &cachedTriggerThreshold);
            metadata->SetBranchAddress("eventThreshold", &cachedEventThreshold);
            metadata->GetEntry(0);
            cacheLoaded = cachedVersion == kTimingCacheVersion &&
                          cachedFiles == static_cast<Int_t>(traceFiles.size()) &&
                          cachedBits == m_eventBits &&
                          abs(cachedTriggerThreshold - m_triggerThreshold) < 1.0e-12 &&
                          abs(cachedEventThreshold - m_eventThreshold) < 1.0e-12;
        }
        if (cacheLoaded) {
            ULong64_t cachedEventID = 0;
            Double_t cachedTrigger = 0.0, cachedFirstBitEdge = 0.0, cachedT0 = 0.0;
            Int_t cachedFile = 0, cachedSegment = 0;
            cacheTree->SetBranchAddress("eventID", &cachedEventID);
            cacheTree->SetBranchAddress("triggerTimeNs", &cachedTrigger);
            cacheTree->SetBranchAddress("firstC3TransitionTimeNs", &cachedFirstBitEdge);
            cacheTree->SetBranchAddress("truthT0Ns", &cachedT0);
            cacheTree->SetBranchAddress("traceFileIndex", &cachedFile);
            cacheTree->SetBranchAddress("segmentIndex", &cachedSegment);
            for (Long64_t entry = 0; entry < cacheTree->GetEntries(); ++entry) {
                cacheTree->GetEntry(entry);
                if (!wantedEventIDs.count(cachedEventID)) continue;
                references.emplace(cachedEventID, ReferenceTime{
                    cachedEventID, cachedTrigger, cachedFirstBitEdge, cachedT0, cachedFile,
                    cachedSegment});
            }
            cout << "[TimeResolution] timing cache=" << cachePath
                 << ", matched IDs=" << references.size() << '\n';
        }
    }

    if (useExternalT0 && !cacheLoaded) for (const auto& [fileIndex, eventPath] : traceFiles) {
        ostringstream number;
        number << setw(5) << setfill('0') << fileIndex;
        const filesystem::path triggerPath = filesystem::path(m_csvDirectory) / ("C4Trace" + number.str() + ".csv");
        if (!filesystem::exists(triggerPath)) {
            cerr << "[TimeResolution] missing C4 partner for trace " << number.str() << '\n';
            continue;
        }

        TraceFile eventCode, trigger;
        string eventError, triggerError;
        auto eventRead = async(launch::async, [&] {
            return ReadTraceFile(eventPath, eventCode, eventError, -20.0e-9, 550.0e-9);
        });
        auto triggerRead = async(launch::async, [&] {
            return ReadTraceFile(triggerPath, trigger, triggerError, -20.0e-9, 30.0e-9);
        });
        const bool eventOK = eventRead.get();
        const bool triggerOK = triggerRead.get();
        if (!eventOK || !triggerOK) {
            const string& error = !eventOK ? eventError : triggerError;
            cerr << "[TimeResolution] " << error << '\n';
            continue;
        }
        const int segments = min(eventCode.segmentCount, trigger.segmentCount);
        TDirectory* traceDiagnosticDirectory = nullptr;
        if (waveformDiagnostics) {
            ostringstream traceName;
            traceName << "Trace" << setw(5) << setfill('0') << fileIndex;
            traceDiagnosticDirectory = waveformDiagnostics->mkdir(traceName.str().c_str());
        }
        for (int segment = 0; segment < segments; ++segment) {
            ReferenceTime reference;
            DecodeDiagnostics diagnostics;
            if (!DecodeReference(eventCode.segments[segment], trigger.segments[segment],
                                 m_triggerThreshold, m_eventThreshold, m_eventBits,
                                 expectedEventID, reference,
                                 traceDiagnosticDirectory ? &diagnostics : nullptr)) {
                if (traceDiagnosticDirectory)
                    WriteWaveformDiagnostic(
                        eventCode.segments[segment], trigger.segments[segment], nullptr, nullptr,
                        fileIndex, segment, m_eventThreshold, m_triggerThreshold,
                        *traceDiagnosticDirectory);
                ++invalidSegments;
                ++expectedEventID;
                continue;
            }
            if (traceDiagnosticDirectory)
                WriteWaveformDiagnostic(
                    eventCode.segments[segment], trigger.segments[segment], &reference,
                    &diagnostics, fileIndex, segment, m_eventThreshold, m_triggerThreshold,
                    *traceDiagnosticDirectory);
            if (decodedSegments == 0) firstDecodedEventID = reference.decodedEventID;
            const int64_t delta = static_cast<int64_t>(reference.decodedEventID) -
                                  static_cast<int64_t>(expectedEventID);
            if (delta == 0) {
                ++continuousSegments;
            } else {
                ostringstream anomaly;
                anomaly << "trace=" << setw(5) << setfill('0') << fileIndex
                        << " segment=" << segment + 1 << " expected=" << expectedEventID
                        << " decoded=" << reference.decodedEventID << " delta=" << delta;
                if (sequenceAnomalies.size() < 100) sequenceAnomalies.push_back(anomaly.str());
                if (delta == -1) {
                    ++repeatedSegments;
                } else if (delta > 0) {
                    ++forwardGapSegments;
                    missingEventCount += static_cast<uint64_t>(delta);
                } else {
                    ++backwardSegments;
                }
            }
            ++decodedSegments;
            lastDecodedEventID = reference.decodedEventID;
            expectedEventID = reference.decodedEventID + 1;
            reference.fileIndex = fileIndex;
            reference.segmentIndex = segment;
            decodedReferences.push_back(reference);
            if (!wantedEventIDs.count(reference.decodedEventID)) continue;
            if (!references.emplace(reference.decodedEventID, reference).second) ++duplicateIDs;
        }
        ++processedFiles;
        if (processedFiles % 25 == 0 || processedFiles == traceFiles.size())
            cout << "\r[TimeResolution] waveform files " << processedFiles << '/' << traceFiles.size()
                 << ", matched IDs=" << references.size() << flush;
        if (!m_writeWaveformDiagnostics && references.size() == wantedEventIDs.size()) break;
    }
    if (waveformDiagnostics) waveformDiagnostics->Flush();
    if (useExternalT0 && !cacheLoaded) cout << "\n[TimeResolution] invalid segments=" << invalidSegments
         << ", duplicate matched IDs=" << duplicateIDs << '\n';
    if (useExternalT0 && !cacheLoaded) cout << "[TimeResolution] C3 sequence: decoded=" << decodedSegments
         << "/" << traceFiles.size() * 10
         << ", first=" << firstDecodedEventID << ", last=" << lastDecodedEventID
         << ", continuous(+1)=" << continuousSegments
         << ", repeated=" << repeatedSegments
         << ", forward gaps=" << forwardGapSegments
         << ", missing IDs=" << missingEventCount
         << ", backward=" << backwardSegments << '\n';
    if (useExternalT0 && !cacheLoaded) for (const string& anomaly : sequenceAnomalies)
        cout << "[TimeResolution] C3 anomaly: " << anomaly << '\n';

    if (useExternalT0 && !cacheLoaded && processedFiles == traceFiles.size()) {
        unique_ptr<TFile> cache(TFile::Open(cachePath.c_str(), "RECREATE"));
        if (cache && !cache->IsZombie()) {
            TTree cacheTree("TimingReferences", "Decoded oscilloscope timing references");
            ULong64_t cachedEventID = 0;
            Double_t cachedTrigger = 0.0, cachedFirstBitEdge = 0.0, cachedT0 = 0.0;
            Int_t cachedFile = 0, cachedSegment = 0;
            cacheTree.Branch("eventID", &cachedEventID);
            cacheTree.Branch("triggerTimeNs", &cachedTrigger);
            cacheTree.Branch("firstC3TransitionTimeNs", &cachedFirstBitEdge);
            cacheTree.Branch("truthT0Ns", &cachedT0);
            cacheTree.Branch("traceFileIndex", &cachedFile);
            cacheTree.Branch("segmentIndex", &cachedSegment);
            for (const ReferenceTime& reference : decodedReferences) {
                cachedEventID = reference.decodedEventID;
                cachedTrigger = reference.triggerTimeNs;
                cachedFirstBitEdge = reference.firstC3TransitionTimeNs;
                cachedT0 = reference.truthT0Ns;
                cachedFile = reference.fileIndex;
                cachedSegment = reference.segmentIndex;
                cacheTree.Fill();
            }
            cacheTree.Write();

            TTree metadata("CacheMetadata", "Oscilloscope timing cache configuration");
            Int_t cachedVersion = kTimingCacheVersion;
            Int_t cachedFiles = static_cast<Int_t>(traceFiles.size());
            Int_t cachedBits = m_eventBits;
            Double_t cachedTriggerThreshold = m_triggerThreshold;
            Double_t cachedEventThreshold = m_eventThreshold;
            metadata.Branch("cacheVersion", &cachedVersion);
            metadata.Branch("traceFileCount", &cachedFiles);
            metadata.Branch("eventBits", &cachedBits);
            metadata.Branch("triggerThreshold", &cachedTriggerThreshold);
            metadata.Branch("eventThreshold", &cachedEventThreshold);
            metadata.Fill();
            metadata.Write();
            cout << "[TimeResolution] wrote timing cache=" << cachePath
                 << " entries=" << decodedReferences.size() << '\n';
        }
    }
    if (useExternalT0 && references.empty())
        cerr << "[TimeResolution] no decoded C3 event ID matched rawEventID; "
             << "external T0 outputs will be skipped\n";

    unique_ptr<TFile> output(TFile::Open(outputPath.c_str(), "RECREATE"));
    if (!output || output->IsZombie()) {
        cerr << "[TimeResolution] cannot create " << outputPath << '\n';
        return false;
    }

    TTree eventTree("EventTimes", "Matched tracker and oscilloscope times");
    ULong64_t outRawEventID = 0;
    Int_t outTrackIndex = 0, outDetID = 0, outFileIndex = 0, outSegmentIndex = 0;
    Double_t triggerTimeNs = 0.0;
    Double_t firstC3TransitionTimeNs = 0.0, truthT0Ns = 0.0;
    array<Double_t, kEstimatorCount> estimatorValues{};
    eventTree.Branch("rawEventID", &outRawEventID);
    eventTree.Branch("trackIndex", &outTrackIndex);
    eventTree.Branch("detectorID", &outDetID);
    eventTree.Branch("traceFileIndex", &outFileIndex);
    eventTree.Branch("segmentIndex", &outSegmentIndex);
    eventTree.Branch("triggerTimeNs", &triggerTimeNs);
    eventTree.Branch("firstC3TransitionTimeNs", &firstC3TransitionTimeNs);
    eventTree.Branch("truthT0Ns", &truthT0Ns);
    for (size_t estimator = 0; estimator < kEstimatorCount; ++estimator)
        eventTree.Branch(kEstimatorNames[estimator].c_str(), &estimatorValues[estimator]);

    map<pair<int, size_t>, vector<TimeSample>> detectorTimeSamples;
    map<tuple<int, int, size_t>, vector<TimeSample>> pairResiduals;
    map<pair<int, size_t>, vector<TimeSample>> truthResiduals;
    size_t tracksWithReference = 0;
    for (const auto& [key, detectors] : trackTimes) {
        const auto referenceIt = references.find(key.first);
        const bool hasReference = referenceIt != references.end();
        const ReferenceTime* reference = hasReference ? &referenceIt->second : nullptr;
        if (hasReference) ++tracksWithReference;
        for (const auto& [detectorID, times] : detectors) {
            outRawEventID = key.first;
            outTrackIndex = key.second;
            outDetID = detectorID;
            outFileIndex = hasReference ? reference->fileIndex : -1;
            outSegmentIndex = hasReference ? reference->segmentIndex : -1;
            triggerTimeNs = hasReference ? reference->triggerTimeNs : numeric_limits<double>::quiet_NaN();
            firstC3TransitionTimeNs = hasReference
                                          ? reference->firstC3TransitionTimeNs
                                          : numeric_limits<double>::quiet_NaN();
            truthT0Ns = hasReference ? reference->truthT0Ns : numeric_limits<double>::quiet_NaN();
            estimatorValues = times.value;
            eventTree.Fill();
            for (size_t estimator = 0; estimator < kEstimatorCount; ++estimator) {
                if (!isfinite(times.value[estimator])) continue;
                detectorTimeSamples[{detectorID, estimator}].push_back(
                    {times.value[estimator], truthT0Ns});
                if (hasReference) {
                    truthResiduals[{detectorID, estimator}].push_back(
                        {times.value[estimator] - reference->truthT0Ns, reference->truthT0Ns});
                }
            }
        }

        for (auto first = detectors.begin(); first != detectors.end(); ++first) {
            for (auto second = next(first); second != detectors.end(); ++second) {
                for (size_t estimator = 0; estimator < kEstimatorCount; ++estimator) {
                    const double a = first->second.value[estimator];
                    const double b = second->second.value[estimator];
                    if (isfinite(a) && isfinite(b))
                        pairResiduals[{first->first, second->first, estimator}].push_back(
                            {a - b});
                }
            }
        }
    }
    eventTree.Write();

    TTree summary("ResolutionSummary", "Gaussian timing-resolution fits");
    string analysis, estimatorName;
    Int_t detectorA = 0, detectorB = -1;
    Long64_t entries = 0;
    Double_t meanNs = 0.0, sigmaNs = 0.0, sigmaErrorNs = 0.0, resolutionNs = 0.0;
    summary.Branch("analysis", &analysis);
    summary.Branch("estimator", &estimatorName);
    summary.Branch("detectorA", &detectorA);
    summary.Branch("detectorB", &detectorB);
    summary.Branch("entries", &entries);
    summary.Branch("meanNs", &meanNs);
    summary.Branch("sigmaNs", &sigmaNs);
    summary.Branch("sigmaErrorNs", &sigmaErrorNs);
    summary.Branch("resolutionNs", &resolutionNs);

    TDirectory* timeDirectory = output->mkdir("TimeDistributions");
    TDirectory* rawTimeDirectory = timeDirectory->mkdir("RawDetectorTime");
    const bool haveExternalT0 = tracksWithReference > 0;
    for (const auto& [key, samples] : detectorTimeSamples) {
        detectorA = key.first;
        detectorB = -1;
        const size_t estimator = key.second;
        estimatorName = kEstimatorNames[estimator];
        entries = static_cast<Long64_t>(samples.size());
        const string name = SafeName("time_d" + to_string(detectorA) + "_" + estimatorName);
        analysis = "detector_time";
        const FitResult rawTimeFit = WriteAndFit(
            SampleValues(samples), *rawTimeDirectory, name,
            "Detector time; t_{detector} [ns];Entries", m_histogramBins);
        entries = rawTimeFit.entries;
        meanNs = rawTimeFit.mean;
        sigmaNs = rawTimeFit.sigma;
        sigmaErrorNs = rawTimeFit.sigmaError;
        resolutionNs = rawTimeFit.sigma;
        summary.Fill();
    }

    TDirectory* pairDirectory = output->mkdir("PairResiduals");
    TDirectory* pairRawDirectory = pairDirectory->mkdir("Raw");
    for (const auto& [key, samples] : pairResiduals) {
        size_t estimator = 0;
        tie(detectorA, detectorB, estimator) = key;
        estimatorName = kEstimatorNames[estimator];
        const string name = SafeName("pair_d" + to_string(detectorA) + "_d" + to_string(detectorB) + "_" + estimatorName);
        const string title = "Tracker pair residual; t_{" + to_string(detectorA) +
                             "} - t_{" + to_string(detectorB) + "} [ns];Entries";
        const vector<double> rawValues = SampleValues(samples);
        analysis = "tracker_pair";
        const FitResult fit = WriteAndFit(rawValues, *pairRawDirectory, name, title,
                                         m_histogramBins);
        entries = fit.entries;
        meanNs = fit.mean;
        sigmaNs = fit.sigma;
        sigmaErrorNs = fit.sigmaError;
        resolutionNs = fit.sigma / sqrt(2.0);
        summary.Fill();
    }

    if (haveExternalT0) {
        TDirectory* truthDirectory = output->mkdir("TruthResiduals");
        TDirectory* truthRawDirectory = truthDirectory->mkdir("Raw");
        for (const auto& [key, samples] : truthResiduals) {
            detectorA = key.first;
            detectorB = -1;
            const size_t estimator = key.second;
            estimatorName = kEstimatorNames[estimator];
            const string name = SafeName("truth_d" + to_string(detectorA) + "_" + estimatorName);
            const vector<double> rawValues = SampleValues(samples);
            analysis = "truth_t0";
            const FitResult fit = WriteAndFit(rawValues, *truthRawDirectory, name,
                "External-T0 corrected time; t_{detector} - T0 [ns];Entries", m_histogramBins);
            entries = fit.entries;
            meanNs = fit.mean;
            sigmaNs = fit.sigma;
            sigmaErrorNs = fit.sigmaError;
            resolutionNs = fit.sigma;
            summary.Fill();
        }
    }
    output->cd();
    summary.Write();
    cout << "[TimeResolution] tracks=" << trackTimes.size()
         << ", with external T0=" << tracksWithReference
         << ", output=" << outputPath << '\n';
    return !trackTimes.empty();
}

REGISTER_SCRIPT("TimeResolution", TimeResolutionScript);
