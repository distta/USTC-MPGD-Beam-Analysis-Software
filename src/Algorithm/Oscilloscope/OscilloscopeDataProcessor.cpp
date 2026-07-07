#include "Algorithm/algorithms/OscilloscopeDataProcessor.h"

#include "Algorithm/AlgorithmFactory.h"

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraph.h>
#include <TLatex.h>
#include <TLine.h>
#include <TNamed.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <regex>
#include <sstream>
#include <utility>
#include <vector>

using namespace std;

namespace {

constexpr int kTimingCacheVersion = 6;
constexpr double kEventBitPeriodSeconds = 25.0e-9;

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

bool ParseTraceHeader(ifstream& input, const filesystem::path& path,
                      int& segmentCount, int& segmentSize, string& error) {
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
        segmentCount = stoi(fields[1]);
        segmentSize = stoi(fields[3]);
    } catch (...) {
        error = "invalid segment counts in " + path.string();
        return false;
    }
    if (segmentCount <= 0 || segmentSize <= 1) {
        error = "non-positive segment dimensions in " + path.string();
        return false;
    }

    while (getline(input, line)) {
        if (line.rfind("Time,Ampl", 0) == 0) return true;
    }
    error = "Time,Ampl header missing in " + path.string();
    return false;
}

}  // namespace

struct OscilloscopeDataProcessor::StreamTraceReader {
    ifstream input;
    filesystem::path path;
    int segmentCount = 0;
    int segmentSize = 0;
    double timeMin = 0.0;
    double timeMax = 0.0;

    bool Open(const filesystem::path& inputPath, double minTime,
              double maxTime, string& error) {
        path = inputPath;
        timeMin = minTime;
        timeMax = maxTime;
        input.open(path);
        if (!input) {
            error = "cannot open " + path.string();
            return false;
        }
        return ParseTraceHeader(input, path, segmentCount, segmentSize, error);
    }

    bool ReadSegment(int segmentIndex, WaveformSegment& segment) {
        segment.time.clear();
        segment.amplitude.clear();
        segment.time.reserve(static_cast<size_t>(segmentSize / 2));
        segment.amplitude.reserve(static_cast<size_t>(segmentSize / 2));

        string line;
        double firstTime = 0.0;
        double pitch = 0.0;
        for (int localSample = 0; localSample < segmentSize; ++localSample) {
            if (!getline(input, line)) {
                return false;
            }

            double time = 0.0;
            double amplitude = 0.0;
            if (localSample < 2) {
                if (!ParsePair(line, time, amplitude)) continue;
                if (localSample == 0) {
                    firstTime = time;
                } else {
                    pitch = time - firstTime;
                }
            } else {
                time = firstTime + localSample * pitch;
                const char* comma = strchr(line.c_str(), ',');
                if (!comma) continue;
                char* end = nullptr;
                amplitude = strtod(comma + 1, &end);
                if (end == comma + 1) continue;
            }

            if (time >= timeMin && time <= timeMax) {
                segment.time.push_back(time);
                segment.amplitude.push_back(amplitude);
            }
        }
        return true;
    }
};

double OscilloscopeDataProcessor::DecodeArriveTime(const WaveformSegment& waveform,
                                                   double threshold, double startTime, bool risingFlag) const {
    const size_t count = min(waveform.time.size(), waveform.amplitude.size());
    if (count < 2) return numeric_limits<double>::quiet_NaN();

    auto timeIt = lower_bound(waveform.time.begin(), waveform.time.begin() + count, startTime);
    size_t begin = timeIt == waveform.time.begin()
                       ? 1
                       : static_cast<size_t>(timeIt - waveform.time.begin());
    for (size_t i = begin; i < count; ++i) {
        const double y0 = waveform.amplitude[i - 1];
        const double y1 = waveform.amplitude[i];

        const bool rising = y0 < threshold && y1 >= threshold;
        const bool falling = y0 >= threshold && y1 < threshold;

        if (risingFlag) {
            if (!rising) continue;
        } else {
            if (!rising && !falling) continue;
        }

        const double fraction = (threshold - y0) / (y1 - y0);
        const double crossing = waveform.time[i - 1] +
                                fraction * (waveform.time[i] - waveform.time[i - 1]);
        if (crossing >= startTime) return crossing;
    }
    return numeric_limits<double>::quiet_NaN();
}

uint64_t OscilloscopeDataProcessor::DecodeEventID(const WaveformSegment& eventCode, double startTime) const {
    uint64_t decodedEventID = 0;
    const size_t count = min(eventCode.time.size(), eventCode.amplitude.size());
    if (count == 0) return numeric_limits<uint64_t>::max();

    for (int bit = 0; bit < m_eventBits; ++bit) {
        const double sampleTime = startTime + (bit + 0.5) * kEventBitPeriodSeconds;
        auto sampleIt = lower_bound(eventCode.time.begin(), eventCode.time.begin() + count,
                                    sampleTime);
        if (sampleIt == eventCode.time.begin() + count || sampleTime < eventCode.time.front())
            return numeric_limits<uint64_t>::max();

        size_t sample = static_cast<size_t>(sampleIt - eventCode.time.begin());
        if (sample > 0 && sampleTime - eventCode.time[sample - 1] <= eventCode.time[sample] - sampleTime)
            --sample;

        const double amplitude = eventCode.amplitude[sample];
        const bool value = amplitude >= m_eventThreshold;
        decodedEventID = (decodedEventID << 1U) | static_cast<uint64_t>(value);
    }
    return decodedEventID;
}

bool OscilloscopeDataProcessor::DecodeOscilloscopeData(
    const WaveformSegment& eventCode, const WaveformSegment& trigger,
    OscilloscopeData& oscilloscopeData) const {
    if (eventCode.time.empty() || trigger.time.empty()) return false;

    const double triggerTime =
        DecodeArriveTime(trigger, m_triggerThreshold, trigger.time.front(), true);

    if (!isfinite(triggerTime)) return false;

    const double firstEventTransition =
        DecodeArriveTime(eventCode, m_eventThreshold, triggerTime, false);
    if (!isfinite(firstEventTransition)) return false;

    const uint64_t decodedEventID =
        DecodeEventID(eventCode, firstEventTransition);
    if (decodedEventID == numeric_limits<uint64_t>::max()) return false;

    oscilloscopeData.decodedEventID = decodedEventID;
    oscilloscopeData.triggerTime = triggerTime * 1.0e9;
    oscilloscopeData.EvnentIDtime = firstEventTransition * 1.0e9;
    oscilloscopeData.truthT0 = oscilloscopeData.triggerTime - oscilloscopeData.EvnentIDtime;
    return true;
}

void OscilloscopeDataProcessor::WriteWaveformDiagnostic(
    const WaveformSegment& eventCode,
    const WaveformSegment& trigger,
    const OscilloscopeData* oscilloscopeData, int fileIndex, int segmentIndex,
    TDirectory& directory) const {
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
    drawThreshold(eventCode, m_eventThreshold, kRed + 1);

    TGraph triggerGraph = makeGraph(trigger);
    triggerGraph.SetLineColor(kOrange + 7);
    triggerGraph.SetTitle("C4 trigger;Time relative to oscilloscope trigger [ns];Amplitude [V]");
    canvas.cd(2);
    triggerGraph.Draw("AL");
    drawThreshold(trigger, m_triggerThreshold, kRed + 1);
    if (oscilloscopeData) {
        const auto [low, high] = amplitudeRange(trigger);
        TLine crossing(oscilloscopeData->triggerTime, low, oscilloscopeData->triggerTime, high);
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
    if (oscilloscopeData) {
        annotation << ", eventID=" << oscilloscopeData->decodedEventID
                   << ", EventID Arrive Time=" << fixed << setprecision(3)
                   << oscilloscopeData->EvnentIDtime << " ns"
                   << ", T0=" << fixed << setprecision(3) << oscilloscopeData->truthT0
                   << " ns";
    } else {
        annotation << ", decode FAILED";
    }
    label.DrawLatex(0.10, 0.92, annotation.str().c_str());
    canvas.Write();
    string metadataText = annotation.str();
    TNamed metadata((baseName.str() + "_decoding").c_str(), metadataText.c_str());
    metadata.Write();
}

bool OscilloscopeDataProcessor::LoadOscilloscopeDataCache(const string& cachePath) {

    unique_ptr<TFile> cache(TFile::Open(cachePath.c_str(), "READ"));
    auto* cacheTree = cache ? dynamic_cast<TTree*>(cache->Get("OscilloscopeData")) : nullptr;
    if (!cacheTree) {
        return false;
    }

    ULong64_t eventID = 0;
    Double_t triggerTime = 0.0, EventIDTime = 0.0, t0 = 0.0;
    Int_t file = 0, segment = 0;
    cacheTree->SetBranchAddress("eventID", &eventID);
    cacheTree->SetBranchAddress("triggerTime", &triggerTime);
    cacheTree->SetBranchAddress("EventIDTime", &EventIDTime);
    cacheTree->SetBranchAddress("truthT0", &t0);
    cacheTree->SetBranchAddress("traceFileIndex", &file);
    cacheTree->SetBranchAddress("segmentIndex", &segment);
    for (Long64_t entry = 0; entry < cacheTree->GetEntries(); ++entry) {
        cacheTree->GetEntry(entry);
        m_dataByEventID.emplace(eventID, OscilloscopeData{eventID, triggerTime, EventIDTime, t0, file, segment});
    }
    cout << "[OscilloscopeData] data cache=" << cachePath << ", entries=" << m_dataByEventID.size() << '\n';
    return true;
}

void OscilloscopeDataProcessor::WriteOscilloscopeDataCache(
    const string& cachePath, const vector<OscilloscopeData>& decodedData) const {
    unique_ptr<TFile> cache(TFile::Open(cachePath.c_str(), "RECREATE"));
    if (!cache || cache->IsZombie()) return;

    TTree cacheTree("OscilloscopeData", "Decoded oscilloscope data");
    ULong64_t eventID = 0;
    Double_t triggerTime = 0.0, EventIDTime = 0.0, t0 = 0.0;
    Int_t file = 0, segment = 0;
    cacheTree.Branch("eventID", &eventID);
    cacheTree.Branch("triggerTime", &triggerTime);
    cacheTree.Branch("EventIDTime", &EventIDTime);
    cacheTree.Branch("truthT0", &t0);
    cacheTree.Branch("traceFileIndex", &file);
    cacheTree.Branch("segmentIndex", &segment);
    for (const OscilloscopeData& oscilloscopeData : decodedData) {
        eventID = oscilloscopeData.decodedEventID;
        triggerTime = oscilloscopeData.triggerTime;
        EventIDTime = oscilloscopeData.EvnentIDtime;
        t0 = oscilloscopeData.truthT0;
        file = oscilloscopeData.fileIndex;
        segment = oscilloscopeData.segmentIndex;
        cacheTree.Fill();
    }
    cacheTree.Write();

    cout << "[OscilloscopeData] wrote data cache=" << cachePath << " entries=" << decodedData.size() << '\n';
}

void OscilloscopeDataProcessor::LoadConfig(const json& config) {
    m_csvDirectory = config.value("csvDirectory", m_csvDirectory);
    m_outputDir = config.value("outputDir", m_outputDir);
    m_dataCacheFile = config.value("dataCacheFile",
                                   config.value("timingCacheFile", m_dataCacheFile));
    m_waveformDiagnosticFile =
        config.value("waveformDiagnosticFile", m_waveformDiagnosticFile);
    m_rebuildTimingCache = config.value("rebuildTimingCache", m_rebuildTimingCache);
    m_writeWaveformDiagnostics =
        config.value("writeWaveformDiagnostics", m_writeWaveformDiagnostics);
    m_triggerThreshold = config.value("triggerThreshold", m_triggerThreshold);
    m_eventThreshold = config.value("eventThreshold", m_eventThreshold);
    m_eventBits = config.value("eventBits", m_eventBits);
    m_maxWaveformFiles = config.value("maxWaveformFiles", m_maxWaveformFiles);
}

void OscilloscopeDataProcessor::Print() const {
    cout << "OscilloscopeDataProcessor: CSV=" << m_csvDirectory
         << ", trigger threshold=" << m_triggerThreshold
         << " V, event bits=" << m_eventBits
         << ", event channel=" << m_eventCodeChannel
         << ", trigger channel=" << m_triggerChannel
         << ", rebuild cache=" << (m_rebuildTimingCache ? "yes" : "no") << '\n';
    if (m_writeWaveformDiagnostics)
        cout << "  waveform diagnostics=" << m_waveformDiagnosticFile << '\n';
}

bool OscilloscopeDataProcessor::Initialize() {
    m_initialized = false;
    m_dataByEventID.clear();
    m_decodedEntries = 0;

    const string cachePath = filesystem::path(m_dataCacheFile).is_absolute()
                                 ? m_dataCacheFile
                                 : m_outputDir + m_dataCacheFile;
    const bool cacheLoaded = !m_rebuildTimingCache &&
                             filesystem::exists(cachePath) &&
                             LoadOscilloscopeDataCache(cachePath);

    if (cacheLoaded) {
        m_initialized = true;
        return true;
    }

    if (m_csvDirectory.empty() || !filesystem::is_directory(m_csvDirectory)) {
        cerr << "[TimeResolution] CSV directory does not exist: " << m_csvDirectory << endl;
        return false;
    }

    map<int, map<int, filesystem::path>> traceFiles;
    const regex filenamePattern(R"(^(C[0-9]+)Trace([0-9]+)\.csv$)");
    for (const auto& item : filesystem::directory_iterator(m_csvDirectory)) {
        smatch match;
        const string name = item.path().filename().string();
        if (item.is_regular_file() && regex_match(name, match, filenamePattern))
            traceFiles[stoi(match[2].str())][stoi(match[1].str().substr(1))] = item.path();
    }

    vector<int> traceIndices;
    traceIndices.reserve(traceFiles.size());
    for (const auto& [traceIndex, channels] : traceFiles)
        traceIndices.push_back(traceIndex);
    if (m_maxWaveformFiles > 0 &&
        traceIndices.size() > static_cast<size_t>(m_maxWaveformFiles))
        traceIndices.resize(m_maxWaveformFiles);

    if (traceIndices.empty()) {
        cerr << "[OscilloscopeData] no C?Trace*.csv files found in " << m_csvDirectory << '\n';
        return false;
    }

    unique_ptr<TFile> waveformDiagnostics;
    if (m_writeWaveformDiagnostics) {
        const string diagnosticPath = filesystem::path(m_waveformDiagnosticFile).is_absolute()
                                          ? m_waveformDiagnosticFile
                                          : m_outputDir + m_waveformDiagnosticFile;
        waveformDiagnostics.reset(TFile::Open(diagnosticPath.c_str(), "RECREATE"));
        if (!waveformDiagnostics || waveformDiagnostics->IsZombie()) {
            cerr << "[OscilloscopeData] cannot create waveform diagnostic file " << diagnosticPath << '\n';
            return false;
        }
        cout << "[OscilloscopeData] writing every oscilloscope segment to " << diagnosticPath << '\n';
    }

    vector<OscilloscopeData> decodedData;
    size_t processedFiles = 0;
    size_t failedFiles = 0;
    size_t invalidSegments = 0;

    for (int fileIndex : traceIndices) {
        
        const auto& channelFiles = traceFiles[fileIndex];

        const auto eventFile = channelFiles.find(m_eventCodeChannel);
        const auto triggerFile = channelFiles.find(m_triggerChannel);
        if (eventFile == channelFiles.end() || triggerFile == channelFiles.end()) {
            cerr << "[OscilloscopeData] missing C" << m_eventCodeChannel
                 << " or C" << m_triggerChannel << " for Trace" << fileIndex << '\n';
            ++failedFiles;
            ++processedFiles;
            continue;
        }

        StreamTraceReader eventReader;
        StreamTraceReader triggerReader;
        string error;
        if (!eventReader.Open(eventFile->second, -20.0e-9, 550.0e-9, error) ||
            !triggerReader.Open(triggerFile->second, -20.0e-9, 30.0e-9, error)) {
            cerr << "[OscilloscopeData] " << error << '\n';
            ++failedFiles;
            ++processedFiles;
            continue;
        }

        const int segments = min(eventReader.segmentCount, triggerReader.segmentCount);
        TDirectory* traceDiagnosticDirectory = nullptr;
        if (waveformDiagnostics) {
            ostringstream traceName;
            traceName << "Trace" << setw(5) << setfill('0') << fileIndex;
            traceDiagnosticDirectory = waveformDiagnostics->mkdir(traceName.str().c_str());
        }
        WaveformSegment eventCode;
        WaveformSegment trigger;
        for (int segment = 0; segment < segments; ++segment) {
            if (!eventReader.ReadSegment(segment, eventCode) ||
                !triggerReader.ReadSegment(segment, trigger)) {
                ++failedFiles;
                break;
            }

            OscilloscopeData oscilloscopeData;
            if (!DecodeOscilloscopeData(eventCode, trigger, oscilloscopeData)) {
                if (traceDiagnosticDirectory)
                    WriteWaveformDiagnostic(eventCode, trigger, nullptr, fileIndex, segment,
                                            *traceDiagnosticDirectory);
                ++invalidSegments;
                continue;
            }
            if (traceDiagnosticDirectory)
                WriteWaveformDiagnostic(eventCode, trigger, &oscilloscopeData, fileIndex, segment,
                                        *traceDiagnosticDirectory);
            oscilloscopeData.fileIndex = fileIndex;
            oscilloscopeData.segmentIndex = segment;
            decodedData.push_back(oscilloscopeData);
            m_dataByEventID.emplace(oscilloscopeData.decodedEventID, oscilloscopeData);
        }
        ++processedFiles;
        if (processedFiles % 25 == 0 || processedFiles == traceIndices.size())
            cout << "\r[OscilloscopeData] waveform files " << processedFiles
                 << '/' << traceIndices.size()
                 << ", decoded IDs=" << m_dataByEventID.size() << flush;
    }
    if (!traceIndices.empty()) cout << '\n';
    cout << "[OscilloscopeData] failed files=" << failedFiles
         << ", invalid segments=" << invalidSegments
         << endl;

    if (waveformDiagnostics) waveformDiagnostics->Flush();

    if (!cacheLoaded && failedFiles == 0)
        WriteOscilloscopeDataCache(cachePath, decodedData);
    m_decodedEntries = cacheLoaded ? m_dataByEventID.size() : decodedData.size();
    m_initialized = true;
    return true;
}

OscilloscopeDataResult OscilloscopeDataProcessor::LoadData(
    const set<uint64_t>& wantedEventIDs) const {
    OscilloscopeDataResult result;
    if (!m_initialized) {
        cerr << "[OscilloscopeData] LoadData called before Initialize; no data loaded\n";
        return result;
    }

    for (uint64_t eventID : wantedEventIDs) {
        auto found = m_dataByEventID.find(eventID);
        if (found != m_dataByEventID.end())
            result.dataByEventID.emplace(eventID, found->second);
    }
    if (result.dataByEventID.empty())
        cerr << "[TimeResolution] no decoded C3 event ID matched rawEventID; "
             << "external T0 outputs will be skipped\n";
    result.decodedEntries = m_decodedEntries;
    return result;
}

REGISTER_ALGORITHM("OscilloscopeDataProcessor", OscilloscopeDataProcessor);
