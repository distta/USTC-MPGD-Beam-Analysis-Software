#include "Script/TimeResolutionScript.h"

#include "Algorithm/Oscilloscope/OscilloscopeDataProcessor.h"
#include "Algorithm/Analyzer/WaveformProcessor.h"
#include "Event/DataModel.h"
#include "Script/Base/RawDataParser.h"
#include "Script/Base/ScriptFactory.h"

#include <TFile.h>
#include <TDirectory.h>
#include <TF1.h>
#include <TH1D.h>
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
const array<string, kEstimatorCount> kEstimatorNames = {
    "XMean", "XMin", "YMean", "YMin", "XYMean", "XYMin"};

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

struct ReconstructionTimingResult {
    map<TrackKey, map<int, DetectorTimes>> trackTimes;
    set<uint64_t> wantedEventIDs;
    size_t missingRawDetectors = 0;
    size_t emptyTimingDetectors = 0;
};

// ==========================
// Reconstruction-based timing
// ==========================

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

bool LoadReconstructionTimes(TTree& validation, RawDataParser& parser,
                             ReconstructionTimingResult& result) {
    for (const char* branch : {"eventID", "rawEventID", "trackIndex", "detID",
                               "clusterIndices", "stripHits", "clusters"}) {
        if (!validation.GetBranch(branch)) {
            cerr << "[TimeResolution] branch " << branch
                 << " is missing; rerun the updated TrackAnalysis\n";
            return false;
        }
    }

    Int_t eventID = 0;
    ULong64_t rawEventID = 0;
    Int_t trackIndex = 0, detID = 0;
    vector<Int_t>* clusterIndices = nullptr;
    vector<StripHit>* stripHits = nullptr;
    vector<Cluster>* clusters = nullptr;
    validation.SetBranchAddress("eventID", &eventID);
    validation.SetBranchAddress("rawEventID", &rawEventID);
    validation.SetBranchAddress("trackIndex", &trackIndex);
    validation.SetBranchAddress("detID", &detID);
    validation.SetBranchAddress("clusterIndices", &clusterIndices);
    validation.SetBranchAddress("stripHits", &stripHits);
    validation.SetBranchAddress("clusters", &clusters);

    WaveformProcessor timingFitter;
    timingFitter.LoadConfig(json{{"mode", "Fit"}});
    int loadedEventID = numeric_limits<int>::min();
    unordered_map<int, vector<RawData>> rawHits;
    for (Long64_t entry = 0; entry < validation.GetEntries(); ++entry) {
        validation.GetEntry(entry);
        if (!clusterIndices || !stripHits || !clusters) continue;
        if (eventID != loadedEventID) {
            rawHits = parser.LoadEvent(eventID);
            loadedEventID = eventID;
        }
        const auto rawDetector = rawHits.find(detID);
        if (rawDetector == rawHits.end()) {
            ++result.missingRawDetectors;
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
        if (!hasAnyTime) ++result.emptyTimingDetectors;
        result.trackTimes[{rawEventID, trackIndex}][detID] = times;
        result.wantedEventIDs.insert(rawEventID);
    }
    return true;
}

// =================
// Output statistics
// =================

string SafeName(string value) {
    for (char& c : value)
        if (!isalnum(static_cast<unsigned char>(c))) c = '_';
    return value;
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

bool WriteTimingOutput(const string& outputPath,
                       const map<TrackKey, map<int, DetectorTimes>>& trackTimes,
                       const map<uint64_t, OscilloscopeData>& dataByEventID,
                       int histogramBins) {
    unique_ptr<TFile> output(TFile::Open(outputPath.c_str(), "RECREATE"));
    if (!output || output->IsZombie()) {
        cerr << "[TimeResolution] cannot create " << outputPath << '\n';
        return false;
    }

    TTree eventTree("EventTimes", "Matched tracker and oscilloscope times");
    ULong64_t outRawEventID = 0;
    Int_t outTrackIndex = 0, outDetID = 0, outFileIndex = 0, outSegmentIndex = 0;
    Double_t triggerTimeNs = 0.0;
    Double_t EventIDTime = 0.0, truthT0 = 0.0;
    array<Double_t, kEstimatorCount> estimatorValues{};
    eventTree.Branch("rawEventID", &outRawEventID);
    eventTree.Branch("trackIndex", &outTrackIndex);
    eventTree.Branch("detectorID", &outDetID);
    eventTree.Branch("traceFileIndex", &outFileIndex);
    eventTree.Branch("segmentIndex", &outSegmentIndex);
    eventTree.Branch("triggerTime", &triggerTimeNs);
    eventTree.Branch("EventIDTime", &EventIDTime);
    eventTree.Branch("truthT0", &truthT0);
    for (size_t estimator = 0; estimator < kEstimatorCount; ++estimator)
        eventTree.Branch(kEstimatorNames[estimator].c_str(), &estimatorValues[estimator]);

    map<pair<int, size_t>, vector<double>> detectorTimeSamples;
    map<tuple<int, int, size_t>, vector<double>> pairResiduals;
    map<pair<int, size_t>, vector<double>> truthResiduals;
    size_t tracksWithOscilloscopeData = 0;
    for (const auto& [key, detectors] : trackTimes) {
        const auto oscilloscopeDataIt = dataByEventID.find(key.first);
        const bool hasOscilloscopeData = oscilloscopeDataIt != dataByEventID.end();
        const OscilloscopeData* oscilloscopeData = hasOscilloscopeData ? &oscilloscopeDataIt->second : nullptr;
        if (hasOscilloscopeData) ++tracksWithOscilloscopeData;
        for (const auto& [detectorID, times] : detectors) {
            outRawEventID = key.first;
            outTrackIndex = key.second;
            outDetID = detectorID;
            outFileIndex = hasOscilloscopeData ? oscilloscopeData->fileIndex : -1;
            outSegmentIndex = hasOscilloscopeData ? oscilloscopeData->segmentIndex : -1;
            triggerTimeNs = hasOscilloscopeData ? oscilloscopeData->triggerTime : numeric_limits<double>::quiet_NaN();
            EventIDTime = hasOscilloscopeData
                                          ? oscilloscopeData->EvnentIDtime
                                          : numeric_limits<double>::quiet_NaN();
            truthT0 = hasOscilloscopeData ? oscilloscopeData->truthT0 : numeric_limits<double>::quiet_NaN();
            estimatorValues = times.value;
            eventTree.Fill();
            for (size_t estimator = 0; estimator < kEstimatorCount; ++estimator) {
                if (!isfinite(times.value[estimator])) continue;
                detectorTimeSamples[{detectorID, estimator}].push_back(times.value[estimator]);
                if (hasOscilloscopeData) {
                    truthResiduals[{detectorID, estimator}].push_back(
                        times.value[estimator] - oscilloscopeData->truthT0);
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
                            a - b);
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

    auto fillSummary = [&](const string& analysisName, const string& estimator,
                           int detA, int detB, const FitResult& fit,
                           double resolution) {
        analysis = analysisName;
        estimatorName = estimator;
        detectorA = detA;
        detectorB = detB;
        entries = fit.entries;
        meanNs = fit.mean;
        sigmaNs = fit.sigma;
        sigmaErrorNs = fit.sigmaError;
        resolutionNs = resolution;
        summary.Fill();
    };

    TDirectory* timeDirectory = output->mkdir("TimeDistributions");
    TDirectory* rawTimeDirectory = timeDirectory->mkdir("RawDetectorTime");
    for (const auto& [key, samples] : detectorTimeSamples) {
        const int detectorID = key.first;
        const size_t estimator = key.second;
        const string estimatorLabel = kEstimatorNames[estimator];
        const string name = SafeName("time_d" + to_string(detectorID) + "_" + estimatorLabel);
        const FitResult fit = WriteAndFit(
            samples, *rawTimeDirectory, name, "Detector time; t_{detector} [ns];Entries",
            histogramBins);
        fillSummary("detector_time", estimatorLabel, detectorID, -1, fit, fit.sigma);
    }

    TDirectory* pairDirectory = output->mkdir("PairResiduals");
    TDirectory* pairRawDirectory = pairDirectory->mkdir("Raw");
    for (const auto& [key, samples] : pairResiduals) {
        size_t estimator = 0;
        tie(detectorA, detectorB, estimator) = key;
        const string estimatorLabel = kEstimatorNames[estimator];
        const string name = SafeName("pair_d" + to_string(detectorA) + "_d" +
                                     to_string(detectorB) + "_" + estimatorLabel);
        const string title = "Tracker pair residual; t_{" + to_string(detectorA) +
                             "} - t_{" + to_string(detectorB) + "} [ns];Entries";
        const FitResult fit = WriteAndFit(samples, *pairRawDirectory, name, title,
                                          histogramBins);
        fillSummary("tracker_pair", estimatorLabel, detectorA, detectorB, fit,
                    fit.sigma / sqrt(2.0));
    }

    if (tracksWithOscilloscopeData > 0) {
        TDirectory* truthDirectory = output->mkdir("TruthResiduals");
        TDirectory* truthRawDirectory = truthDirectory->mkdir("Raw");
        for (const auto& [key, samples] : truthResiduals) {
            const int detectorID = key.first;
            const size_t estimator = key.second;
            const string estimatorLabel = kEstimatorNames[estimator];
            const string name = SafeName("truth_d" + to_string(detectorID) + "_" + estimatorLabel);
            const FitResult fit = WriteAndFit(
                samples, *truthRawDirectory, name,
                "External-T0 corrected time; t_{detector} - T0 [ns];Entries",
                histogramBins);
            fillSummary("truth_t0", estimatorLabel, detectorID, -1, fit, fit.sigma);
        }
    }

    output->cd();
    summary.Write();
    cout << "[TimeResolution] tracks=" << trackTimes.size()
         << ", with external T0=" << tracksWithOscilloscopeData
         << ", output=" << outputPath << '\n';
    return !trackTimes.empty();
}

}  // namespace

void TimeResolutionScript::LoadConfig(const json& config) {
    m_csvDirectory = config.value("csvDirectory", m_csvDirectory);
    m_trackFile = config.value("trackFile", m_trackFile);
    m_outputFile = config.value("outputFile", m_outputFile);
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
    
    ReconstructionTimingResult reconstruction;
    if (!LoadReconstructionTimes(*validation, *parser, reconstruction)) return false;
    cout << "[TimeResolution] loaded " << reconstruction.trackTimes.size() << " tracks and "
         << reconstruction.wantedEventIDs.size() << " raw event IDs\n";
    if (reconstruction.missingRawDetectors > 0 || reconstruction.emptyTimingDetectors > 0)
        cout << "[TimeResolution] Fit timing diagnostics: missing raw detector entries="
             << reconstruction.missingRawDetectors
             << ", empty detector times=" << reconstruction.emptyTimingDetectors << '\n';

    OscilloscopeDataProcessor oscilloscopeProcessor;
    oscilloscopeProcessor.LoadConfig(json{
        {"csvDirectory", m_csvDirectory},
        {"outputDir", GetOutputDir()},
        {"dataCacheFile", m_dataCacheFile},
        {"waveformDiagnosticFile", m_waveformDiagnosticFile},
        {"rebuildTimingCache", m_rebuildTimingCache},
        {"writeWaveformDiagnostics", m_writeWaveformDiagnostics},
        {"triggerThreshold", m_triggerThreshold},
        {"eventThreshold", m_eventThreshold},
        {"eventBits", m_eventBits},
        {"maxWaveformFiles", m_maxWaveformFiles},
    });
    OscilloscopeDataResult oscilloscope;
    if (oscilloscopeProcessor.Initialize())
        oscilloscope = oscilloscopeProcessor.LoadData(reconstruction.wantedEventIDs);
    const auto& trackTimes = reconstruction.trackTimes;
    const auto& dataByEventID = oscilloscope.dataByEventID;

    return WriteTimingOutput(outputPath, trackTimes, dataByEventID, m_histogramBins);
}

REGISTER_SCRIPT("TimeResolution", TimeResolutionScript);
