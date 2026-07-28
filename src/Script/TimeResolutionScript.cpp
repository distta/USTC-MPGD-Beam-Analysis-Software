#include "Script/TimeResolutionScript.h"

#include "Algorithm/Analyzer/HitProcessor.h"
#include "Detector/DetectorFactory.h"
#include "Event/DataModel.h"
#include "Script/Base/RawDataParser.h"
#include "Script/Base/ScriptFactory.h"
#include "Terminal.h"

#include <TFile.h>
#include <TCanvas.h>
#include <TDirectory.h>
#include <TF1.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLine.h>
#include <TMarker.h>
#include <TProfile.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <future>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace std;

namespace {

constexpr size_t kEstimatorCount = 6;
constexpr size_t kTrackerTimingEstimator = 4;
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
    map<TrackKey, int> eventIDs;
    set<uint64_t> wantedEventIDs;
    size_t missingRawDetectors = 0;
    size_t emptyTimingDetectors = 0;
};

struct TrackTimeReference {
    uint64_t rawEventID = 0;
    int trackIndex = -1;
    double time = numeric_limits<double>::quiet_NaN();
};

struct DUTTimingSample {
    int eventID = -1;
    uint64_t rawEventID = 0;
    int trackIndex = -1;
    int detectorID = -1;
    double amplitude = numeric_limits<double>::quiet_NaN();
    double dutTime = numeric_limits<double>::quiet_NaN();
    double trackTime = numeric_limits<double>::quiet_NaN();
    double residual = numeric_limits<double>::quiet_NaN();
    double meanAmplitude = numeric_limits<double>::quiet_NaN();
    double clusterCharge = numeric_limits<double>::quiet_NaN();
    double clusterMaxAmplitude = numeric_limits<double>::quiet_NaN();
    double clusterSize = numeric_limits<double>::quiet_NaN();
    double clusterCentroid = numeric_limits<double>::quiet_NaN();
    double clusterLocalX = numeric_limits<double>::quiet_NaN();
    double clusterLocalY = numeric_limits<double>::quiet_NaN();
    double predictedX = numeric_limits<double>::quiet_NaN();
    double predictedY = numeric_limits<double>::quiet_NaN();
    bool insideActiveArea = false;
};

struct DUTTimingResult {
    map<int, vector<DUTTimingSample>> samplesByDetector;
    map<int, set<TrackKey>> activeAreaTrackCases;
    map<int, set<TrackKey>> activeAreaMatchedHitCases;
    size_t unmatchedTrackTimes = 0;
    size_t invalidDUTTimes = 0;
};

struct TrackTimeWeights {
    array<double, 3> value{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
    array<double, 3> detectorVariance{
        numeric_limits<double>::quiet_NaN(),
        numeric_limits<double>::quiet_NaN(),
        numeric_limits<double>::quiet_NaN()};
    double resolution = numeric_limits<double>::quiet_NaN();
    bool valid = false;
};

struct OscilloscopeT0Reference {
    double time = numeric_limits<double>::quiet_NaN();
    double resolution = numeric_limits<double>::quiet_NaN();
    int validChannels = 0;
    int traceFileIndex = -1;
    int segmentIndex = -1;
};

struct OscilloscopeT0Result {
    map<uint64_t, OscilloscopeT0Reference> references;
    array<double, 3> channelResolution{
        numeric_limits<double>::quiet_NaN(),
        numeric_limits<double>::quiet_NaN(),
        numeric_limits<double>::quiet_NaN()};
    array<double, 3> channelResolutionError{
        numeric_limits<double>::quiet_NaN(),
        numeric_limits<double>::quiet_NaN(),
        numeric_limits<double>::quiet_NaN()};
    array<double, 3> channelOffset{0.0, 0.0, 0.0};
    array<double, 3> channelWeight{
        numeric_limits<double>::quiet_NaN(),
        numeric_limits<double>::quiet_NaN(),
        numeric_limits<double>::quiet_NaN()};
    double resolution = numeric_limits<double>::quiet_NaN();
    size_t inputEntries = 0;
    size_t invalidEventIDs = 0;
    size_t duplicateEventIDs = 0;
};

struct T0ResidualSamples {
    vector<double> residual;
    double referenceVarianceSum = 0.0;

    void Add(double value, double referenceResolution) {
        if (!isfinite(value) || !isfinite(referenceResolution) ||
            referenceResolution < 0.0)
            return;
        residual.push_back(value);
        referenceVarianceSum +=
            referenceResolution * referenceResolution;
    }

    double MeanReferenceVariance() const {
        return residual.empty()
                   ? numeric_limits<double>::quiet_NaN()
                   : referenceVarianceSum /
                         static_cast<double>(residual.size());
    }
};

struct TimeMeasurement {
    double time = numeric_limits<double>::quiet_NaN();
    double error = numeric_limits<double>::quiet_NaN();
};

class TimingProgress {
   public:
    TimingProgress(string label, Long64_t total,
                   const size_t& fitCount)
        : m_label(move(label)), m_total(total), m_fitCount(fitCount) {}

    void Update(Long64_t completed) {
        if (!Terminal::Interactive()) return;
        const int percent =
            m_total > 0
                ? min(100, static_cast<int>(
                               100 * completed / m_total))
                : 100;
        const int displayedPercent =
            completed >= m_total ? 100 : 5 * (percent / 5);
        if (displayedPercent == m_lastPercent && completed < m_total) return;
        m_lastPercent = displayedPercent;

        ostringstream line;
        line << '\r' << "      "
             << left << setw(15) << m_label
             << right << setw(3) << displayedPercent << "%  "
             << completed << '/' << m_total
             << "  fits=" << m_fitCount;
        cout << line.str() << flush;
        if (completed >= m_total) Terminal::ClearProgress();
    }

   private:
    string m_label;
    Long64_t m_total = 0;
    const size_t& m_fitCount;
    int m_lastPercent = -1;
};

double ErrorWeightedMean(const vector<TimeMeasurement>& measurements) {
    double weightedTime = 0.0;
    double weightSum = 0.0;
    for (const TimeMeasurement& measurement : measurements) {
        if (!isfinite(measurement.time) || !isfinite(measurement.error) ||
            measurement.error <= 0.0)
            continue;
        const double weight =
            1.0 / (measurement.error * measurement.error);
        if (!isfinite(weight)) continue;
        weightedTime += weight * measurement.time;
        weightSum += weight;
    }
    if (weightSum > 0.0) return weightedTime / weightSum;

    double timeSum = 0.0;
    size_t validTimes = 0;
    for (const TimeMeasurement& measurement : measurements) {
        if (!isfinite(measurement.time)) continue;
        timeSum += measurement.time;
        ++validTimes;
    }
    return validTimes > 0
               ? timeSum / static_cast<double>(validTimes)
               : numeric_limits<double>::quiet_NaN();
}

// ==========================
// Reconstruction-based timing
// ==========================

array<double, kEstimatorCount> ExtractTimes(const vector<int>& clusterIndices,
                                            const vector<Cluster>& clusters,
                                            const vector<ChannelHit>& channelHits,
                                            const vector<RawData>& detectorRawData,
                                            HitProcessor& timingFitter,
                                            size_t& fitCount) {
    array<double, kEstimatorCount> result;
    result.fill(numeric_limits<double>::quiet_NaN());
    array<vector<TimeMeasurement>, 2> planeStripTimes;
    map<int, TimeMeasurement> fittedStripTimes;

    for (int clusterIndex : clusterIndices) {
        if (clusterIndex < 0 || clusterIndex >= static_cast<int>(clusters.size())) continue;
        const Cluster& cluster = clusters[clusterIndex];
        if (cluster.type < 0 || cluster.type > 1) continue;

        vector<TimeMeasurement> stripTimes;
        for (int stripIndex : cluster.channelHitIndices) {
            if (stripIndex < 0 || stripIndex >= static_cast<int>(channelHits.size())) continue;
            const int rawIndex = channelHits[stripIndex].rawIndices;
            if (rawIndex < 0 || rawIndex >= static_cast<int>(detectorRawData.size())) continue;

            auto cached = fittedStripTimes.find(rawIndex);
            if (cached == fittedStripTimes.end()) {
                const ChannelHit fitted =
                    timingFitter.ProcessHit(detectorRawData[rawIndex]);
                ++fitCount;
                TimeMeasurement measurement;
                if (fitted.isValid && isfinite(fitted.time)) {
                    measurement.time = fitted.time;
                    measurement.error = fitted.timeError;
                }
                cached =
                    fittedStripTimes.emplace(rawIndex, measurement).first;
            }
            if (isfinite(cached->second.time))
                stripTimes.push_back(cached->second);
        }
        if (stripTimes.empty()) continue;
        planeStripTimes[cluster.type].insert(
            planeStripTimes[cluster.type].end(), stripTimes.begin(), stripTimes.end());
    }

    vector<TimeMeasurement> combined;
    for (size_t plane = 0; plane < 2; ++plane) {
        const auto& times = planeStripTimes[plane];
        if (times.empty()) continue;
        result[plane * 2] = ErrorWeightedMean(times);
        result[plane * 2 + 1] =
            min_element(
                times.begin(), times.end(),
                [](const TimeMeasurement& first,
                   const TimeMeasurement& second) {
                    return first.time < second.time;
                })
                ->time;
        combined.insert(combined.end(), times.begin(), times.end());
    }
    if (!combined.empty()) {
        result[4] = ErrorWeightedMean(combined);
        result[5] =
            min_element(
                combined.begin(), combined.end(),
                [](const TimeMeasurement& first,
                   const TimeMeasurement& second) {
                    return first.time < second.time;
                })
                ->time;
    }
    return result;
}

bool LoadReconstructionTimes(TTree& validation, RawDataParser& parser,
                             const json& timingWaveformConfig,
                             size_t& fitCount,
                             ReconstructionTimingResult& result) {
    for (const char* branch : {"eventID", "rawEventID", "trackIndex", "detID",
                               "clusterIndices", "channelHits", "clusters"}) {
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
    vector<ChannelHit>* channelHits = nullptr;
    vector<Cluster>* clusters = nullptr;
    validation.SetBranchAddress("eventID", &eventID);
    validation.SetBranchAddress("rawEventID", &rawEventID);
    validation.SetBranchAddress("trackIndex", &trackIndex);
    validation.SetBranchAddress("detID", &detID);
    validation.SetBranchAddress("clusterIndices", &clusterIndices);
    validation.SetBranchAddress("channelHits", &channelHits);
    validation.SetBranchAddress("clusters", &clusters);

    HitProcessor timingFitter;
    timingFitter.LoadConfig(timingWaveformConfig);
    int loadedEventID = numeric_limits<int>::min();
    unordered_map<int, vector<RawData>> rawHits;
    const Long64_t entries = validation.GetEntries();
    TimingProgress progress("Tracker timing", entries, fitCount);
    for (Long64_t entry = 0; entry < entries; ++entry) {
        progress.Update(entry);
        validation.GetEntry(entry);
        if (!clusterIndices || !channelHits || !clusters) continue;
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
        times.value = ExtractTimes(*clusterIndices, *clusters, *channelHits,
                                   rawDetector->second, timingFitter,
                                   fitCount);
        bool hasAnyTime = false;
        for (double value : times.value) {
            if (isfinite(value)) {
                hasAnyTime = true;
                break;
            }
        }
        if (!hasAnyTime) ++result.emptyTimingDetectors;
        result.trackTimes[{rawEventID, trackIndex}][detID] = times;
        result.eventIDs[{rawEventID, trackIndex}] = eventID;
        result.wantedEventIDs.insert(rawEventID);
    }
    progress.Update(entries);
    return true;
}

map<int, TrackTimeReference> BuildTrackTimeReferences(
    const map<TrackKey, map<int, DetectorTimes>>& trackTimes,
    const map<TrackKey, int>& eventIDs,
    const array<int, 3>& trackerIDs,
    const TrackTimeWeights& weights) {
    map<int, TrackTimeReference> references;
    set<int> ambiguousEventIDs;
    for (const auto& [key, detectors] : trackTimes) {
        array<double, 3> times{};
        bool valid = true;
        for (size_t i = 0; i < trackerIDs.size(); ++i) {
            const auto detector = detectors.find(trackerIDs[i]);
            if (detector == detectors.end() ||
                !isfinite(detector->second.value[kTrackerTimingEstimator])) {
                valid = false;
                break;
            }
            times[i] =
                detector->second.value[kTrackerTimingEstimator];
        }
        const auto event = eventIDs.find(key);
        if (!valid || event == eventIDs.end()) continue;
        if (references.count(event->second)) {
            ambiguousEventIDs.insert(event->second);
            continue;
        }
        double trackTime = 0.0;
        for (size_t i = 0; i < trackerIDs.size(); ++i) {
            trackTime += weights.value[i] * times[i];
        }
        references[event->second] = {key.first, key.second, trackTime};
    }
    for (int eventID : ambiguousEventIDs) references.erase(eventID);
    return references;
}

DUTTimingResult LoadDUTTiming(
    TTree& dutTree, RawDataParser& parser,
    const map<int, TrackTimeReference>& trackReferences,
    const json& timingWaveformConfig,
    size_t& fitCount) {
    DUTTimingResult result;
    for (const char* branch :
         {"eventID", "dutID", "hitFlag", "predX", "predY",
          "selectedCluster", "selectedChannelHits"}) {
        if (!dutTree.GetBranch(branch)) {
            cerr << "[TimeResolution] PadDUTTree branch " << branch
                 << " is missing\n";
            return result;
        }
    }

    Int_t eventID = 0, dutID = 0, hitFlag = 0;
    Double_t predictedX = numeric_limits<double>::quiet_NaN();
    Double_t predictedY = numeric_limits<double>::quiet_NaN();
    Cluster* selectedCluster = nullptr;
    vector<ChannelHit>* selectedChannelHits = nullptr;
    dutTree.SetBranchAddress("eventID", &eventID);
    dutTree.SetBranchAddress("dutID", &dutID);
    dutTree.SetBranchAddress("hitFlag", &hitFlag);
    dutTree.SetBranchAddress("predX", &predictedX);
    dutTree.SetBranchAddress("predY", &predictedY);
    dutTree.SetBranchAddress("selectedCluster", &selectedCluster);
    dutTree.SetBranchAddress("selectedChannelHits", &selectedChannelHits);

    HitProcessor timingFitter;
    timingFitter.LoadConfig(timingWaveformConfig);
    int loadedEventID = numeric_limits<int>::min();
    unordered_map<int, vector<RawData>> rawHits;
    const Long64_t entries = dutTree.GetEntries();
    TimingProgress progress("DUT timing", entries, fitCount);
    for (Long64_t entry = 0; entry < entries; ++entry) {
        progress.Update(entry);
        dutTree.GetEntry(entry);
        const bool hasMatchedHit =
            hitFlag != 0 && selectedCluster && selectedChannelHits &&
            !selectedChannelHits->empty();
        const auto reference = trackReferences.find(eventID);
        if (reference == trackReferences.end()) {
            if (hasMatchedHit) ++result.unmatchedTrackTimes;
            continue;
        }
        bool insideActiveArea = false;
        const auto detector =
            DetectorFactory::GetInstance().GetDetector(dutID);
        const auto* pad = detector
                              ? detector->GetPlanarPadConfig()
                              : nullptr;
        if (pad && isfinite(predictedX) && isfinite(predictedY)) {
            const double xMinimum = -0.5 * pad->pitchX;
            const double xMaximum =
                (pad->columns - 0.5) * pad->pitchX;
            const double yMinimum = -0.5 * pad->pitchY;
            const double yMaximum =
                (pad->rows - 0.5) * pad->pitchY;
            insideActiveArea =
                predictedX >= xMinimum && predictedX < xMaximum &&
                predictedY >= yMinimum && predictedY < yMaximum;
        }
        const TrackKey trackCase{
            reference->second.rawEventID,
            reference->second.trackIndex};
        if (insideActiveArea) {
            result.activeAreaTrackCases[dutID].insert(trackCase);
            if (hasMatchedHit)
                result.activeAreaMatchedHitCases[dutID].insert(
                    trackCase);
        }
        if (!hasMatchedHit) continue;
        if (eventID != loadedEventID) {
            rawHits = parser.LoadEvent(eventID);
            loadedEventID = eventID;
        }
        const auto rawDetector = rawHits.find(dutID);
        if (rawDetector == rawHits.end()) {
            ++result.invalidDUTTimes;
            continue;
        }

        vector<TimeMeasurement> dutMeasurements;
        double amplitudeSum = 0.0;
        double maximumAmplitude = -numeric_limits<double>::infinity();
        size_t validDUTChannels = 0;
        for (const ChannelHit& selectedHit : *selectedChannelHits) {
            const int rawIndex = selectedHit.rawIndices;
            if (rawIndex < 0 ||
                rawIndex >= static_cast<int>(rawDetector->second.size())) {
                continue;
            }
            const ChannelHit fitted =
                timingFitter.ProcessHit(rawDetector->second[rawIndex]);
            ++fitCount;
            if (fitted.isValid && isfinite(fitted.time)) {
                dutMeasurements.push_back(
                    {fitted.time, fitted.timeError});
                amplitudeSum += fitted.amp;
                maximumAmplitude = max(maximumAmplitude, fitted.amp);
                ++validDUTChannels;
            }
        }
        if (validDUTChannels == 0) {
            ++result.invalidDUTTimes;
            continue;
        }
        const double dutTime = ErrorWeightedMean(dutMeasurements);
        DUTTimingSample sample;
        sample.eventID = eventID;
        sample.rawEventID = reference->second.rawEventID;
        sample.trackIndex = reference->second.trackIndex;
        sample.detectorID = dutID;
        sample.amplitude = maximumAmplitude;
        sample.dutTime = dutTime;
        sample.trackTime = reference->second.time;
        sample.residual = dutTime - reference->second.time;
        sample.meanAmplitude =
            amplitudeSum / static_cast<double>(validDUTChannels);
        sample.clusterCharge = selectedCluster->charge;
        sample.clusterMaxAmplitude = maximumAmplitude;
        sample.clusterSize = selectedCluster->size;
        sample.clusterCentroid = selectedCluster->centroid;
        sample.clusterLocalX = selectedCluster->hasLocalPosition
                                   ? selectedCluster->localPosition.X()
                                   : selectedCluster->pos;
        sample.clusterLocalY = selectedCluster->hasLocalPosition
                                   ? selectedCluster->localPosition.Y()
                                   : numeric_limits<double>::quiet_NaN();
        sample.predictedX = predictedX;
        sample.predictedY = predictedY;
        sample.insideActiveArea = insideActiveArea;
        result.samplesByDetector[dutID].push_back(move(sample));
    }
    progress.Update(entries);
    return result;
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

FitResult FitGaussianWithoutWriting(const vector<double>& values) {
    FitResult result;
    result.entries = static_cast<long long>(values.size());
    if (values.size() < 3) return result;
    const double mean =
        accumulate(values.begin(), values.end(), 0.0) / values.size();
    double variance = 0.0;
    for (double value : values)
        variance += (value - mean) * (value - mean);
    const double sigma = sqrt(variance / (values.size() - 1));
    if (!isfinite(sigma) || sigma <= 0.0) return result;

    static unsigned long fitCounter = 0;
    const string suffix = to_string(fitCounter++);
    TH1D histogram(
        ("weight_fit_hist_" + suffix).c_str(), "",
        160, mean - 6.0 * sigma, mean + 6.0 * sigma);
    histogram.SetDirectory(nullptr);
    for (double value : values) histogram.Fill(value);
    TF1 gaussian(
        ("weight_fit_gaus_" + suffix).c_str(), "gaus",
        mean - 2.0 * sigma, mean + 2.0 * sigma);
    gaussian.SetParameters(histogram.GetMaximum(), mean, sigma);
    const int status = histogram.Fit(&gaussian, "Q");
    if (status == 0 && isfinite(gaussian.GetParameter(2))) {
        result.mean = gaussian.GetParameter(1);
        result.sigma = abs(gaussian.GetParameter(2));
        result.sigmaError = gaussian.GetParError(2);
    } else {
        result.mean = mean;
        result.sigma = sigma;
    }
    return result;
}

double Median(vector<double> values) {
    values.erase(remove_if(values.begin(), values.end(),
                           [](double value) { return !isfinite(value); }),
                 values.end());
    if (values.empty())
        return numeric_limits<double>::quiet_NaN();
    const size_t middle = values.size() / 2;
    nth_element(values.begin(), values.begin() + middle, values.end());
    const double upper = values[middle];
    if (values.size() % 2 != 0) return upper;
    nth_element(values.begin(), values.begin() + middle - 1,
                values.begin() + middle);
    return 0.5 * (values[middle - 1] + upper);
}

OscilloscopeT0Result LoadOscilloscopeT0(
    const string& inputPath, const set<uint64_t>& wantedEventIDs) {
    OscilloscopeT0Result result;
    unique_ptr<TFile> input(TFile::Open(inputPath.c_str(), "READ"));
    if (!input || input->IsZombie()) {
        cerr << "[TimeResolution] cannot open oscilloscope file "
             << inputPath << "; oscilloscope T0 analysis will be skipped\n";
        return result;
    }

    auto* compactEvents = dynamic_cast<TTree*>(input->Get("Events"));
    if (compactEvents) {
        for (const char* branch :
             {"eventID", "referenceTime", "C1Time", "C2Time", "C3Time"}) {
            if (!compactEvents->GetBranch(branch)) {
                cerr << "[TimeResolution] compact oscilloscope branch "
                     << branch
                     << " is missing; T0 analysis will be skipped\n";
                return result;
            }
        }

        const array<string, 3> histogramNames = {
            "hDeltaT_C1_C2", "hDeltaT_C1_C3", "hDeltaT_C2_C3"};
        array<FitResult, 3> pairFits;
        for (size_t pair = 0; pair < histogramNames.size(); ++pair) {
            auto* histogram =
                dynamic_cast<TH1D*>(input->Get(histogramNames[pair].c_str()));
            const string fitName = "fit_" + histogramNames[pair];
            TF1* fit = histogram
                           ? histogram->GetFunction(fitName.c_str())
                           : nullptr;
            if (!fit || !isfinite(fit->GetParameter(2)) ||
                fit->GetParameter(2) == 0.0) {
                cerr << "[TimeResolution] Gaussian fit " << fitName
                     << " is missing; T0 analysis will be skipped\n";
                return result;
            }
            pairFits[pair].mean = fit->GetParameter(1);
            pairFits[pair].sigma = abs(fit->GetParameter(2));
            pairFits[pair].sigmaError = fit->GetParError(2);
        }

        const array<double, 3> variance = {
            0.5 * (pairFits[0].sigma * pairFits[0].sigma +
                   pairFits[1].sigma * pairFits[1].sigma -
                   pairFits[2].sigma * pairFits[2].sigma),
            0.5 * (pairFits[0].sigma * pairFits[0].sigma +
                   pairFits[2].sigma * pairFits[2].sigma -
                   pairFits[1].sigma * pairFits[1].sigma),
            0.5 * (pairFits[1].sigma * pairFits[1].sigma +
                   pairFits[2].sigma * pairFits[2].sigma -
                   pairFits[0].sigma * pairFits[0].sigma)};
        double inverseVarianceSum = 0.0;
        for (size_t channel = 0; channel < variance.size(); ++channel) {
            if (!isfinite(variance[channel]) || variance[channel] <= 0.0) {
                cerr << "[TimeResolution] invalid compact oscilloscope C"
                     << channel + 1
                     << " resolution; T0 analysis will be skipped\n";
                return result;
            }
            result.channelResolution[channel] = sqrt(variance[channel]);
            inverseVarianceSum += 1.0 / variance[channel];

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
                    pairFits[term].sigma * pairFits[term].sigmaError;
                varianceErrorSquared += contribution * contribution;
            }
            result.channelResolutionError[channel] =
                sqrt(varianceErrorSquared) /
                (2.0 * result.channelResolution[channel]);
        }
        result.resolution = 1.0 / sqrt(inverseVarianceSum);
        for (size_t channel = 0; channel < variance.size(); ++channel)
            result.channelWeight[channel] =
                (1.0 / variance[channel]) / inverseVarianceSum;

        ULong64_t eventID = 0;
        Double_t referenceTime =
            numeric_limits<double>::quiet_NaN();
        compactEvents->SetBranchAddress("eventID", &eventID);
        compactEvents->SetBranchAddress(
            "referenceTime", &referenceTime);
        result.inputEntries =
            static_cast<size_t>(compactEvents->GetEntries());
        set<uint64_t> ambiguousEventIDs;
        for (Long64_t entry = 0;
             entry < compactEvents->GetEntries(); ++entry) {
            compactEvents->GetEntry(entry);
            if (!wantedEventIDs.count(eventID) ||
                !isfinite(referenceTime) ||
                ambiguousEventIDs.count(eventID))
                continue;
            OscilloscopeT0Reference reference;
            reference.time = referenceTime;
            reference.resolution = result.resolution;
            reference.validChannels = 3;
            const auto inserted =
                result.references.emplace(eventID, reference);
            if (!inserted.second) {
                result.references.erase(eventID);
                ambiguousEventIDs.insert(eventID);
            }
        }
        result.duplicateEventIDs = ambiguousEventIDs.size();
        return result;
    }

    auto* events =
        dynamic_cast<TTree*>(input->Get("OscilloscopeEvents"));
    auto* resolutionSummary =
        dynamic_cast<TTree*>(input->Get("TimingResolution"));
    if (!events || !resolutionSummary) {
        cerr << "[TimeResolution] OscilloscopeEvents or TimingResolution is "
                "missing in "
             << inputPath
             << "; oscilloscope T0 analysis will be skipped\n";
        return result;
    }

    array<Double_t, 3> channelResolution{};
    array<Double_t, 3> channelResolutionError{};
    array<Bool_t, 3> channelResolutionValid{};
    for (size_t channel = 0; channel < channelResolution.size(); ++channel) {
        const string suffix = "C" + to_string(channel + 1);
        const string resolutionBranch = "resolution" + suffix;
        const string errorBranch = "resolutionError" + suffix;
        const string validBranch = "resolution" + suffix + "Valid";
        if (!resolutionSummary->GetBranch(resolutionBranch.c_str()) ||
            !resolutionSummary->GetBranch(errorBranch.c_str()) ||
            !resolutionSummary->GetBranch(validBranch.c_str())) {
            cerr << "[TimeResolution] oscilloscope resolution branch for "
                 << suffix << " is missing; T0 analysis will be skipped\n";
            return result;
        }
        resolutionSummary->SetBranchAddress(
            resolutionBranch.c_str(), &channelResolution[channel]);
        resolutionSummary->SetBranchAddress(
            errorBranch.c_str(), &channelResolutionError[channel]);
        resolutionSummary->SetBranchAddress(
            validBranch.c_str(), &channelResolutionValid[channel]);
    }
    if (resolutionSummary->GetEntries() < 1) {
        cerr << "[TimeResolution] oscilloscope TimingResolution is empty; "
                "T0 analysis will be skipped\n";
        return result;
    }
    resolutionSummary->GetEntry(0);

    double inverseVarianceSum = 0.0;
    for (size_t channel = 0; channel < channelResolution.size(); ++channel) {
        if (!channelResolutionValid[channel] ||
            !isfinite(channelResolution[channel]) ||
            channelResolution[channel] <= 0.0) {
            cerr << "[TimeResolution] invalid C" << channel + 1
                 << " oscilloscope resolution; T0 analysis will be skipped\n";
            return result;
        }
        result.channelResolution[channel] = channelResolution[channel];
        result.channelResolutionError[channel] =
            channelResolutionError[channel];
        inverseVarianceSum +=
            1.0 / (channelResolution[channel] * channelResolution[channel]);
    }
    result.resolution = 1.0 / sqrt(inverseVarianceSum);
    for (size_t channel = 0; channel < channelResolution.size(); ++channel)
        result.channelWeight[channel] =
            (1.0 / (channelResolution[channel] *
                    channelResolution[channel])) /
            inverseVarianceSum;

    for (const char* branch :
         {"eventID", "eventIDValid", "traceFileIndex", "segmentIndex",
          "C1TimeValid", "C1Time", "C2TimeValid", "C2Time",
          "C3TimeValid", "C3Time"}) {
        if (!events->GetBranch(branch)) {
            cerr << "[TimeResolution] oscilloscope branch " << branch
                 << " is missing; T0 analysis will be skipped\n";
            return OscilloscopeT0Result{};
        }
    }

    struct ScopeEvent {
        uint64_t eventID = 0;
        bool eventIDValid = false;
        int traceFileIndex = -1;
        int segmentIndex = -1;
        array<bool, 3> timeValid{};
        array<double, 3> time{};
    };
    vector<ScopeEvent> scopeEvents;
    scopeEvents.reserve(static_cast<size_t>(events->GetEntries()));
    array<vector<double>, 2> pairDifferences;

    ULong64_t eventID = 0;
    Bool_t eventIDValid = false;
    Int_t traceFileIndex = -1, segmentIndex = -1;
    array<Bool_t, 3> timeValid{};
    array<Double_t, 3> time{};
    events->SetBranchAddress("eventID", &eventID);
    events->SetBranchAddress("eventIDValid", &eventIDValid);
    events->SetBranchAddress("traceFileIndex", &traceFileIndex);
    events->SetBranchAddress("segmentIndex", &segmentIndex);
    for (size_t channel = 0; channel < time.size(); ++channel) {
        const string prefix = "C" + to_string(channel + 1);
        events->SetBranchAddress(
            (prefix + "TimeValid").c_str(), &timeValid[channel]);
        events->SetBranchAddress(
            (prefix + "Time").c_str(), &time[channel]);
    }

    result.inputEntries = static_cast<size_t>(events->GetEntries());
    for (Long64_t entry = 0; entry < events->GetEntries(); ++entry) {
        events->GetEntry(entry);
        ScopeEvent scopeEvent;
        scopeEvent.eventID = eventID;
        scopeEvent.eventIDValid = eventIDValid;
        scopeEvent.traceFileIndex = traceFileIndex;
        scopeEvent.segmentIndex = segmentIndex;
        for (size_t channel = 0; channel < time.size(); ++channel) {
            scopeEvent.timeValid[channel] =
                timeValid[channel] && isfinite(time[channel]);
            scopeEvent.time[channel] = time[channel];
        }
        if (!scopeEvent.eventIDValid) ++result.invalidEventIDs;
        if (scopeEvent.timeValid[0] && scopeEvent.timeValid[1])
            pairDifferences[0].push_back(
                scopeEvent.time[0] - scopeEvent.time[1]);
        if (scopeEvent.timeValid[0] && scopeEvent.timeValid[2])
            pairDifferences[1].push_back(
                scopeEvent.time[0] - scopeEvent.time[2]);
        scopeEvents.push_back(scopeEvent);
    }

    result.channelOffset[1] = Median(move(pairDifferences[0]));
    result.channelOffset[2] = Median(move(pairDifferences[1]));
    if (!isfinite(result.channelOffset[1]) ||
        !isfinite(result.channelOffset[2])) {
        cerr << "[TimeResolution] cannot align oscilloscope channel offsets; "
                "T0 analysis will be skipped\n";
        return OscilloscopeT0Result{};
    }

    set<uint64_t> ambiguousEventIDs;
    for (const ScopeEvent& scopeEvent : scopeEvents) {
        if (!scopeEvent.eventIDValid ||
            !wantedEventIDs.count(scopeEvent.eventID) ||
            ambiguousEventIDs.count(scopeEvent.eventID))
            continue;
        double weightedTime = 0.0;
        double eventInverseVarianceSum = 0.0;
        int validChannels = 0;
        for (size_t channel = 0; channel < scopeEvent.time.size();
             ++channel) {
            if (!scopeEvent.timeValid[channel]) continue;
            const double inverseVariance =
                1.0 / (result.channelResolution[channel] *
                       result.channelResolution[channel]);
            weightedTime +=
                inverseVariance *
                (scopeEvent.time[channel] +
                 result.channelOffset[channel]);
            eventInverseVarianceSum += inverseVariance;
            ++validChannels;
        }
        if (eventInverseVarianceSum <= 0.0) continue;
        OscilloscopeT0Reference reference;
        reference.time = weightedTime / eventInverseVarianceSum;
        reference.resolution = 1.0 / sqrt(eventInverseVarianceSum);
        reference.validChannels = validChannels;
        reference.traceFileIndex = scopeEvent.traceFileIndex;
        reference.segmentIndex = scopeEvent.segmentIndex;
        const auto inserted =
            result.references.emplace(scopeEvent.eventID, reference);
        if (!inserted.second) {
            result.references.erase(scopeEvent.eventID);
            ambiguousEventIDs.insert(scopeEvent.eventID);
        }
    }
    result.duplicateEventIDs = ambiguousEventIDs.size();
    return result;
}


TrackTimeWeights CalculateTrackTimeWeights(
    const map<TrackKey, map<int, DetectorTimes>>& trackTimes,
    const array<int, 3>& trackerIDs) {
    array<vector<double>, 3> pairResiduals;
    for (const auto& [key, detectors] : trackTimes) {
        (void)key;
        array<double, 3> times{};
        bool complete = true;
        for (size_t i = 0; i < trackerIDs.size(); ++i) {
            const auto detector = detectors.find(trackerIDs[i]);
            if (detector == detectors.end()) {
                complete = false;
                break;
            }
            times[i] =
                detector->second.value[kTrackerTimingEstimator];
            if (!isfinite(times[i])) {
                complete = false;
                break;
            }
        }
        if (!complete) continue;
        pairResiduals[0].push_back(times[0] - times[1]);
        pairResiduals[1].push_back(times[0] - times[2]);
        pairResiduals[2].push_back(times[1] - times[2]);
    }
    const FitResult fit12 = FitGaussianWithoutWriting(pairResiduals[0]);
    const FitResult fit13 = FitGaussianWithoutWriting(pairResiduals[1]);
    const FitResult fit23 = FitGaussianWithoutWriting(pairResiduals[2]);
    TrackTimeWeights result;
    if (!isfinite(fit12.sigma) || !isfinite(fit13.sigma) ||
        !isfinite(fit23.sigma)) {
        return result;
    }
    result.detectorVariance[0] =
        0.5 * (fit12.sigma * fit12.sigma +
               fit13.sigma * fit13.sigma -
               fit23.sigma * fit23.sigma);
    result.detectorVariance[1] =
        0.5 * (fit12.sigma * fit12.sigma +
               fit23.sigma * fit23.sigma -
               fit13.sigma * fit13.sigma);
    result.detectorVariance[2] =
        0.5 * (fit13.sigma * fit13.sigma +
               fit23.sigma * fit23.sigma -
               fit12.sigma * fit12.sigma);
    double inverseVarianceSum = 0.0;
    for (double varianceValue : result.detectorVariance) {
        if (!isfinite(varianceValue) || varianceValue <= 0.0)
            return result;
        inverseVarianceSum += 1.0 / varianceValue;
    }
    for (size_t i = 0; i < result.value.size(); ++i) {
        result.value[i] =
            (1.0 / result.detectorVariance[i]) / inverseVarianceSum;
    }
    result.resolution = 1.0 / sqrt(inverseVarianceSum);
    result.valid = true;
    return result;
}

pair<double, double> RobustRange(vector<double> values) {
    values.erase(remove_if(values.begin(), values.end(),
                           [](double value) { return !isfinite(value); }),
                 values.end());
    if (values.empty()) return {0.0, 1.0};
    sort(values.begin(), values.end());
    const size_t lowIndex = static_cast<size_t>(
        0.005 * static_cast<double>(values.size() - 1));
    const size_t highIndex = static_cast<size_t>(
        0.995 * static_cast<double>(values.size() - 1));
    double low = values[lowIndex];
    double high = values[highIndex];
    if (!(high > low)) {
        low -= 0.5;
        high += 0.5;
    }
    return {low, high};
}

struct TimingEfficiencyResult {
    bool valid = false;
    size_t denominator = 0;
    size_t validTimes = 0;
    size_t bestCount = 0;
    double bestWindowStart = numeric_limits<double>::quiet_NaN();
    double bestWindowEnd = numeric_limits<double>::quiet_NaN();
    double bestWindowCenter = numeric_limits<double>::quiet_NaN();
    double bestEfficiency = numeric_limits<double>::quiet_NaN();
};

TimingEfficiencyResult WriteTimingEfficiencyCurve(
    const vector<double>& residuals, size_t denominator,
    double windowWidthNs,
    double stepNs, int histogramBins, TDirectory& directory,
    const string& name, const string& title,
    bool compactPlot = false) {
    TimingEfficiencyResult result;
    vector<double> sorted = residuals;
    sorted.erase(
        remove_if(sorted.begin(), sorted.end(),
                  [](double value) { return !isfinite(value); }),
        sorted.end());
    sort(sorted.begin(), sorted.end());
    result.denominator = denominator;
    result.validTimes = sorted.size();
    if (sorted.empty() || result.denominator == 0 ||
        sorted.size() > result.denominator ||
        !isfinite(windowWidthNs) || windowWidthNs <= 0.0 ||
        !isfinite(stepNs) || stepNs <= 0.0)
        return result;

    size_t right = 0;
    size_t bestLeft = 0;
    size_t bestRight = 0;
    for (size_t left = 0; left < sorted.size(); ++left) {
        if (right < left) right = left;
        while (right < sorted.size() &&
               sorted[right] - sorted[left] <= windowWidthNs)
            ++right;
        const size_t count = right - left;
        if (count > result.bestCount) {
            result.bestCount = count;
            bestLeft = left;
            bestRight = right;
        }
    }
    const double earliestStart =
        sorted[bestRight - 1] - windowWidthNs;
    const double latestStart = sorted[bestLeft];
    result.bestWindowStart =
        0.5 * (earliestStart + latestStart);
    result.bestWindowEnd =
        result.bestWindowStart + windowWidthNs;
    result.bestWindowCenter =
        result.bestWindowStart + 0.5 * windowWidthNs;
    result.bestEfficiency =
        static_cast<double>(result.bestCount) /
        static_cast<double>(result.denominator);
    result.valid = true;

    const auto [robustLow, robustHigh] = RobustRange(sorted);
    double scanLow = robustLow - windowWidthNs;
    double scanHigh = robustHigh;
    const double span = max(stepNs, scanHigh - scanLow);
    const double effectiveStep =
        max(stepNs, span / 20000.0);
    vector<double> windowStarts;
    for (double start = scanLow; start <= scanHigh;
         start += effectiveStep)
        windowStarts.push_back(start);
    windowStarts.push_back(result.bestWindowStart);
    sort(windowStarts.begin(), windowStarts.end());
    windowStarts.erase(
        unique(windowStarts.begin(), windowStarts.end(),
               [](double first, double second) {
                   return abs(first - second) < 1.0e-9;
               }),
        windowStarts.end());

    TDirectory::TContext context(&directory);
    TGraph curve(static_cast<int>(windowStarts.size()));
    curve.SetName(("g_" + name).c_str());
    ostringstream graphTitle;
    graphTitle << ";Window start t_{start} [ns];Efficiency";
    curve.SetTitle(graphTitle.str().c_str());
    for (size_t point = 0; point < windowStarts.size(); ++point) {
        const double start = windowStarts[point];
        const auto first =
            lower_bound(sorted.begin(), sorted.end(), start);
        const auto last =
            upper_bound(sorted.begin(), sorted.end(),
                        start + windowWidthNs);
        const double efficiency =
            static_cast<double>(last - first) /
            static_cast<double>(result.denominator);
        curve.SetPoint(
            static_cast<int>(point), start, efficiency);
    }
    curve.SetLineColor(kBlue + 1);
    curve.SetLineWidth(2);
    curve.SetMinimum(0.0);
    curve.SetMaximum(min(1.05, max(0.05, 1.08 * result.bestEfficiency)));

    TCanvas canvas(
        ("c_" + name).c_str(), title.c_str(), 1000, 700);
    canvas.SetTopMargin(0.20);
    curve.Draw("AL");
    TLatex heading;
    heading.SetNDC();
    heading.SetTextAlign(22);
    heading.SetTextSize(0.040);
    heading.DrawLatex(0.5, 0.96, title.c_str());
    TLine bestLine(
        result.bestWindowStart, 0.0,
        result.bestWindowStart, result.bestEfficiency);
    bestLine.SetLineColor(kRed + 1);
    bestLine.SetLineStyle(2);
    bestLine.SetLineWidth(2);
    bestLine.Draw("SAME");
    TMarker bestMarker(
        result.bestWindowStart, result.bestEfficiency, 29);
    bestMarker.SetMarkerColor(kRed + 1);
    bestMarker.SetMarkerSize(2.0);
    bestMarker.Draw("SAME");
    ostringstream maximumAnnotation;
    maximumAnnotation << fixed << setprecision(3)
                      << "max = "
                      << 100.0 * result.bestEfficiency
                      << "% (" << result.bestCount << '/'
                      << result.denominator << ')';
    ostringstream windowAnnotation;
    windowAnnotation << fixed << setprecision(3)
                     << "best " << windowWidthNs
                     << " ns window = ["
                     << result.bestWindowStart << ", "
                     << result.bestWindowEnd << "] ns";
    TLatex label;
    label.SetNDC();
    label.SetTextSize(0.030);
    label.SetTextColor(kRed + 1);
    label.DrawLatex(
        0.14, 0.90, maximumAnnotation.str().c_str());
    label.DrawLatex(
        0.14, 0.85, windowAnnotation.str().c_str());
    canvas.Modified();
    canvas.Update();
    canvas.Write();

    const auto [residualLow, residualHigh] = RobustRange(sorted);
    const double plotLow = compactPlot
                               ? min(residualLow,
                                     result.bestWindowStart)
                               : min(residualLow,
                                     result.bestWindowStart) -
                                     0.05 * windowWidthNs;
    const double plotHigh = compactPlot
                                ? max(residualHigh,
                                      result.bestWindowEnd)
                                : max(residualHigh,
                                      result.bestWindowEnd) +
                                      0.05 * windowWidthNs;
    TH1D residualHistogram(
        ("h_residual_" + name).c_str(),
        ";#Deltat=t-T0_{scope} [ns];Events",
        histogramBins, plotLow, plotHigh);
    residualHistogram.SetDirectory(nullptr);
    for (double value : sorted) residualHistogram.Fill(value);
    residualHistogram.SetLineColor(kBlue + 1);
    residualHistogram.SetLineWidth(2);
    residualHistogram.SetStats(false);

    vector<double> fitValues;
    fitValues.reserve(result.bestCount);
    copy_if(
        sorted.begin(), sorted.end(), back_inserter(fitValues),
        [&](double value) {
            return value >= result.bestWindowStart &&
                   value <= result.bestWindowEnd;
        });
    const FitResult fit = FitGaussianWithoutWriting(fitValues);
    unique_ptr<TF1> gaussian;
    if (isfinite(fit.mean) && isfinite(fit.sigma) &&
        fit.sigma > 0.0) {
        const double fitLow =
            max(plotLow, fit.mean - 2.5 * fit.sigma);
        const double fitHigh =
            min(plotHigh, fit.mean + 2.5 * fit.sigma);
        gaussian = make_unique<TF1>(
            ("fit_residual_" + name).c_str(), "gaus",
            fitLow, fitHigh);
        gaussian->SetParameters(
            residualHistogram.GetMaximum(), fit.mean, fit.sigma);
        residualHistogram.Fit(
            gaussian.get(), "RQ0N", "", fitLow, fitHigh);
        gaussian->SetLineColor(kGreen + 2);
        gaussian->SetLineWidth(3);
    }

    TCanvas residualCanvas(
        ("c_residual_" + name).c_str(),
        (title + " residual and best window").c_str(), 1000, 700);
    residualCanvas.SetTopMargin(0.24);
    const double yMaximum = 1.08 * residualHistogram.GetMaximum();
    residualHistogram.SetMaximum(yMaximum);
    residualHistogram.Draw("HIST");
    heading.DrawLatex(0.5, 0.96, title.c_str());
    TLine windowStart(
        result.bestWindowStart, 0.0,
        result.bestWindowStart, yMaximum);
    TLine windowEnd(
        result.bestWindowEnd, 0.0,
        result.bestWindowEnd, yMaximum);
    windowStart.SetLineColor(kRed + 1);
    windowEnd.SetLineColor(kRed + 1);
    windowStart.SetLineStyle(2);
    windowEnd.SetLineStyle(2);
    windowStart.SetLineWidth(2);
    windowEnd.SetLineWidth(2);
    windowStart.Draw("SAME");
    windowEnd.Draw("SAME");
    if (gaussian) gaussian->Draw("SAME");
    label.SetTextColor(kRed + 1);
    label.DrawLatex(
        0.14, 0.90, maximumAnnotation.str().c_str());
    label.DrawLatex(
        0.14, 0.85, windowAnnotation.str().c_str());
    if (gaussian) {
        ostringstream fitAnnotation;
        fitAnnotation << fixed << setprecision(3)
                      << "Gaussian: #mu = "
                      << gaussian->GetParameter(1)
                      << " ns, #sigma = "
                      << abs(gaussian->GetParameter(2))
                      << " ns";
        label.SetTextColor(kGreen + 2);
        label.DrawLatex(
            0.14, 0.80, fitAnnotation.str().c_str());
    }
    residualCanvas.Modified();
    residualCanvas.Update();
    residualCanvas.Write();
    return result;
}

double PearsonCorrelation(const vector<double>& x, const vector<double>& y) {
    double sumX = 0.0, sumY = 0.0;
    size_t entries = 0;
    for (size_t i = 0; i < min(x.size(), y.size()); ++i) {
        if (!isfinite(x[i]) || !isfinite(y[i])) continue;
        sumX += x[i];
        sumY += y[i];
        ++entries;
    }
    if (entries < 2) return numeric_limits<double>::quiet_NaN();
    const double meanX = sumX / entries;
    const double meanY = sumY / entries;
    double covariance = 0.0, varianceX = 0.0, varianceY = 0.0;
    for (size_t i = 0; i < min(x.size(), y.size()); ++i) {
        if (!isfinite(x[i]) || !isfinite(y[i])) continue;
        const double dx = x[i] - meanX;
        const double dy = y[i] - meanY;
        covariance += dx * dy;
        varianceX += dx * dx;
        varianceY += dy * dy;
    }
    if (varianceX <= 0.0 || varianceY <= 0.0)
        return numeric_limits<double>::quiet_NaN();
    return covariance / sqrt(varianceX * varianceY);
}

double WriteDUTCorrelation(
    const vector<double>& factor, const vector<double>& residual,
    TDirectory& directory,
    const string& name, const string& axisTitle, int bins) {
    const auto [xLow, xHigh] = RobustRange(factor);
    const auto [yLow, yHigh] = RobustRange(residual);
    TDirectory::TContext context(&directory);
    TH2D histogram(
        ("residual_vs_" + name).c_str(),
        ("DUT-track residual vs " + axisTitle + ";" + axisTitle +
         ";t_{DUT}-t_{track} [ns]")
            .c_str(),
        bins, xLow, xHigh, bins, yLow, yHigh);
    TProfile profile(
        ("profile_vs_" + name).c_str(),
        ("Residual profile vs " + axisTitle + ";" + axisTitle +
         ";<#Deltat> [ns]")
            .c_str(),
        bins, xLow, xHigh);
    for (size_t i = 0; i < factor.size(); ++i) {
        if (i < residual.size() && isfinite(factor[i]) &&
            isfinite(residual[i])) {
            histogram.Fill(factor[i], residual[i]);
            profile.Fill(factor[i], residual[i]);
        }
    }
    histogram.Write();
    profile.Write();
    return PearsonCorrelation(factor, residual);
}

bool WriteTimingOutput(const string& outputPath,
                       const map<TrackKey, map<int, DetectorTimes>>& trackTimes,
                       const map<TrackKey, int>& eventIDs,
                       const OscilloscopeT0Result& oscilloscopeT0,
                       const array<int, 3>& trackerIDs,
                       const DUTTimingResult& dutTiming,
                       const TrackTimeWeights& trackWeights,
                       double timingEfficiencyWindowNs,
                       double timingEfficiencyStepNs,
                       int histogramBins) {
    unique_ptr<TFile> output(TFile::Open(outputPath.c_str(), "RECREATE"));
    if (!output || output->IsZombie()) {
        cerr << "[TimeResolution] cannot create " << outputPath << '\n';
        return false;
    }
    TTree eventTree(
        "EventTimes", "Tracker detector and matched oscilloscope T0 times");
    ULong64_t outRawEventID = 0;
    Int_t outTrackIndex = 0, outDetID = 0;
    Bool_t outHasOscilloscopeT0 = false;
    Double_t outOscilloscopeT0 =
        numeric_limits<double>::quiet_NaN();
    Double_t outOscilloscopeT0Resolution =
        numeric_limits<double>::quiet_NaN();
    Int_t outOscilloscopeValidChannels = 0;
    Int_t outOscilloscopeTraceFileIndex = -1;
    Int_t outOscilloscopeSegmentIndex = -1;
    array<Double_t, kEstimatorCount> estimatorValues{};
    eventTree.Branch("rawEventID", &outRawEventID);
    eventTree.Branch("trackIndex", &outTrackIndex);
    eventTree.Branch("detectorID", &outDetID);
    eventTree.Branch("hasOscilloscopeT0", &outHasOscilloscopeT0);
    eventTree.Branch("oscilloscopeT0", &outOscilloscopeT0);
    eventTree.Branch(
        "oscilloscopeT0Resolution", &outOscilloscopeT0Resolution);
    eventTree.Branch(
        "oscilloscopeT0ValidChannels", &outOscilloscopeValidChannels);
    eventTree.Branch(
        "oscilloscopeTraceFileIndex", &outOscilloscopeTraceFileIndex);
    eventTree.Branch(
        "oscilloscopeSegmentIndex", &outOscilloscopeSegmentIndex);
    for (size_t estimator = 0; estimator < kEstimatorCount; ++estimator)
        eventTree.Branch(kEstimatorNames[estimator].c_str(), &estimatorValues[estimator]);

    TTree trackTimeTree("TrackTimes", "Inverse-variance weighted tracker time");
    Int_t trackEventID = 0;
    ULong64_t trackRawEventID = 0;
    Int_t trackIndex = 0;
    Int_t tracker1ID = trackerIDs[0], tracker2ID = trackerIDs[1],
          tracker3ID = trackerIDs[2];
    Double_t tracker1Time = numeric_limits<double>::quiet_NaN();
    Double_t tracker2Time = numeric_limits<double>::quiet_NaN();
    Double_t tracker3Time = numeric_limits<double>::quiet_NaN();
    Double_t trackTime = numeric_limits<double>::quiet_NaN();
    Bool_t trackHasOscilloscopeT0 = false;
    Double_t trackOscilloscopeT0 =
        numeric_limits<double>::quiet_NaN();
    Double_t trackOscilloscopeT0Resolution =
        numeric_limits<double>::quiet_NaN();
    Double_t trackMinusOscilloscopeT0 =
        numeric_limits<double>::quiet_NaN();
    Double_t tracker1Weight = trackWeights.value[0];
    Double_t tracker2Weight = trackWeights.value[1];
    Double_t tracker3Weight = trackWeights.value[2];
    trackTimeTree.Branch("eventID", &trackEventID);
    trackTimeTree.Branch("rawEventID", &trackRawEventID);
    trackTimeTree.Branch("trackIndex", &trackIndex);
    trackTimeTree.Branch("tracker1ID", &tracker1ID);
    trackTimeTree.Branch("tracker2ID", &tracker2ID);
    trackTimeTree.Branch("tracker3ID", &tracker3ID);
    trackTimeTree.Branch("tracker1Time", &tracker1Time);
    trackTimeTree.Branch("tracker2Time", &tracker2Time);
    trackTimeTree.Branch("tracker3Time", &tracker3Time);
    trackTimeTree.Branch("trackTime", &trackTime);
    trackTimeTree.Branch(
        "hasOscilloscopeT0", &trackHasOscilloscopeT0);
    trackTimeTree.Branch("oscilloscopeT0", &trackOscilloscopeT0);
    trackTimeTree.Branch(
        "oscilloscopeT0Resolution", &trackOscilloscopeT0Resolution);
    trackTimeTree.Branch(
        "trackMinusOscilloscopeT0", &trackMinusOscilloscopeT0);
    trackTimeTree.Branch("tracker1Weight", &tracker1Weight);
    trackTimeTree.Branch("tracker2Weight", &tracker2Weight);
    trackTimeTree.Branch("tracker3Weight", &tracker3Weight);

    map<pair<int, size_t>, vector<double>> detectorTimeSamples;
    map<tuple<int, int, size_t>, vector<double>> pairResiduals;
    map<pair<int, size_t>, T0ResidualSamples> detectorT0Residuals;
    vector<double> trackTimeSamples;
    T0ResidualSamples trackT0Residuals;
    array<vector<double>, 3> threeTrackerPairResiduals;
    size_t tracksWithOscilloscopeT0 = 0;
    size_t completeThreeTrackerTimes = 0;
    for (const auto& [key, detectors] : trackTimes) {
        const auto scopeReference =
            oscilloscopeT0.references.find(key.first);
        const bool hasOscilloscopeT0 =
            scopeReference != oscilloscopeT0.references.end();
        if (hasOscilloscopeT0) ++tracksWithOscilloscopeT0;
        const OscilloscopeT0Reference* t0 =
            hasOscilloscopeT0 ? &scopeReference->second : nullptr;
        for (const auto& [detectorID, times] : detectors) {
            outRawEventID = key.first;
            outTrackIndex = key.second;
            outDetID = detectorID;
            outHasOscilloscopeT0 = hasOscilloscopeT0;
            outOscilloscopeT0 =
                t0 ? t0->time : numeric_limits<double>::quiet_NaN();
            outOscilloscopeT0Resolution =
                t0 ? t0->resolution
                   : numeric_limits<double>::quiet_NaN();
            outOscilloscopeValidChannels =
                t0 ? t0->validChannels : 0;
            outOscilloscopeTraceFileIndex =
                t0 ? t0->traceFileIndex : -1;
            outOscilloscopeSegmentIndex =
                t0 ? t0->segmentIndex : -1;
            estimatorValues = times.value;
            eventTree.Fill();
            for (size_t estimator = 0; estimator < kEstimatorCount; ++estimator) {
                if (!isfinite(times.value[estimator])) continue;
                detectorTimeSamples[{detectorID, estimator}].push_back(times.value[estimator]);
                if (t0)
                    detectorT0Residuals[{detectorID, estimator}].Add(
                        times.value[estimator] - t0->time,
                        t0->resolution);
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

        array<double, 3> selectedTimes{};
        bool hasAllTrackerTimes = true;
        for (size_t i = 0; i < trackerIDs.size(); ++i) {
            const auto detector = detectors.find(trackerIDs[i]);
            if (detector == detectors.end() ||
                !isfinite(detector->second.value[kTrackerTimingEstimator])) {
                hasAllTrackerTimes = false;
                break;
            }
            selectedTimes[i] =
                detector->second.value[kTrackerTimingEstimator];
        }
        if (hasAllTrackerTimes) {
            const auto event = eventIDs.find(key);
            trackEventID = event == eventIDs.end() ? -1 : event->second;
            trackRawEventID = key.first;
            trackIndex = key.second;
            tracker1Time = selectedTimes[0];
            tracker2Time = selectedTimes[1];
            tracker3Time = selectedTimes[2];
            trackTime =
                tracker1Weight * tracker1Time +
                tracker2Weight * tracker2Time +
                tracker3Weight * tracker3Time;
            trackHasOscilloscopeT0 = hasOscilloscopeT0;
            trackOscilloscopeT0 =
                t0 ? t0->time : numeric_limits<double>::quiet_NaN();
            trackOscilloscopeT0Resolution =
                t0 ? t0->resolution
                   : numeric_limits<double>::quiet_NaN();
            trackMinusOscilloscopeT0 =
                t0 ? trackTime - t0->time
                   : numeric_limits<double>::quiet_NaN();
            trackTimeTree.Fill();
            trackTimeSamples.push_back(trackTime);
            if (t0)
                trackT0Residuals.Add(
                    trackMinusOscilloscopeT0, t0->resolution);
            threeTrackerPairResiduals[0].push_back(tracker1Time - tracker2Time);
            threeTrackerPairResiduals[1].push_back(tracker1Time - tracker3Time);
            threeTrackerPairResiduals[2].push_back(tracker2Time - tracker3Time);
            ++completeThreeTrackerTimes;
        }
    }
    eventTree.Write();
    trackTimeTree.Write();

    TTree dutTimeTree("DUTTimes", "DUT time relative to weighted tracker time");
    Int_t dutEventID = 0, dutTrackIndex = -1, dutDetectorID = -1;
    ULong64_t dutRawEventID = 0;
    Double_t dutTime = numeric_limits<double>::quiet_NaN();
    Double_t dutTrackTime = numeric_limits<double>::quiet_NaN();
    Double_t dutMinusTrackTime = numeric_limits<double>::quiet_NaN();
    Bool_t dutHasOscilloscopeT0 = false;
    Double_t dutOscilloscopeT0 =
        numeric_limits<double>::quiet_NaN();
    Double_t dutOscilloscopeT0Resolution =
        numeric_limits<double>::quiet_NaN();
    Double_t dutMinusOscilloscopeT0 =
        numeric_limits<double>::quiet_NaN();
    Double_t dutAmplitude = numeric_limits<double>::quiet_NaN();
    Double_t dutMeanAmplitude = numeric_limits<double>::quiet_NaN();
    Double_t dutClusterCharge = numeric_limits<double>::quiet_NaN();
    Double_t dutClusterMaxAmplitude = numeric_limits<double>::quiet_NaN();
    Double_t dutClusterSize = numeric_limits<double>::quiet_NaN();
    Double_t dutClusterCentroid = numeric_limits<double>::quiet_NaN();
    Double_t dutClusterLocalX = numeric_limits<double>::quiet_NaN();
    Double_t dutClusterLocalY = numeric_limits<double>::quiet_NaN();
    Double_t dutPredictedX = numeric_limits<double>::quiet_NaN();
    Double_t dutPredictedY = numeric_limits<double>::quiet_NaN();
    dutTimeTree.Branch("eventID", &dutEventID);
    dutTimeTree.Branch("rawEventID", &dutRawEventID);
    dutTimeTree.Branch("trackIndex", &dutTrackIndex);
    dutTimeTree.Branch("dutID", &dutDetectorID);
    dutTimeTree.Branch("dutTime", &dutTime);
    dutTimeTree.Branch("trackTime", &dutTrackTime);
    dutTimeTree.Branch("dutMinusTrackTime", &dutMinusTrackTime);
    dutTimeTree.Branch(
        "hasOscilloscopeT0", &dutHasOscilloscopeT0);
    dutTimeTree.Branch("oscilloscopeT0", &dutOscilloscopeT0);
    dutTimeTree.Branch(
        "oscilloscopeT0Resolution", &dutOscilloscopeT0Resolution);
    dutTimeTree.Branch(
        "dutMinusOscilloscopeT0", &dutMinusOscilloscopeT0);
    dutTimeTree.Branch("amplitude", &dutAmplitude);
    dutTimeTree.Branch("meanAmplitude", &dutMeanAmplitude);
    dutTimeTree.Branch("clusterCharge", &dutClusterCharge);
    dutTimeTree.Branch("clusterMaxAmplitude", &dutClusterMaxAmplitude);
    dutTimeTree.Branch("clusterSize", &dutClusterSize);
    dutTimeTree.Branch("clusterCentroid", &dutClusterCentroid);
    dutTimeTree.Branch("clusterLocalX", &dutClusterLocalX);
    dutTimeTree.Branch("clusterLocalY", &dutClusterLocalY);
    dutTimeTree.Branch("predictedX", &dutPredictedX);
    dutTimeTree.Branch("predictedY", &dutPredictedY);
    map<int, T0ResidualSamples> dutT0Residuals;
    map<int, T0ResidualSamples> dutEfficiencyT0Residuals;
    for (const auto& [detectorID, samples] : dutTiming.samplesByDetector) {
        for (const DUTTimingSample& sample : samples) {
            dutEventID = sample.eventID;
            dutRawEventID = sample.rawEventID;
            dutTrackIndex = sample.trackIndex;
            dutDetectorID = detectorID;
            dutTime = sample.dutTime;
            dutTrackTime = sample.trackTime;
            dutMinusTrackTime = sample.residual;
            const auto scopeReference =
                oscilloscopeT0.references.find(sample.rawEventID);
            dutHasOscilloscopeT0 =
                scopeReference != oscilloscopeT0.references.end();
            if (dutHasOscilloscopeT0) {
                dutOscilloscopeT0 = scopeReference->second.time;
                dutOscilloscopeT0Resolution =
                    scopeReference->second.resolution;
                dutMinusOscilloscopeT0 =
                    sample.dutTime - dutOscilloscopeT0;
                dutT0Residuals[detectorID].Add(
                    dutMinusOscilloscopeT0,
                    dutOscilloscopeT0Resolution);
                if (sample.insideActiveArea)
                    dutEfficiencyT0Residuals[detectorID].Add(
                        dutMinusOscilloscopeT0,
                        dutOscilloscopeT0Resolution);
            } else {
                dutOscilloscopeT0 =
                    numeric_limits<double>::quiet_NaN();
                dutOscilloscopeT0Resolution =
                    numeric_limits<double>::quiet_NaN();
                dutMinusOscilloscopeT0 =
                    numeric_limits<double>::quiet_NaN();
            }
            dutAmplitude = sample.amplitude;
            dutMeanAmplitude = sample.meanAmplitude;
            dutClusterCharge = sample.clusterCharge;
            dutClusterMaxAmplitude = sample.clusterMaxAmplitude;
            dutClusterSize = sample.clusterSize;
            dutClusterCentroid = sample.clusterCentroid;
            dutClusterLocalX = sample.clusterLocalX;
            dutClusterLocalY = sample.clusterLocalY;
            dutPredictedX = sample.predictedX;
            dutPredictedY = sample.predictedY;
            dutTimeTree.Fill();
        }
    }
    dutTimeTree.Write();

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

    TDirectory* threeTrackerDirectory = output->mkdir("ThreeTracker");
    WriteAndFit(
        trackTimeSamples, *threeTrackerDirectory, "track_time_weighted",
        "Inverse-variance weighted tracker time;t_{track} [ns];Tracks",
        histogramBins);
    const array<pair<size_t, size_t>, 3> pairIndices = {
        pair<size_t, size_t>{0, 1}, {0, 2}, {1, 2}};
    array<FitResult, 3> threeTrackerPairFits;
    for (size_t pairIndex = 0; pairIndex < pairIndices.size(); ++pairIndex) {
        const auto [first, second] = pairIndices[pairIndex];
        const string pairLabel =
            "d" + to_string(trackerIDs[first]) + "_d" + to_string(trackerIDs[second]);
        threeTrackerPairFits[pairIndex] = WriteAndFit(
            threeTrackerPairResiduals[pairIndex], *threeTrackerDirectory,
            "three_tracker_pair_" + pairLabel,
            "Three-tracker common-sample residual; t_{" +
                to_string(trackerIDs[first]) + "} - t_{" +
                to_string(trackerIDs[second]) + "} [ns];Tracks",
            histogramBins);
    }
    const FitResult* fit12 = &threeTrackerPairFits[0];
    const FitResult* fit13 = &threeTrackerPairFits[1];
    const FitResult* fit23 = &threeTrackerPairFits[2];

    TTree threeTrackerSummary(
        "ThreeTrackerResolution",
        "Individual tracker resolutions and inverse-variance weighted track-time resolution");
    Int_t resolutionValid = 0;
    Double_t sigma12Ns = numeric_limits<double>::quiet_NaN();
    Double_t sigma13Ns = numeric_limits<double>::quiet_NaN();
    Double_t sigma23Ns = numeric_limits<double>::quiet_NaN();
    Double_t tracker1VarianceNs2 = numeric_limits<double>::quiet_NaN();
    Double_t tracker2VarianceNs2 = numeric_limits<double>::quiet_NaN();
    Double_t tracker3VarianceNs2 = numeric_limits<double>::quiet_NaN();
    Double_t tracker1ResolutionNs = numeric_limits<double>::quiet_NaN();
    Double_t tracker2ResolutionNs = numeric_limits<double>::quiet_NaN();
    Double_t tracker3ResolutionNs = numeric_limits<double>::quiet_NaN();
    Double_t tracker1ResolutionErrorNs = numeric_limits<double>::quiet_NaN();
    Double_t tracker2ResolutionErrorNs = numeric_limits<double>::quiet_NaN();
    Double_t tracker3ResolutionErrorNs = numeric_limits<double>::quiet_NaN();
    Double_t trackResolutionNs = numeric_limits<double>::quiet_NaN();
    Double_t trackResolutionErrorNs = numeric_limits<double>::quiet_NaN();
    ULong64_t completeTracks = static_cast<ULong64_t>(completeThreeTrackerTimes);
    threeTrackerSummary.Branch("valid", &resolutionValid);
    threeTrackerSummary.Branch("tracker1ID", &tracker1ID);
    threeTrackerSummary.Branch("tracker2ID", &tracker2ID);
    threeTrackerSummary.Branch("tracker3ID", &tracker3ID);
    threeTrackerSummary.Branch("completeTracks", &completeTracks);
    threeTrackerSummary.Branch("sigma12Ns", &sigma12Ns);
    threeTrackerSummary.Branch("sigma13Ns", &sigma13Ns);
    threeTrackerSummary.Branch("sigma23Ns", &sigma23Ns);
    threeTrackerSummary.Branch("tracker1VarianceNs2", &tracker1VarianceNs2);
    threeTrackerSummary.Branch("tracker2VarianceNs2", &tracker2VarianceNs2);
    threeTrackerSummary.Branch("tracker3VarianceNs2", &tracker3VarianceNs2);
    threeTrackerSummary.Branch("tracker1ResolutionNs", &tracker1ResolutionNs);
    threeTrackerSummary.Branch("tracker2ResolutionNs", &tracker2ResolutionNs);
    threeTrackerSummary.Branch("tracker3ResolutionNs", &tracker3ResolutionNs);
    threeTrackerSummary.Branch("tracker1ResolutionErrorNs", &tracker1ResolutionErrorNs);
    threeTrackerSummary.Branch("tracker2ResolutionErrorNs", &tracker2ResolutionErrorNs);
    threeTrackerSummary.Branch("tracker3ResolutionErrorNs", &tracker3ResolutionErrorNs);
    threeTrackerSummary.Branch("tracker1Weight", &tracker1Weight);
    threeTrackerSummary.Branch("tracker2Weight", &tracker2Weight);
    threeTrackerSummary.Branch("tracker3Weight", &tracker3Weight);
    threeTrackerSummary.Branch("trackResolutionNs", &trackResolutionNs);
    threeTrackerSummary.Branch("trackResolutionErrorNs", &trackResolutionErrorNs);

    const bool hasPairFits =
        fit12 && fit13 && fit23 && isfinite(fit12->sigma) &&
        isfinite(fit13->sigma) && isfinite(fit23->sigma);
    if (hasPairFits) {
        sigma12Ns = fit12->sigma;
        sigma13Ns = fit13->sigma;
        sigma23Ns = fit23->sigma;
        tracker1VarianceNs2 =
            0.5 * (sigma12Ns * sigma12Ns + sigma13Ns * sigma13Ns -
                   sigma23Ns * sigma23Ns);
        tracker2VarianceNs2 =
            0.5 * (sigma12Ns * sigma12Ns + sigma23Ns * sigma23Ns -
                   sigma13Ns * sigma13Ns);
        tracker3VarianceNs2 =
            0.5 * (sigma13Ns * sigma13Ns + sigma23Ns * sigma23Ns -
                   sigma12Ns * sigma12Ns);
        resolutionValid = trackWeights.valid &&
                          tracker1VarianceNs2 > 0.0 &&
                          tracker2VarianceNs2 > 0.0 &&
                          tracker3VarianceNs2 > 0.0;
        if (resolutionValid) {
            tracker1ResolutionNs = sqrt(tracker1VarianceNs2);
            tracker2ResolutionNs = sqrt(tracker2VarianceNs2);
            tracker3ResolutionNs = sqrt(tracker3VarianceNs2);
            trackResolutionNs = trackWeights.resolution;

            auto resolutionError = [](double resolution,
                                      const FitResult& first,
                                      const FitResult& second,
                                      const FitResult& third) {
                if (resolution <= 0.0 || !isfinite(first.sigmaError) ||
                    !isfinite(second.sigmaError) || !isfinite(third.sigmaError))
                    return numeric_limits<double>::quiet_NaN();
                const double varianceError = sqrt(
                    pow(first.sigma * first.sigmaError, 2) +
                    pow(second.sigma * second.sigmaError, 2) +
                    pow(third.sigma * third.sigmaError, 2));
                return varianceError / (2.0 * resolution);
            };
            tracker1ResolutionErrorNs =
                resolutionError(tracker1ResolutionNs, *fit12, *fit13, *fit23);
            tracker2ResolutionErrorNs =
                resolutionError(tracker2ResolutionNs, *fit12, *fit23, *fit13);
            tracker3ResolutionErrorNs =
                resolutionError(tracker3ResolutionNs, *fit13, *fit23, *fit12);
            if (trackResolutionNs > 0.0 &&
                isfinite(tracker1ResolutionErrorNs) &&
                isfinite(tracker2ResolutionErrorNs) &&
                isfinite(tracker3ResolutionErrorNs)) {
                trackResolutionErrorNs =
                    sqrt(pow(tracker1Weight * tracker1Weight *
                                 tracker1ResolutionNs *
                                 tracker1ResolutionErrorNs,
                             2) +
                         pow(tracker2Weight * tracker2Weight *
                                 tracker2ResolutionNs *
                                 tracker2ResolutionErrorNs,
                             2) +
                         pow(tracker3Weight * tracker3Weight *
                                 tracker3ResolutionNs *
                                 tracker3ResolutionErrorNs,
                             2)) /
                    trackResolutionNs;
            }

            const array<double, 3> resolutions = {
                tracker1ResolutionNs, tracker2ResolutionNs, tracker3ResolutionNs};
            const array<double, 3> errors = {
                tracker1ResolutionErrorNs, tracker2ResolutionErrorNs,
                tracker3ResolutionErrorNs};
            for (size_t i = 0; i < trackerIDs.size(); ++i) {
                FitResult derived;
                derived.entries = min({fit12->entries, fit13->entries, fit23->entries});
                derived.sigma = resolutions[i];
                derived.sigmaError = errors[i];
                fillSummary("three_tracker_detector", kEstimatorNames[kTrackerTimingEstimator],
                            trackerIDs[i], -1, derived, resolutions[i]);
            }
            FitResult trackDerived;
            trackDerived.entries = static_cast<long long>(completeThreeTrackerTimes);
            trackDerived.sigma = trackResolutionNs;
            trackDerived.sigmaError = trackResolutionErrorNs;
            fillSummary("track_time_weighted", kEstimatorNames[kTrackerTimingEstimator],
                        -1, -1, trackDerived, trackResolutionNs);
        }
    }
    threeTrackerDirectory->cd();
    threeTrackerSummary.Fill();
    threeTrackerSummary.Write();

    ostringstream trackerSummary;
    trackerSummary << fixed << setprecision(1)
                   << "track reference ";
    if (resolutionValid) {
        trackerSummary << trackResolutionNs << " ns · trackers "
                       << tracker1ResolutionNs << '/'
                       << tracker2ResolutionNs << '/'
                       << tracker3ResolutionNs << " ns · "
                       << Terminal::Count(completeThreeTrackerTimes)
                       << " tracks";
    } else {
        trackerSummary << "invalid";
    }
    Terminal::Detail(trackerSummary.str());
    if (Terminal::Verbose() && resolutionValid) {
        ostringstream weights;
        weights << fixed << setprecision(2)
                << "tracker IDs " << trackerIDs[0] << '/'
                << trackerIDs[1] << '/' << trackerIDs[2]
                << " · weights " << tracker1Weight << '/'
                << tracker2Weight << '/' << tracker3Weight;
        Terminal::Detail(Terminal::Muted(weights.str()));
    }

    if (!oscilloscopeT0.references.empty()) {
        TDirectory* t0Directory = output->mkdir("OscilloscopeT0");
        TDirectory* detectorDirectory =
            t0Directory->mkdir("DetectorResiduals");
        auto subtractT0Resolution =
            [](const FitResult& fit,
               const T0ResidualSamples& samples) {
                const double referenceVariance =
                    samples.MeanReferenceVariance();
                if (!isfinite(fit.sigma) ||
                    !isfinite(referenceVariance))
                    return numeric_limits<double>::quiet_NaN();
                const double intrinsicVariance =
                    fit.sigma * fit.sigma - referenceVariance;
                return intrinsicVariance >= 0.0
                           ? sqrt(intrinsicVariance)
                           : numeric_limits<double>::quiet_NaN();
            };

        for (const auto& [key, samples] : detectorT0Residuals) {
            const int detectorID = key.first;
            const size_t estimator = key.second;
            const string estimatorLabel = kEstimatorNames[estimator];
            const string name =
                SafeName("detector_minus_t0_d" +
                         to_string(detectorID) + "_" + estimatorLabel);
            const FitResult fit = WriteAndFit(
                samples.residual, *detectorDirectory, name,
                "Detector time relative to oscilloscope T0;"
                "t_{detector}-T0_{scope} [ns];Entries",
                histogramBins);
            fillSummary(
                "oscilloscope_t0_detector", estimatorLabel,
                detectorID, -1, fit,
                subtractT0Resolution(fit, samples));
        }

        const FitResult trackT0Fit = WriteAndFit(
            trackT0Residuals.residual, *t0Directory,
            "track_minus_oscilloscope_t0",
            "Weighted track time relative to oscilloscope T0;"
            "t_{track}-T0_{scope} [ns];Tracks",
            histogramBins);
        const double trackT0Resolution =
            subtractT0Resolution(trackT0Fit, trackT0Residuals);
        fillSummary(
            "oscilloscope_t0_track",
            kEstimatorNames[kTrackerTimingEstimator],
            -1, -1, trackT0Fit, trackT0Resolution);

        TTree t0Summary(
            "OscilloscopeT0Summary",
            "Oscilloscope channel alignment, weighting and matched T0 summary");
        Double_t channelResolution[3] = {
            oscilloscopeT0.channelResolution[0],
            oscilloscopeT0.channelResolution[1],
            oscilloscopeT0.channelResolution[2]};
        Double_t channelResolutionError[3] = {
            oscilloscopeT0.channelResolutionError[0],
            oscilloscopeT0.channelResolutionError[1],
            oscilloscopeT0.channelResolutionError[2]};
        Double_t channelOffset[3] = {
            oscilloscopeT0.channelOffset[0],
            oscilloscopeT0.channelOffset[1],
            oscilloscopeT0.channelOffset[2]};
        Double_t channelWeight[3] = {
            oscilloscopeT0.channelWeight[0],
            oscilloscopeT0.channelWeight[1],
            oscilloscopeT0.channelWeight[2]};
        Double_t combinedResolutionNs = oscilloscopeT0.resolution;
        Double_t measuredTrackMinusT0SigmaNs = trackT0Fit.sigma;
        Double_t intrinsicTrackResolutionNs = trackT0Resolution;
        ULong64_t inputEntries = oscilloscopeT0.inputEntries;
        ULong64_t matchedRawEventIDs =
            oscilloscopeT0.references.size();
        ULong64_t matchedTracks = tracksWithOscilloscopeT0;
        ULong64_t invalidEventIDs =
            oscilloscopeT0.invalidEventIDs;
        ULong64_t duplicateEventIDs =
            oscilloscopeT0.duplicateEventIDs;
        t0Summary.Branch(
            "channelResolutionNs", channelResolution,
            "channelResolutionNs[3]/D");
        t0Summary.Branch(
            "channelResolutionErrorNs", channelResolutionError,
            "channelResolutionErrorNs[3]/D");
        t0Summary.Branch(
            "channelOffsetNs", channelOffset,
            "channelOffsetNs[3]/D");
        t0Summary.Branch(
            "channelWeight", channelWeight, "channelWeight[3]/D");
        t0Summary.Branch(
            "combinedResolutionNs", &combinedResolutionNs);
        t0Summary.Branch(
            "measuredTrackMinusT0SigmaNs",
            &measuredTrackMinusT0SigmaNs);
        t0Summary.Branch(
            "intrinsicTrackResolutionNs",
            &intrinsicTrackResolutionNs);
        t0Summary.Branch("inputEntries", &inputEntries);
        t0Summary.Branch("matchedRawEventIDs", &matchedRawEventIDs);
        t0Summary.Branch("matchedTracks", &matchedTracks);
        t0Summary.Branch("invalidEventIDs", &invalidEventIDs);
        t0Summary.Branch("duplicateEventIDs", &duplicateEventIDs);
        t0Directory->cd();
        t0Summary.Fill();
        t0Summary.Write();

        struct EfficiencyRow {
            string source;
            TimingEfficiencyResult result;
        };
        vector<EfficiencyRow> efficiencyRows;
        TDirectory* efficiencyDirectory =
            t0Directory->mkdir("TimingEfficiency");
        map<int, size_t> dutEfficiencyDenominators;
        map<int, size_t> dutMatchedHitCounts;
        map<int, TimingEfficiencyResult> dutEfficiencyResults;
        auto countCasesWithT0 =
            [&](const set<TrackKey>& cases) {
                return static_cast<size_t>(count_if(
                    cases.begin(), cases.end(),
                    [&](const TrackKey& key) {
                        return oscilloscopeT0.references.count(
                                   key.first) != 0;
                    }));
            };
        for (const auto& [dutID, cases] :
             dutTiming.activeAreaTrackCases) {
            dutEfficiencyDenominators[dutID] =
                countCasesWithT0(cases);
        }
        for (const auto& [dutID, cases] :
             dutTiming.activeAreaMatchedHitCases) {
            dutMatchedHitCounts[dutID] =
                countCasesWithT0(cases);
        }

        for (const auto& [dutID, samples] :
             dutEfficiencyT0Residuals) {
            const string source = "DUT " + to_string(dutID);
            const size_t denominator =
                dutEfficiencyDenominators[dutID];
            const TimingEfficiencyResult efficiency =
                WriteTimingEfficiencyCurve(
                    samples.residual, denominator,
                    timingEfficiencyWindowNs, timingEfficiencyStepNs,
                    histogramBins,
                    *efficiencyDirectory,
                    SafeName("dut_d" + to_string(dutID)),
                    source + " timing-window efficiency", true);
            if (efficiency.valid)
                efficiencyRows.push_back({source, efficiency});
            if (efficiency.valid)
                dutEfficiencyResults[dutID] = efficiency;
        }

        if (!efficiencyRows.empty()) {
            cout << "\n[TimeResolution] " << timingEfficiencyWindowNs
                 << " ns timing-window efficiency relative to "
                    "oscilloscope T0 (DUT only)\n"
                 << "  " << left << setw(20) << "Source"
                 << right << setw(12) << "Best/All"
                 << setw(14) << "Maximum"
                 << setw(27) << "Best window [ns]" << '\n';
            for (const EfficiencyRow& row : efficiencyRows) {
                ostringstream counts;
                counts << row.result.bestCount << '/'
                       << row.result.denominator;
                ostringstream maximum;
                maximum << fixed << setprecision(3)
                        << 100.0 * row.result.bestEfficiency << '%';
                ostringstream window;
                window << fixed << setprecision(6) << '['
                       << row.result.bestWindowStart << ", "
                       << row.result.bestWindowEnd << ']';
                cout << "  " << left << setw(20) << row.source
                     << right << setw(12) << counts.str()
                     << setw(14) << maximum.str()
                     << setw(27) << window.str() << '\n';
            }
            for (const auto& [dutID, denominator] :
                 dutEfficiencyDenominators) {
                const size_t matchedHits =
                    dutMatchedHitCounts[dutID];
                const auto residuals =
                    dutEfficiencyT0Residuals.find(dutID);
                const size_t validTimes =
                    residuals ==
                            dutEfficiencyT0Residuals.end()
                        ? 0
                        : residuals->second.residual.size();
                const auto efficiency =
                    dutEfficiencyResults.find(dutID);
                const size_t bestCount =
                    efficiency == dutEfficiencyResults.end()
                        ? 0
                        : efficiency->second.bestCount;
                const size_t spatialMisses =
                    denominator > matchedHits
                        ? denominator - matchedHits
                        : 0;
                const size_t invalidTimes =
                    matchedHits > validTimes
                        ? matchedHits - validTimes
                        : 0;
                const size_t outsideWindow =
                    validTimes > bestCount
                        ? validTimes - bestCount
                        : 0;
                cout << "  DUT " << dutID
                     << " DUTEfficiency selection: "
                     << denominator
                     << " active-area track+T0; losses = "
                     << spatialMisses << " no spatial match + "
                     << invalidTimes << " invalid time + "
                     << outsideWindow
                     << " outside best window\n";
            }
        }

        cout << "[TimeResolution] oscilloscope T0 matched tracks="
             << tracksWithOscilloscopeT0
             << ", measured track-T0 sigma=" << trackT0Fit.sigma
             << " ns, T0-subtracted track resolution="
             << trackT0Resolution << " ns\n";
    }

    if (!dutTiming.samplesByDetector.empty()) {
        TDirectory* dutDirectory = output->mkdir("DUTTimeResolution");
        for (const auto& [dutID, samples] : dutTiming.samplesByDetector) {
            vector<double> dutTimes;
            vector<double> referenceTimes;
            vector<double> residuals;
            const array<string, 9> factorNames = {
                "max_amplitude", "mean_amplitude", "cluster_charge",
                "cluster_size", "cluster_centroid", "cluster_local_x",
                "cluster_local_y", "predicted_x", "predicted_y"};
            const array<string, 9> factorTitles = {
                "Maximum fitted channel amplitude [ADC]",
                "Mean fitted channel amplitude [ADC]",
                "Cluster charge", "Cluster size", "Cluster centroid",
                "Cluster local X", "Cluster local Y",
                "Track predicted X", "Track predicted Y"};
            array<vector<double>, 9> factors;
            dutTimes.reserve(samples.size());
            referenceTimes.reserve(samples.size());
            residuals.reserve(samples.size());
            for (const DUTTimingSample& sample : samples) {
                dutTimes.push_back(sample.dutTime);
                referenceTimes.push_back(sample.trackTime);
                residuals.push_back(sample.residual);
                factors[0].push_back(sample.amplitude);
                factors[1].push_back(sample.meanAmplitude);
                factors[2].push_back(sample.clusterCharge);
                factors[3].push_back(sample.clusterSize);
                factors[4].push_back(sample.clusterCentroid);
                factors[5].push_back(sample.clusterLocalX);
                factors[6].push_back(sample.clusterLocalY);
                factors[7].push_back(sample.predictedX);
                factors[8].push_back(sample.predictedY);
            }
            TDirectory* detectorDirectory =
                dutDirectory->mkdir(("DUT_" + to_string(dutID)).c_str());
            WriteAndFit(
                dutTimes, *detectorDirectory, "dut_time",
                "DUT matched-cluster mean fitted time;t_{DUT} [ns];Events",
                histogramBins);
            WriteAndFit(
                referenceTimes, *detectorDirectory, "track_time",
                "Weighted tracker reference time;t_{track} [ns];Events",
                histogramBins);
            const FitResult measured = WriteAndFit(
                residuals, *detectorDirectory, "dut_minus_track_time",
                "DUT time residual;t_{DUT}-t_{track} [ns];Events",
                histogramBins);
            FitResult dutT0Fit;
            double dutT0Resolution =
                numeric_limits<double>::quiet_NaN();
            const auto dutT0Samples = dutT0Residuals.find(dutID);
            if (dutT0Samples != dutT0Residuals.end()) {
                dutT0Fit = WriteAndFit(
                    dutT0Samples->second.residual,
                    *detectorDirectory,
                    "dut_minus_oscilloscope_t0",
                    "DUT time relative to oscilloscope T0;"
                    "t_{DUT}-T0_{scope} [ns];Events",
                    histogramBins);
                const double referenceVariance =
                    dutT0Samples->second.MeanReferenceVariance();
                const double dutVariance =
                    dutT0Fit.sigma * dutT0Fit.sigma -
                    referenceVariance;
                if (isfinite(dutT0Fit.sigma) &&
                    isfinite(referenceVariance) &&
                    dutVariance >= 0.0)
                    dutT0Resolution = sqrt(dutVariance);
                fillSummary(
                    "oscilloscope_t0_dut", "XYMean",
                    dutID, -1, dutT0Fit, dutT0Resolution);
            }

            TDirectory* correlationDirectory =
                detectorDirectory->mkdir("Correlations");
            TTree correlationSummary(
                "DUTCorrelationSummary",
                "DUT-track residual correlations with cluster factors");
            string factorName;
            Long64_t correlationEntries =
                static_cast<Long64_t>(samples.size());
            Double_t pearson =
                numeric_limits<double>::quiet_NaN();
            correlationSummary.Branch("factor", &factorName);
            correlationSummary.Branch("entries", &correlationEntries);
            correlationSummary.Branch("pearson", &pearson);
            for (size_t factor = 0; factor < factors.size(); ++factor) {
                pearson = WriteDUTCorrelation(
                    factors[factor], residuals,
                    *correlationDirectory, factorNames[factor],
                    factorTitles[factor], min(histogramBins, 80));
                factorName = factorNames[factor];
                correlationSummary.Fill();
            }
            correlationDirectory->cd();
            correlationSummary.Write();

            TTree dutSummary(
                "DUTResolution",
                "DUT resolution after subtracting tracker reference in quadrature");
            Int_t valid = 0;
            Int_t summaryDUTID = dutID;
            Long64_t matchedEntries = measured.entries;
            Double_t residualMeanNs = measured.mean;
            Double_t measuredSigmaNs = measured.sigma;
            Double_t measuredSigmaErrorNs = measured.sigmaError;
            Double_t referenceResolutionNs = trackResolutionNs;
            Double_t referenceResolutionErrorNs = trackResolutionErrorNs;
            Double_t dutVarianceNs2 = numeric_limits<double>::quiet_NaN();
            Double_t dutResolutionNs = numeric_limits<double>::quiet_NaN();
            Double_t dutResolutionErrorNs = numeric_limits<double>::quiet_NaN();
            if (resolutionValid && isfinite(measuredSigmaNs)) {
                dutVarianceNs2 =
                    measuredSigmaNs * measuredSigmaNs -
                    referenceResolutionNs * referenceResolutionNs;
                valid = dutVarianceNs2 >= 0.0;
                if (valid) {
                    dutResolutionNs = sqrt(dutVarianceNs2);
                    if (dutResolutionNs > 0.0 &&
                        isfinite(measuredSigmaErrorNs) &&
                        isfinite(referenceResolutionErrorNs)) {
                        dutResolutionErrorNs =
                            sqrt(pow(measuredSigmaNs * measuredSigmaErrorNs, 2) +
                                 pow(referenceResolutionNs *
                                         referenceResolutionErrorNs,
                                     2)) /
                            dutResolutionNs;
                    }
                }
            }
            dutSummary.Branch("valid", &valid);
            dutSummary.Branch("dutID", &summaryDUTID);
            dutSummary.Branch("entries", &matchedEntries);
            dutSummary.Branch("residualMeanNs", &residualMeanNs);
            dutSummary.Branch("measuredSigmaNs", &measuredSigmaNs);
            dutSummary.Branch("measuredSigmaErrorNs", &measuredSigmaErrorNs);
            dutSummary.Branch("trackResolutionNs", &referenceResolutionNs);
            dutSummary.Branch("trackResolutionErrorNs",
                              &referenceResolutionErrorNs);
            dutSummary.Branch("dutVarianceNs2", &dutVarianceNs2);
            dutSummary.Branch("dutResolutionNs", &dutResolutionNs);
            dutSummary.Branch("dutResolutionErrorNs", &dutResolutionErrorNs);
            detectorDirectory->cd();
            dutSummary.Fill();
            dutSummary.Write();

            ostringstream dutSummaryLine;
            dutSummaryLine << fixed << setprecision(1)
                           << "DUT" << dutID << ' ';
            if (valid) {
                dutSummaryLine << dutResolutionNs << " ns · measured "
                               << measuredSigmaNs << " ns · "
                               << Terminal::Count(samples.size())
                               << " matches";
            } else {
                dutSummaryLine << "resolution invalid · "
                               << Terminal::Count(samples.size())
                               << " matches";
            }
            Terminal::Detail(dutSummaryLine.str());
            if (dutT0Samples != dutT0Residuals.end())
                cout << "[TimeResolution] DUT " << dutID
                     << " vs oscilloscope T0: matched="
                     << dutT0Samples->second.residual.size()
                     << ", measured sigma=" << dutT0Fit.sigma
                     << " ns, T0-subtracted resolution="
                     << dutT0Resolution << " ns\n";
        }
    }

    output->cd();
    summary.Write();
    if (Terminal::Verbose()) {
        Terminal::Detail(Terminal::Muted(
            Terminal::Count(trackTimes.size()) + " tracks · output " +
            outputPath));
    }
    return !trackTimes.empty();
}

void WriteRawTimeHistogram(const vector<double>& values, TDirectory& directory,
                           const string& name, const string& title, int bins) {
    vector<double> finiteValues;
    copy_if(values.begin(), values.end(), back_inserter(finiteValues),
            [](double value) { return isfinite(value); });
    if (finiteValues.empty()) return;
    const auto [robustLow, robustHigh] = RobustRange(finiteValues);
    double low = robustLow;
    double high = robustHigh;
    const double padding = 0.05 * max(1.0e-6, high - low);
    low -= padding;
    high += padding;
    TDirectory::TContext context(&directory);
    TH1D histogram(name.c_str(), title.c_str(), bins, low, high);
    for (double value : finiteValues) histogram.Fill(value);
    histogram.Write();
}

TimingEfficiencyResult WriteCompactTimingEfficiency(
    const vector<double>& residuals, size_t denominator,
    double windowWidthNs, double stepNs, int histogramBins,
    TDirectory& directory) {
    TimingEfficiencyResult result;
    vector<double> sorted;
    copy_if(residuals.begin(), residuals.end(), back_inserter(sorted),
            [](double value) { return isfinite(value); });
    sort(sorted.begin(), sorted.end());
    result.denominator = denominator;
    result.validTimes = sorted.size();
    if (sorted.empty() || denominator == 0 || windowWidthNs <= 0.0 ||
        stepNs <= 0.0)
        return result;

    const auto [low, high] = RobustRange(sorted);
    const double scanLow = low - windowWidthNs;
    const double scanHigh = high;
    const int scanBins = max(
        1, min(20000, static_cast<int>(
                          ceil((scanHigh - scanLow) / stepNs))));
    const double effectiveStep = (scanHigh - scanLow) / scanBins;
    vector<double> efficiencies(static_cast<size_t>(scanBins), 0.0);
    for (int bin = 0; bin < scanBins; ++bin) {
        const double start = scanLow + (bin + 0.5) * effectiveStep;
        const auto first = lower_bound(sorted.begin(), sorted.end(), start);
        const auto last =
            upper_bound(sorted.begin(), sorted.end(), start + windowWidthNs);
        const size_t count = static_cast<size_t>(last - first);
        efficiencies[static_cast<size_t>(bin)] =
            static_cast<double>(count) / denominator;
        if (count > result.bestCount) {
            result.bestCount = count;
            result.bestWindowStart = start;
        }
    }
    result.bestWindowEnd = result.bestWindowStart + windowWidthNs;
    result.bestWindowCenter =
        result.bestWindowStart + 0.5 * windowWidthNs;
    result.bestEfficiency =
        static_cast<double>(result.bestCount) / denominator;
    result.valid = true;

    WriteAndFit(
        sorted, directory, "hTimeResidual",
        "DUT cluster time relative to external T0;"
        "t_{DUT}-T0 [ns];Entries",
        histogramBins);
    TDirectory::TContext context(&directory);
    TH1D efficiency(
        "hWindowEfficiency25ns",
        "DUT 25 ns timing-window efficiency;"
        "Window start [ns];Efficiency",
        scanBins, scanLow, scanHigh);
    efficiency.SetMinimum(0.0);
    efficiency.SetMaximum(1.05);
    for (int bin = 1; bin <= scanBins; ++bin)
        efficiency.SetBinContent(
            bin, efficiencies[static_cast<size_t>(bin - 1)]);
    efficiency.Write();
    return result;
}

bool WriteCompactTimingOutput(
    const string& outputPath,
    const map<TrackKey, map<int, DetectorTimes>>& trackTimes,
    const map<TrackKey, int>& eventIDs,
    const OscilloscopeT0Result& oscilloscopeT0,
    const array<int, 3>& trackerIDs,
    const DUTTimingResult& dutTiming,
    const TrackTimeWeights& trackWeights,
    double timingEfficiencyWindowNs,
    double timingEfficiencyStepNs,
    int histogramBins,
    chrono::steady_clock::time_point analysisStarted) {
    unique_ptr<TFile> output(TFile::Open(outputPath.c_str(), "RECREATE"));
    if (!output || output->IsZombie()) {
        cerr << "[TimeResolution] cannot create " << outputPath << '\n';
        return false;
    }

    Int_t treeEventID = -1;
    ULong64_t treeRawEventID = 0;
    Int_t treeTrackIndex = -1;
    Int_t treeTracker1ID = trackerIDs[0];
    Int_t treeTracker2ID = trackerIDs[1];
    Int_t treeTracker3ID = trackerIDs[2];
    Double_t treeTracker1Time = numeric_limits<double>::quiet_NaN();
    Double_t treeTracker2Time = numeric_limits<double>::quiet_NaN();
    Double_t treeTracker3Time = numeric_limits<double>::quiet_NaN();
    Double_t treeTrackTime = numeric_limits<double>::quiet_NaN();
    Int_t treeDUTID = -1;
    Double_t treeDUTTime = numeric_limits<double>::quiet_NaN();
    Double_t treeDUTMinusTrackTime = numeric_limits<double>::quiet_NaN();
    Bool_t treeHasExternalT0 = false;
    Double_t treeExternalT0 = numeric_limits<double>::quiet_NaN();
    Double_t treeTrackMinusExternalT0 =
        numeric_limits<double>::quiet_NaN();
    Double_t treeDUTMinusExternalT0 =
        numeric_limits<double>::quiet_NaN();
    TTree timingTree("TimingTree", "Per-track detector timing values");
    timingTree.Branch("eventID", &treeEventID);
    timingTree.Branch("rawEventID", &treeRawEventID);
    timingTree.Branch("trackIndex", &treeTrackIndex);
    timingTree.Branch("tracker1ID", &treeTracker1ID);
    timingTree.Branch("tracker2ID", &treeTracker2ID);
    timingTree.Branch("tracker3ID", &treeTracker3ID);
    timingTree.Branch("tracker1Time", &treeTracker1Time);
    timingTree.Branch("tracker2Time", &treeTracker2Time);
    timingTree.Branch("tracker3Time", &treeTracker3Time);
    timingTree.Branch("trackTime", &treeTrackTime);
    timingTree.Branch("dutID", &treeDUTID);
    timingTree.Branch("dutTime", &treeDUTTime);
    timingTree.Branch("dutMinusTrackTime", &treeDUTMinusTrackTime);
    timingTree.Branch("hasExternalT0", &treeHasExternalT0);
    timingTree.Branch("externalT0", &treeExternalT0);
    timingTree.Branch(
        "trackMinusExternalT0", &treeTrackMinusExternalT0);
    timingTree.Branch(
        "dutMinusExternalT0", &treeDUTMinusExternalT0);

    map<TrackKey, vector<const DUTTimingSample*>> dutSamplesByTrack;
    for (const auto& [dutID, samples] : dutTiming.samplesByDetector) {
        (void)dutID;
        for (const DUTTimingSample& sample : samples)
            dutSamplesByTrack[
                {sample.rawEventID, sample.trackIndex}].push_back(&sample);
    }

    map<int, vector<double>> rawDetectorTimes;
    map<int, vector<double>> externalResiduals;
    array<vector<double>, 3> trackerPairResiduals;
    vector<double> trackTimeSamples;
    vector<double> trackExternalResiduals;
    const array<pair<size_t, size_t>, 3> pairIndices = {
        pair<size_t, size_t>{0, 1}, {0, 2}, {1, 2}};

    for (const auto& [key, detectors] : trackTimes) {
        const auto external = oscilloscopeT0.references.find(key.first);
        for (const auto& [detectorID, times] : detectors) {
            const double time = times.value[kTrackerTimingEstimator];
            if (!isfinite(time)) continue;
            rawDetectorTimes[detectorID].push_back(time);
            if (external != oscilloscopeT0.references.end())
                externalResiduals[detectorID].push_back(
                    time - external->second.time);
        }

        array<double, 3> selectedTimes{};
        bool complete = true;
        for (size_t index = 0; index < trackerIDs.size(); ++index) {
            const auto detector = detectors.find(trackerIDs[index]);
            if (detector == detectors.end() ||
                !isfinite(detector->second.value[kTrackerTimingEstimator])) {
                complete = false;
                break;
            }
            selectedTimes[index] =
                detector->second.value[kTrackerTimingEstimator];
        }
        if (!complete) continue;
        for (size_t pairIndex = 0; pairIndex < pairIndices.size();
             ++pairIndex) {
            const auto [first, second] = pairIndices[pairIndex];
            trackerPairResiduals[pairIndex].push_back(
                selectedTimes[first] - selectedTimes[second]);
        }
        const double trackTime =
            inner_product(selectedTimes.begin(), selectedTimes.end(),
                          trackWeights.value.begin(), 0.0);
        trackTimeSamples.push_back(trackTime);
        if (external != oscilloscopeT0.references.end())
            trackExternalResiduals.push_back(
                trackTime - external->second.time);

        treeRawEventID = key.first;
        treeTrackIndex = key.second;
        const auto event = eventIDs.find(key);
        treeEventID = event == eventIDs.end() ? -1 : event->second;
        treeTracker1Time = selectedTimes[0];
        treeTracker2Time = selectedTimes[1];
        treeTracker3Time = selectedTimes[2];
        treeTrackTime = trackTime;
        treeHasExternalT0 =
            external != oscilloscopeT0.references.end();
        treeExternalT0 =
            treeHasExternalT0
                ? external->second.time
                : numeric_limits<double>::quiet_NaN();
        treeTrackMinusExternalT0 =
            treeHasExternalT0
                ? treeTrackTime - treeExternalT0
                : numeric_limits<double>::quiet_NaN();
        const auto dutSamples = dutSamplesByTrack.find(key);
        if (dutSamples == dutSamplesByTrack.end() ||
            dutSamples->second.empty()) {
            treeDUTID = -1;
            treeDUTTime = numeric_limits<double>::quiet_NaN();
            treeDUTMinusTrackTime =
                numeric_limits<double>::quiet_NaN();
            treeDUTMinusExternalT0 =
                numeric_limits<double>::quiet_NaN();
            timingTree.Fill();
        } else {
            for (const DUTTimingSample* sample : dutSamples->second) {
                treeDUTID = sample->detectorID;
                treeDUTTime = sample->dutTime;
                treeDUTMinusTrackTime = sample->residual;
                treeDUTMinusExternalT0 =
                    treeHasExternalT0 && isfinite(treeDUTTime)
                        ? treeDUTTime - treeExternalT0
                        : numeric_limits<double>::quiet_NaN();
                timingTree.Fill();
            }
        }
    }

    map<int, vector<double>> dutTimes;
    map<int, vector<double>> dutMinusTrack;
    map<int, vector<double>> dutExternalResiduals;
    map<int, vector<double>> dutEfficiencyExternalResiduals;
    for (const auto& [dutID, samples] : dutTiming.samplesByDetector) {
        for (const DUTTimingSample& sample : samples) {
            if (isfinite(sample.dutTime))
                dutTimes[dutID].push_back(sample.dutTime);
            if (isfinite(sample.residual))
                dutMinusTrack[dutID].push_back(sample.residual);
            const auto external =
                oscilloscopeT0.references.find(sample.rawEventID);
            if (external == oscilloscopeT0.references.end() ||
                !isfinite(sample.dutTime))
                continue;
            const double residual =
                sample.dutTime - external->second.time;
            dutExternalResiduals[dutID].push_back(residual);
            if (sample.insideActiveArea)
                dutEfficiencyExternalResiduals[dutID].push_back(residual);
        }
    }

    TDirectory* rawDirectory = output->mkdir("RawTimeDistributions");
    for (const auto& [detectorID, samples] : rawDetectorTimes) {
        const auto detector =
            DetectorFactory::GetInstance().GetDetector(detectorID);
        const string label =
            detector ? detector->GetName()
                     : "Detector" + to_string(detectorID);
        WriteRawTimeHistogram(
            samples, *rawDirectory,
            SafeName("hClusterTime_" + label),
            label + " cluster time;t [ns];Entries", histogramBins);
    }
    for (const auto& [dutID, samples] : dutTimes) {
        const string label = "DUT" + to_string(dutID);
        WriteRawTimeHistogram(
            samples, *rawDirectory,
            "hClusterTime_" + label,
            label + " cluster time;t [ns];Entries", histogramBins);
    }
    WriteRawTimeHistogram(
        trackTimeSamples, *rawDirectory, "hTrackTime",
        "Weighted tracker time;t_{track} [ns];Entries", histogramBins);

    TDirectory* internalDirectory =
        output->mkdir("WithoutExternalT0");
    TDirectory* pairDirectory =
        internalDirectory->mkdir("TrackerPairDifferences");
    array<FitResult, 3> pairFits;
    for (size_t pairIndex = 0; pairIndex < pairIndices.size(); ++pairIndex) {
        const auto [first, second] = pairIndices[pairIndex];
        const string firstName =
            DetectorFactory::GetInstance().GetDetector(
                trackerIDs[first])->GetName();
        const string secondName =
            DetectorFactory::GetInstance().GetDetector(
                trackerIDs[second])->GetName();
        pairFits[pairIndex] = WriteAndFit(
            trackerPairResiduals[pairIndex], *pairDirectory,
            SafeName("hDeltaT_" + firstName + "_" + secondName),
            firstName + " - " + secondName +
                ";#Deltat [ns];Entries",
            histogramBins);
    }
    TDirectory* internalDUTDirectory =
        internalDirectory->mkdir("DUT");
    map<int, double> internalDUTResolutions;
    map<int, FitResult> internalDUTFits;
    for (const auto& [dutID, samples] : dutMinusTrack) {
        TDirectory* detectorDirectory =
            internalDUTDirectory->mkdir(
                ("DUT" + to_string(dutID)).c_str());
        const FitResult fit = WriteAndFit(
            samples, *detectorDirectory, "hDUTMinusTrackTime",
            "DUT cluster time relative to track time;"
            "t_{DUT}-t_{track} [ns];Entries",
            histogramBins);
        internalDUTFits[dutID] = fit;
        const double variance =
            fit.sigma * fit.sigma -
            trackWeights.resolution * trackWeights.resolution;
        const double resolution =
            isfinite(variance) && variance >= 0.0
                ? sqrt(variance)
                : numeric_limits<double>::quiet_NaN();
        internalDUTResolutions[dutID] = resolution;
    }

    double externalTrackResolution =
        numeric_limits<double>::quiet_NaN();
    FitResult externalTrackFit;
    map<int, double> externalDUTResolutions;
    map<int, FitResult> externalDetectorFits;
    map<int, FitResult> externalDUTFits;
    map<int, TimingEfficiencyResult> efficiencyResults;
    map<int, size_t> efficiencyDenominators;
    map<int, size_t> spatiallyMatchedTracks;
    if (!oscilloscopeT0.references.empty()) {
        TDirectory* externalDirectory =
            output->mkdir("WithExternalT0");
        for (const auto& [detectorID, samples] : externalResiduals) {
            const auto detector =
                DetectorFactory::GetInstance().GetDetector(detectorID);
            const string label =
                detector ? detector->GetName()
                         : "Detector" + to_string(detectorID);
            TDirectory* detectorDirectory =
                externalDirectory->mkdir(SafeName(label).c_str());
            const FitResult fit = WriteAndFit(
                samples, *detectorDirectory,
                "hClusterTimeMinusT0",
                label + " cluster time relative to external T0;"
                        "t-T0 [ns];Entries",
                histogramBins);
            externalDetectorFits[detectorID] = fit;
        }
        TDirectory* trackDirectory =
            externalDirectory->mkdir("TrackTime");
        externalTrackFit = WriteAndFit(
            trackExternalResiduals, *trackDirectory,
            "hTrackTimeMinusT0",
            "Track time relative to external T0;"
            "t_{track}-T0 [ns];Entries",
            histogramBins);
        externalTrackResolution = externalTrackFit.sigma;

        for (const auto& [dutID, samples] : dutExternalResiduals) {
            TDirectory* detectorDirectory =
                externalDirectory->mkdir(
                    ("DUT" + to_string(dutID)).c_str());
            const FitResult fit = WriteAndFit(
                samples, *detectorDirectory,
                "hClusterTimeMinusT0",
                "DUT cluster time relative to external T0;"
                "t_{DUT}-T0 [ns];Entries",
                histogramBins);
            externalDUTFits[dutID] = fit;
            externalDUTResolutions[dutID] = fit.sigma;
        }

        TDirectory* efficiencyRoot =
            externalDirectory->mkdir("DUTTimingEfficiency");
        for (const auto& [dutID, residuals] :
             dutEfficiencyExternalResiduals) {
            size_t denominator = 0;
            const auto cases = dutTiming.activeAreaTrackCases.find(dutID);
            if (cases != dutTiming.activeAreaTrackCases.end()) {
                denominator = static_cast<size_t>(count_if(
                    cases->second.begin(), cases->second.end(),
                    [&](const TrackKey& key) {
                        return oscilloscopeT0.references.count(
                                   key.first) != 0;
                    }));
            }
            efficiencyDenominators[dutID] = denominator;
            const auto matchedCases =
                dutTiming.activeAreaMatchedHitCases.find(dutID);
            if (matchedCases != dutTiming.activeAreaMatchedHitCases.end()) {
                spatiallyMatchedTracks[dutID] =
                    static_cast<size_t>(count_if(
                        matchedCases->second.begin(),
                        matchedCases->second.end(),
                        [&](const TrackKey& key) {
                            return oscilloscopeT0.references.count(
                                       key.first) != 0;
                        }));
            }
            TDirectory* detectorDirectory =
                efficiencyRoot->mkdir(
                    ("DUT" + to_string(dutID)).c_str());
            const TimingEfficiencyResult result =
                WriteCompactTimingEfficiency(
                    residuals, denominator,
                    timingEfficiencyWindowNs, timingEfficiencyStepNs,
                    histogramBins, *detectorDirectory);
            if (result.valid) efficiencyResults[dutID] = result;
        }
    }

    string resultAnalysis;
    string resultObject;
    Int_t resultDetectorA = -1;
    Int_t resultDetectorB = -1;
    Long64_t resultEntries = 0;
    Double_t resultMean = numeric_limits<double>::quiet_NaN();
    Double_t resultSigma = numeric_limits<double>::quiet_NaN();
    Double_t resultSigmaError = numeric_limits<double>::quiet_NaN();
    Double_t resultResolution = numeric_limits<double>::quiet_NaN();
    Long64_t resultNumerator = 0;
    Long64_t resultDenominator = 0;
    Double_t resultEfficiency = numeric_limits<double>::quiet_NaN();
    Double_t resultWindowStart = numeric_limits<double>::quiet_NaN();
    Double_t resultWindowEnd = numeric_limits<double>::quiet_NaN();
    TTree resultsTree("TimingResults", "Timing fit and efficiency results");
    resultsTree.Branch("analysis", &resultAnalysis);
    resultsTree.Branch("object", &resultObject);
    resultsTree.Branch("detectorA", &resultDetectorA);
    resultsTree.Branch("detectorB", &resultDetectorB);
    resultsTree.Branch("entries", &resultEntries);
    resultsTree.Branch("meanNs", &resultMean);
    resultsTree.Branch("sigmaNs", &resultSigma);
    resultsTree.Branch("sigmaErrorNs", &resultSigmaError);
    resultsTree.Branch("resolutionNs", &resultResolution);
    resultsTree.Branch("numerator", &resultNumerator);
    resultsTree.Branch("denominator", &resultDenominator);
    resultsTree.Branch("efficiency", &resultEfficiency);
    resultsTree.Branch("windowStartNs", &resultWindowStart);
    resultsTree.Branch("windowEndNs", &resultWindowEnd);
    const auto fillResult = [&](const string& analysis,
                                const string& object, int detectorA,
                                int detectorB, const FitResult& fit,
                                double resolution) {
        resultAnalysis = analysis;
        resultObject = object;
        resultDetectorA = detectorA;
        resultDetectorB = detectorB;
        resultEntries = fit.entries;
        resultMean = fit.mean;
        resultSigma = fit.sigma;
        resultSigmaError = fit.sigmaError;
        resultResolution = resolution;
        resultNumerator = 0;
        resultDenominator = 0;
        resultEfficiency = numeric_limits<double>::quiet_NaN();
        resultWindowStart = numeric_limits<double>::quiet_NaN();
        resultWindowEnd = numeric_limits<double>::quiet_NaN();
        resultsTree.Fill();
    };
    for (size_t index = 0; index < pairFits.size(); ++index) {
        const auto [first, second] = pairIndices[index];
        fillResult("without_external_t0", "tracker_pair",
                   trackerIDs[first], trackerIDs[second],
                   pairFits[index], pairFits[index].sigma);
    }
    for (size_t index = 0; index < trackerIDs.size(); ++index) {
        FitResult derived;
        derived.entries = static_cast<long long>(
            min({trackerPairResiduals[0].size(),
                 trackerPairResiduals[1].size(),
                 trackerPairResiduals[2].size()}));
        derived.sigma = sqrt(trackWeights.detectorVariance[index]);
        fillResult("without_external_t0", "tracker_resolution",
                   trackerIDs[index], -1, derived, derived.sigma);
    }
    {
        FitResult derived;
        derived.entries = static_cast<long long>(trackTimeSamples.size());
        derived.sigma = trackWeights.resolution;
        fillResult("without_external_t0", "track_time", -1, -1,
                   derived, trackWeights.resolution);
    }
    for (const auto& [dutID, fit] : internalDUTFits)
        fillResult("without_external_t0", "dut", dutID, -1, fit,
                   internalDUTResolutions[dutID]);
    for (const auto& [detectorID, fit] : externalDetectorFits)
        fillResult("with_external_t0", "detector", detectorID, -1,
                   fit, fit.sigma);
    if (isfinite(externalTrackResolution)) {
        fillResult("with_external_t0", "track_time", -1, -1,
                   externalTrackFit,
                   externalTrackResolution);
    }
    for (const auto& [dutID, fit] : externalDUTFits)
        fillResult("with_external_t0", "dut", dutID, -1, fit,
                   fit.sigma);
    for (const auto& [dutID, efficiency] : efficiencyResults) {
        resultAnalysis = "with_external_t0";
        resultObject = "dut_timing_window";
        resultDetectorA = dutID;
        resultDetectorB = -1;
        resultEntries = static_cast<Long64_t>(efficiency.validTimes);
        resultMean = numeric_limits<double>::quiet_NaN();
        resultSigma = numeric_limits<double>::quiet_NaN();
        resultSigmaError = numeric_limits<double>::quiet_NaN();
        resultResolution = numeric_limits<double>::quiet_NaN();
        resultNumerator = static_cast<Long64_t>(efficiency.bestCount);
        resultDenominator =
            static_cast<Long64_t>(efficiency.denominator);
        resultEfficiency = efficiency.bestEfficiency;
        resultWindowStart = efficiency.bestWindowStart;
        resultWindowEnd = efficiency.bestWindowEnd;
        resultsTree.Fill();
    }
    output->cd();
    timingTree.Write();
    resultsTree.Write();
    output->Close();

    const auto row = [](const string& label, const string& value) {
        cout << "  " << left << setw(38) << label << right << setw(39)
             << value << '\n';
    };
    const auto separator = [] { cout << string(100, '-') << '\n'; };
    const auto count = [](size_t value) { return Terminal::Count(value); };
    const auto fixedValue = [](double value, int precision,
                               const string& suffix = "") {
        if (!isfinite(value)) return string("—");
        ostringstream text;
        text << fixed << setprecision(precision) << value << suffix;
        return text.str();
    };
    cout << string(100, '=') << '\n'
         << Terminal::Accent("Configuration") << '\n';
    row("External reference",
        oscilloscopeT0.references.empty() ? "disabled"
                                          : "oscilloscope T0");
    row("Timing-window width",
        fixedValue(timingEfficiencyWindowNs, 2, " ns"));
    if (!oscilloscopeT0.references.empty()) {
        separator();
        cout << Terminal::Accent("Oscilloscope T0 Matching") << '\n';
        {
            const double fraction =
                oscilloscopeT0.inputEntries > 0
                    ? 100.0 * oscilloscopeT0.references.size() /
                          oscilloscopeT0.inputEntries
                    : 0.0;
            ostringstream value;
            value << count(oscilloscopeT0.inputEntries) << " / "
                  << count(oscilloscopeT0.references.size()) << "  ("
                  << fixed << setprecision(2) << fraction << "%)";
            row("Input / matched raw event IDs", value.str());
        }
        row("Ambiguous event IDs",
            count(oscilloscopeT0.duplicateEventIDs));
        row("Associated tracks", count(trackExternalResiduals.size()));
        {
            ostringstream weights;
            weights << fixed << setprecision(2)
                    << oscilloscopeT0.channelWeight[0] << " / "
                    << oscilloscopeT0.channelWeight[1] << " / "
                    << oscilloscopeT0.channelWeight[2];
            row("T0 channel weights", weights.str());
        }
        row("Combined T0 resolution",
            fixedValue(oscilloscopeT0.resolution, 2, " ns"));

        for (const auto& [dutID, result] : efficiencyResults) {
            separator();
            cout << Terminal::Accent(
                        "DUT" + to_string(dutID) +
                        " Timing-Window Efficiency")
                 << '\n';
            const size_t denominator = efficiencyDenominators[dutID];
            const size_t matched = spatiallyMatchedTracks[dutID];
            row("Active-area tracks with valid T0", count(denominator));
            row("Spatially matched tracks", count(matched));
            row("Tracks inside best " +
                    fixedValue(timingEfficiencyWindowNs, 0) +
                    " ns window",
                count(result.bestCount));
            cout << '\n';
            const auto efficiencyText = [&](size_t passed, size_t total) {
                ostringstream text;
                const double efficiency =
                    total > 0 ? 100.0 * passed / total : 0.0;
                text << count(passed) << " / " << count(total) << " = "
                     << fixed << setprecision(2) << efficiency << '%';
                return text.str();
            };
            row("Spatial-match efficiency",
                efficiencyText(matched, denominator));
            row("Conditional timing efficiency",
                efficiencyText(result.bestCount, matched));
            row("Combined efficiency",
                efficiencyText(result.bestCount, denominator));
            cout << '\n';
            row("Losses",
                count(denominator > matched ? denominator - matched : 0) +
                    " no spatial match");
            row("",
                count(matched > result.bestCount
                          ? matched - result.bestCount
                          : 0) +
                    " outside timing window");
            row("Best " + fixedValue(timingEfficiencyWindowNs, 0) +
                    " ns window",
                "[" + fixedValue(result.bestWindowStart, 2) + ", " +
                    fixedValue(result.bestWindowEnd, 2) + "] ns");
        }
    }

    separator();
    cout << Terminal::Accent("Time Resolution") << '\n'
         << "  " << left << setw(31) << "Detector / Reference" << right
         << setw(23) << "Self-calibration" << setw(20) << "External T0"
         << '\n';
    const auto resolutionRow = [&](const string& label,
                                   const string& self,
                                   const string& external) {
        cout << "  " << left << setw(31) << label << right << setw(23)
             << self << setw(20) << external << '\n';
    };
    for (size_t index = 0; index < trackerIDs.size(); ++index) {
        const auto detector =
            DetectorFactory::GetInstance().GetDetector(trackerIDs[index]);
        resolutionRow(
            detector ? detector->GetName()
                     : "Tracker" + to_string(index + 1),
            fixedValue(sqrt(trackWeights.detectorVariance[index]), 2,
                       " ns"),
            "—");
    }
    resolutionRow(
        "Combined tracker reference",
        fixedValue(trackWeights.resolution, 2, " ns"),
        fixedValue(externalTrackResolution, 2, " ns"));
    for (const auto& [dutID, selfResolution] : internalDUTResolutions) {
        const auto external = externalDUTResolutions.find(dutID);
        resolutionRow(
            "DUT" + to_string(dutID),
            fixedValue(selfResolution, 2, " ns"),
            external == externalDUTResolutions.end()
                ? "—"
                : fixedValue(external->second, 2, " ns"));
    }
    separator();
    cout << Terminal::Accent("Status") << '\n';
    row(Terminal::Success("[PASS] Time Resolution completed"), "");
    const double runtime = chrono::duration<double>(
                               chrono::steady_clock::now() -
                               analysisStarted)
                               .count();
    row("Runtime", fixedValue(runtime, 1, " s"));
    cout << string(100, '=') << '\n';
    return !trackTimes.empty();
}

}  // namespace

void TimeResolutionScript::LoadConfig(const json& config) {
    m_trackFile = config.value("trackFile", m_trackFile);
    m_dutFile = config.value("dutFile", m_dutFile);
    m_oscilloscopeFile =
        config.value("oscilloscopeFile", m_oscilloscopeFile);
    m_outputFile = config.value("outputFile", m_outputFile);
    m_histogramBins = config.value("histogramBins", m_histogramBins);
    m_trackerIDs = config.value("trackerIDs", m_trackerIDs);
    m_analyzeDUTTiming =
        config.value("analyzeDUTTiming", m_analyzeDUTTiming);
    m_timingEfficiencyWindowNs = config.value(
        "timingEfficiencyWindowNs", m_timingEfficiencyWindowNs);
    m_timingEfficiencyStepNs = config.value(
        "timingEfficiencyStepNs", m_timingEfficiencyStepNs);
    if (config.contains("timingWaveform"))
        m_timingWaveformConfig = config["timingWaveform"];
}

void TimeResolutionScript::Print() const {
    if (!Terminal::Verbose()) return;
    cout << "TimeResolutionScript:\n";
    if (!m_trackerIDs.empty()) {
        cout << "  three-tracker IDs=";
        for (size_t i = 0; i < m_trackerIDs.size(); ++i)
            cout << (i == 0 ? "" : ",") << m_trackerIDs[i];
        cout << '\n';
    }
    cout << "  DUT timing=" << (m_analyzeDUTTiming ? "enabled" : "disabled")
         << '\n';
    cout << "  oscilloscope T0 file="
         << (m_oscilloscopeFile.empty() ? "disabled"
                                        : m_oscilloscopeFile)
         << '\n';
    cout << "  DUT timing-efficiency scan=" << m_timingEfficiencyWindowNs
         << " ns window, " << m_timingEfficiencyStepNs << " ns step\n";
    cout << "  timing waveform=" << m_timingWaveformConfig.dump() << '\n';
}

bool TimeResolutionScript::Validate() const {
    return m_histogramBins > 0 &&
           isfinite(m_timingEfficiencyWindowNs) &&
           m_timingEfficiencyWindowNs > 0.0 &&
           isfinite(m_timingEfficiencyStepNs) &&
           m_timingEfficiencyStepNs > 0.0 &&
           (m_trackerIDs.empty() || m_trackerIDs.size() == 3) &&
           set<int>(m_trackerIDs.begin(), m_trackerIDs.end()).size() ==
               m_trackerIDs.size();
}

bool TimeResolutionScript::Execute() {
    const auto analysisStarted = chrono::steady_clock::now();
    const string trackPath = m_trackFile.empty() ? GetOutputDir() + "TrackInfo.root" : m_trackFile;
    const string outputPath = filesystem::path(m_outputFile).is_absolute()
                                  ? m_outputFile
                                  : GetOutputDir() + m_outputFile;
    vector<int> configuredTrackerIDs = m_trackerIDs;
    if (configuredTrackerIDs.empty()) {
        configuredTrackerIDs =
            DetectorFactory::GetInstance().GetDetectorIDsByRole(Detector::Role::Tracker);
    }
    if (configuredTrackerIDs.size() != 3) {
        cerr << "[TimeResolution] exactly three tracker IDs are required; found "
             << configuredTrackerIDs.size()
             << ". Set scripts[].config.trackerIDs explicitly.\n";
        return false;
    }
    for (int detectorID : configuredTrackerIDs) {
        const auto detector = DetectorFactory::GetInstance().GetDetector(detectorID);
        if (!detector || !detector->isTracker()) {
            cerr << "[TimeResolution] detector " << detectorID
                 << " is missing or is not configured with role Tracker\n";
            return false;
        }
    }
    const array<int, 3> trackerIDs = {
        configuredTrackerIDs[0], configuredTrackerIDs[1], configuredTrackerIDs[2]};

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
    size_t waveformFits = 0;
    if (!LoadReconstructionTimes(
            *validation, *parser, m_timingWaveformConfig,
            waveformFits, reconstruction))
        return false;
    if (Terminal::Verbose()) {
        Terminal::Detail(
            "loaded " + Terminal::Count(reconstruction.trackTimes.size()) +
            " tracks · " +
            Terminal::Count(reconstruction.wantedEventIDs.size()) +
            " raw event IDs");
    }
    const TrackTimeWeights trackWeights =
        CalculateTrackTimeWeights(reconstruction.trackTimes, trackerIDs);
    if (!trackWeights.valid) {
        cerr << "[TimeResolution] cannot derive positive tracker variances "
                "for inverse-variance track-time weights\n";
        return false;
    }
    OscilloscopeT0Result oscilloscopeT0;
    if (!m_oscilloscopeFile.empty()) {
        const string oscilloscopePath =
            filesystem::path(m_oscilloscopeFile).is_absolute()
                ? m_oscilloscopeFile
                : GetOutputDir() + m_oscilloscopeFile;
        if (filesystem::exists(oscilloscopePath))
            oscilloscopeT0 = LoadOscilloscopeT0(
                oscilloscopePath, reconstruction.wantedEventIDs);
        else if (Terminal::Verbose())
            Terminal::Detail(
                Terminal::Muted("T0 skipped · no oscilloscope result"));
    }
    DUTTimingResult dutTiming;
    if (m_analyzeDUTTiming) {
        const string dutPath = filesystem::path(m_dutFile).is_absolute()
                                   ? m_dutFile
                                   : GetOutputDir() + m_dutFile;
        unique_ptr<TFile> dutFile(TFile::Open(dutPath.c_str(), "READ"));
        if (!dutFile || dutFile->IsZombie()) {
            cerr << "[TimeResolution] cannot open DUT file " << dutPath << '\n';
            return false;
        }
        auto* dutTree = dynamic_cast<TTree*>(dutFile->Get("PadDUTTree"));
        if (!dutTree) {
            cerr << "[TimeResolution] PadDUTTree is missing in " << dutPath
                 << '\n';
            return false;
        }
        const auto references = BuildTrackTimeReferences(
            reconstruction.trackTimes, reconstruction.eventIDs, trackerIDs,
            trackWeights);
        dutTiming = LoadDUTTiming(
            *dutTree, *parser, references, m_timingWaveformConfig,
            waveformFits);
    }

    const auto& trackTimes = reconstruction.trackTimes;

    const bool wroteOutput =
        WriteCompactTimingOutput(
            outputPath, trackTimes, reconstruction.eventIDs,
            oscilloscopeT0, trackerIDs, dutTiming, trackWeights,
            m_timingEfficiencyWindowNs,
            m_timingEfficiencyStepNs, m_histogramBins, analysisStarted);
    if (Terminal::Verbose()) {
        Terminal::Detail(Terminal::Muted(
            Terminal::Count(waveformFits) + " waveform fits"));
    }
    return wroteOutput;
}

REGISTER_SCRIPT("TimeResolution", TimeResolutionScript);
