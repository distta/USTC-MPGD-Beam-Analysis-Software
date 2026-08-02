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
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLine.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <filesystem>
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
#include <unordered_set>
#include <utility>
#include <vector>

using namespace std;

namespace {

const array<string, 7> kTimeMethods = {
    "XMean", "XMin", "YMean", "YMin", "XYMean", "XYMin",
    "XYWeighted"};

string CanonicalTimeMethod(string method) {
    transform(method.begin(), method.end(), method.begin(),
              [](unsigned char value) {
                  return static_cast<char>(tolower(value));
              });
    for (const string& candidate : kTimeMethods) {
        string lowered = candidate;
        transform(lowered.begin(), lowered.end(), lowered.begin(),
                  [](unsigned char value) {
                      return static_cast<char>(tolower(value));
                  });
        if (method == lowered) return candidate;
    }
    return {};
}

double CalculateConfiguredTime(
    const vector<double>& times,
    const vector<double>& errors,
    const vector<double>& amplitudes,
    const vector<int>& planes,
    const string& method,
    double amplitudeWeightPower) {
    const bool useX = method == "XMean" || method == "XMin" ||
                      method == "XYMean" || method == "XYMin" ||
                      method == "XYWeighted";
    const bool useY = method == "YMean" || method == "YMin" ||
                      method == "XYMean" || method == "XYMin" ||
                      method == "XYWeighted";
    const bool useMinimum = method == "XMin" || method == "YMin" ||
                            method == "XYMin";
    const bool useAmplitudeWeight = method == "XYWeighted";
    const size_t entries = min(times.size(), planes.size());
    double minimumTime = numeric_limits<double>::infinity();
    double timeSum = 0.0;
    size_t validTimes = 0;
    double weightedTime = 0.0;
    double weightSum = 0.0;
    for (size_t index = 0; index < entries; ++index) {
        const bool selected =
            (planes[index] == 0 && useX) ||
            (planes[index] == 1 && useY);
        if (!selected || !isfinite(times[index])) continue;
        minimumTime = min(minimumTime, times[index]);
        timeSum += times[index];
        ++validTimes;
        if (useAmplitudeWeight) continue;
        if (index >= errors.size() || !isfinite(errors[index]) ||
            errors[index] <= 0.0)
            continue;
        const double weight = 1.0 / (errors[index] * errors[index]);
        if (!isfinite(weight)) continue;
        weightedTime += weight * times[index];
        weightSum += weight;
    }
    if (validTimes == 0)
        return numeric_limits<double>::quiet_NaN();
    if (useMinimum) return minimumTime;
    if (useAmplitudeWeight) {
        weightedTime = 0.0;
        weightSum = 0.0;
        for (size_t index = 0; index < entries; ++index) {
            const bool selected =
                (planes[index] == 0 && useX) ||
                (planes[index] == 1 && useY);
            if (!selected || !isfinite(times[index]) ||
                index >= amplitudes.size() ||
                !isfinite(amplitudes[index]) || amplitudes[index] <= 0.0)
                continue;
            const double weight = pow(amplitudes[index],
                                      amplitudeWeightPower);
            if (!isfinite(weight)) continue;
            weightedTime += weight * times[index];
            weightSum += weight;
        }
        if (weightSum > 0.0) return weightedTime / weightSum;
    }
    if (weightSum > 0.0) return weightedTime / weightSum;
    return timeSum / static_cast<double>(validTimes);
}

// Reconstruction event IDs are unique within the analysis input, while raw
// event IDs can wrap or repeat.  Keep the latter only as external-T0 lookup
// metadata so equal raw IDs never merge two reconstructed tracks.
using TrackKey = pair<int, int>;

struct DetectorTimes {
    double clusterTime = numeric_limits<double>::quiet_NaN();
    vector<double> stripTimes;
    vector<double> stripTimeErrors;
    vector<double> stripAmplitudes;
    vector<double> stripRiseTimes;
    vector<double> stripWidths;
    vector<int> stripIDs;
    vector<int> stripPlanes;
};

struct FitResult {
    long long entries = 0;
    double mean = numeric_limits<double>::quiet_NaN();
    double sigma = numeric_limits<double>::quiet_NaN();
    double sigmaError = numeric_limits<double>::quiet_NaN();
};

struct ReconstructionTimingResult {
    map<TrackKey, map<int, DetectorTimes>> trackTimes;
    map<TrackKey, uint64_t> rawEventIDs;
    set<uint64_t> wantedEventIDs;
    set<uint64_t> ambiguousRawEventIDs;
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
    double dutTime = numeric_limits<double>::quiet_NaN();
    double residual = numeric_limits<double>::quiet_NaN();
    bool insideActiveArea = false;
    vector<double> stripTimes;
    vector<double> stripTimeErrors;
    vector<double> stripAmplitudes;
    vector<double> stripRiseTimes;
    vector<double> stripWidths;
    vector<int> stripIDs;
    vector<int> stripPlanes;
};

struct DUTTimingResult {
    map<int, vector<DUTTimingSample>> samplesByDetector;
    map<int, set<TrackKey>> activeAreaTrackCases;
    map<int, set<TrackKey>> activeAreaMatchedHitCases;
};

enum class DUTTreeSchema {
    Planar,
    Pad
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

double SubtractReferenceResolution(
    const FitResult& fit, const T0ResidualSamples& samples) {
    const double referenceVariance = samples.MeanReferenceVariance();
    if (!isfinite(fit.sigma) || !isfinite(referenceVariance))
        return numeric_limits<double>::quiet_NaN();
    const double intrinsicVariance =
        fit.sigma * fit.sigma - referenceVariance;
    return intrinsicVariance >= 0.0
               ? sqrt(intrinsicVariance)
               : numeric_limits<double>::quiet_NaN();
}

struct TimeWalkCorrection {
    bool valid = false;
    size_t entries = 0;
    double amplitudeLow = numeric_limits<double>::quiet_NaN();
    double amplitudeHigh = numeric_limits<double>::quiet_NaN();
    array<double, 4> parameter{};

    double Evaluate(double amplitude) const {
        if (!valid || !isfinite(amplitude))
            return numeric_limits<double>::quiet_NaN();
        const double x = clamp(amplitude, amplitudeLow, amplitudeHigh);
        return parameter[0] +
               x * (parameter[1] +
                    x * (parameter[2] + x * parameter[3]));
    }
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

// ==========================
// Reconstruction-based timing
// ==========================

DetectorTimes ExtractTimes(const vector<int>& clusterIndices,
                           const vector<Cluster>& clusters,
                           const vector<ChannelHit>& channelHits,
                           const vector<RawData>& detectorRawData,
                           HitProcessor& timingFitter,
                           const string& timeMethod,
                           double amplitudeWeightPower,
                           size_t& fitCount) {
    DetectorTimes result;
    map<int, ChannelHit> fittedStripHits;

    for (int clusterIndex : clusterIndices) {
        if (clusterIndex < 0 || clusterIndex >= static_cast<int>(clusters.size())) continue;
        const Cluster& cluster = clusters[clusterIndex];

        for (int stripIndex : cluster.channelHitIndices) {
            if (stripIndex < 0 || stripIndex >= static_cast<int>(channelHits.size())) continue;
            const int rawIndex = channelHits[stripIndex].rawIndices;
            if (rawIndex < 0 || rawIndex >= static_cast<int>(detectorRawData.size())) continue;

            auto cached = fittedStripHits.find(rawIndex);
            if (cached == fittedStripHits.end()) {
                ChannelHit fitted = timingFitter.ProcessHit(
                    detectorRawData[rawIndex]);
                ++fitCount;
                cached = fittedStripHits.emplace(rawIndex, fitted).first;
            }
        }
    }

    result.stripTimes.reserve(fittedStripHits.size());
    result.stripTimeErrors.reserve(fittedStripHits.size());
    result.stripAmplitudes.reserve(fittedStripHits.size());
    result.stripRiseTimes.reserve(fittedStripHits.size());
    result.stripWidths.reserve(fittedStripHits.size());
    result.stripIDs.reserve(fittedStripHits.size());
    result.stripPlanes.reserve(fittedStripHits.size());
    for (const auto& [rawIndex, hit] : fittedStripHits) {
        (void)rawIndex;
        if (!hit.isValid || !isfinite(hit.time)) continue;
        result.stripTimes.push_back(hit.time);
        result.stripTimeErrors.push_back(hit.timeError);
        result.stripAmplitudes.push_back(hit.amp);
        result.stripRiseTimes.push_back(hit.riseTime);
        result.stripWidths.push_back(hit.width);
        result.stripIDs.push_back(hit.id0);
        result.stripPlanes.push_back(hit.type);
    }
    result.clusterTime = CalculateConfiguredTime(
        result.stripTimes, result.stripTimeErrors,
        result.stripAmplitudes, result.stripPlanes, timeMethod,
        amplitudeWeightPower);
    return result;
}

bool LoadReconstructionTimes(TTree& validation, RawDataParser& parser,
                             const json& timingWaveformConfig,
                             const string& timeMethod,
                             double amplitudeWeightPower,
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
    map<uint64_t, int> rawToEventID;
    const Long64_t entries = validation.GetEntries();
    TimingProgress progress("Tracker timing", entries, fitCount);
    for (Long64_t entry = 0; entry < entries; ++entry) {
        progress.Update(entry);
        validation.GetEntry(entry);
        const TrackKey key{eventID, trackIndex};
        const auto [rawEntry, inserted] =
            result.rawEventIDs.emplace(key, rawEventID);
        if (!inserted && rawEntry->second != rawEventID) {
            cerr << "[TimeResolution] inconsistent rawEventID for event "
                 << eventID << ", track " << trackIndex << '\n';
            return false;
        }
        result.wantedEventIDs.insert(rawEventID);
        const auto [eventEntry, newRawID] =
            rawToEventID.emplace(rawEventID, eventID);
        if (!newRawID && eventEntry->second != eventID)
            result.ambiguousRawEventIDs.insert(rawEventID);
        if (!clusterIndices || !channelHits || !clusters) continue;
        if (eventID != loadedEventID) {
            rawHits = parser.LoadEvent(eventID);
            loadedEventID = eventID;
        }
        const auto rawDetector = rawHits.find(detID);
        if (rawDetector == rawHits.end()) continue;
        DetectorTimes times =
            ExtractTimes(*clusterIndices, *clusters, *channelHits,
                         rawDetector->second, timingFitter,
                         timeMethod, amplitudeWeightPower, fitCount);
        result.trackTimes[key][detID] = times;
    }
    for (uint64_t rawEventID : result.ambiguousRawEventIDs)
        result.wantedEventIDs.erase(rawEventID);
    progress.Update(entries);
    return true;
}

map<int, TrackTimeReference> BuildTrackTimeReferences(
    const map<TrackKey, map<int, DetectorTimes>>& trackTimes,
    const map<TrackKey, uint64_t>& rawEventIDs,
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
                !isfinite(detector->second.clusterTime)) {
                valid = false;
                break;
            }
            times[i] = detector->second.clusterTime;
        }
        const auto rawEvent = rawEventIDs.find(key);
        if (!valid || rawEvent == rawEventIDs.end()) continue;
        if (references.count(key.first)) {
            ambiguousEventIDs.insert(key.first);
            continue;
        }
        double trackTime = 0.0;
        for (size_t i = 0; i < trackerIDs.size(); ++i) {
            trackTime += weights.value[i] * times[i];
        }
        references[key.first] = {
            rawEvent->second, key.second, trackTime};
    }
    for (int eventID : ambiguousEventIDs) references.erase(eventID);
    return references;
}

bool IsInsideDUTActiveArea(
    const shared_ptr<Detector>& detector,
    double predictedX, double predictedY) {
    if (!detector || !isfinite(predictedX) || !isfinite(predictedY))
        return false;

    if (const auto* pad = detector->GetPlanarPadConfig()) {
        const double xMinimum = -0.5 * pad->pitchX;
        const double xMaximum =
            (pad->columns - 0.5) * pad->pitchX;
        const double yMinimum = -0.5 * pad->pitchY;
        const double yMaximum =
            (pad->rows - 0.5) * pad->pitchY;
        return predictedX >= xMinimum && predictedX < xMaximum &&
               predictedY >= yMinimum && predictedY < yMaximum;
    }

    const auto* planar = detector->GetPlanarConfig();
    if (!planar || planar->readoutPlaneType.empty()) return false;
    constexpr double degreesToRadians =
        3.14159265358979323846 / 180.0;
    for (int type : planar->readoutPlaneType) {
        const auto angle = planar->readoutPlaneAngle.find(type);
        const auto pitch = planar->readoutPlanePitch.find(type);
        const auto strips = planar->readoutPlaneStripNumber.find(type);
        if (angle == planar->readoutPlaneAngle.end() ||
            pitch == planar->readoutPlanePitch.end() ||
            strips == planar->readoutPlaneStripNumber.end() ||
            !isfinite(pitch->second) || pitch->second <= 0.0 ||
            strips->second <= 0) {
            return false;
        }
        const double radians = angle->second * degreesToRadians;
        const double coordinate =
            predictedX * cos(radians) + predictedY * sin(radians);
        const double minimum = -0.5 * pitch->second;
        const double maximum =
            (strips->second - 0.5) * pitch->second;
        if (coordinate < minimum || coordinate >= maximum) return false;
    }
    return true;
}

DUTTimingResult LoadDUTTiming(
    TTree& dutTree, RawDataParser& parser,
    const map<int, TrackTimeReference>& trackReferences,
    const json& timingWaveformConfig,
    const string& timeMethod,
    double amplitudeWeightPower,
    size_t& fitCount) {
    DUTTimingResult result;
    for (const char* branch :
         {"eventID", "dutID", "hitFlag", "predX", "predY"}) {
        if (!dutTree.GetBranch(branch)) {
            cerr << "[TimeResolution] DUT tree branch " << branch
                 << " is missing\n";
            return result;
        }
    }

    const auto hasBranches = [&](initializer_list<const char*> names) {
        return all_of(
            names.begin(), names.end(),
            [&](const char* name) { return dutTree.GetBranch(name); });
    };
    const bool hasPlanarSchema = hasBranches(
        {"clusterX", "clusterY", "channelHitsX", "channelHitsY"});
    const bool hasPadSchema =
        hasBranches({"selectedCluster", "selectedChannelHits"});
    if (!hasPlanarSchema && !hasPadSchema) {
        cerr << "[TimeResolution] DUT tree has neither the planar nor the "
                "planar_pad timing branches\n";
        return result;
    }
    const DUTTreeSchema schema =
        hasPlanarSchema ? DUTTreeSchema::Planar : DUTTreeSchema::Pad;

    Int_t eventID = 0, dutID = 0, hitFlag = 0;
    Double_t predictedX = numeric_limits<double>::quiet_NaN();
    Double_t predictedY = numeric_limits<double>::quiet_NaN();
    Cluster* selectedCluster = nullptr;
    vector<ChannelHit>* selectedChannelHits = nullptr;
    Cluster* clusterX = nullptr;
    Cluster* clusterY = nullptr;
    vector<ChannelHit>* channelHitsX = nullptr;
    vector<ChannelHit>* channelHitsY = nullptr;
    dutTree.SetBranchAddress("eventID", &eventID);
    dutTree.SetBranchAddress("dutID", &dutID);
    dutTree.SetBranchAddress("hitFlag", &hitFlag);
    dutTree.SetBranchAddress("predX", &predictedX);
    dutTree.SetBranchAddress("predY", &predictedY);
    if (schema == DUTTreeSchema::Planar) {
        dutTree.SetBranchAddress("clusterX", &clusterX);
        dutTree.SetBranchAddress("clusterY", &clusterY);
        dutTree.SetBranchAddress("channelHitsX", &channelHitsX);
        dutTree.SetBranchAddress("channelHitsY", &channelHitsY);
    } else {
        dutTree.SetBranchAddress("selectedCluster", &selectedCluster);
        dutTree.SetBranchAddress(
            "selectedChannelHits", &selectedChannelHits);
    }

    HitProcessor timingFitter;
    timingFitter.LoadConfig(timingWaveformConfig);
    int loadedEventID = numeric_limits<int>::min();
    unordered_map<int, vector<RawData>> rawHits;
    unordered_set<int> warnedDetectorIDs;
    const Long64_t entries = dutTree.GetEntries();
    TimingProgress progress("DUT timing", entries, fitCount);
    for (Long64_t entry = 0; entry < entries; ++entry) {
        progress.Update(entry);
        dutTree.GetEntry(entry);

        vector<const ChannelHit*> matchedHits;
        bool hasMatchedCluster = false;
        const auto appendMatch =
            [&](const Cluster* cluster,
                const vector<ChannelHit>* hits) {
                if (!cluster || !hits || hits->empty()) return;
                hasMatchedCluster = true;
                for (const auto& hit : *hits)
                    matchedHits.push_back(&hit);
            };
        if (schema == DUTTreeSchema::Planar) {
            if (hitFlag & 1) appendMatch(clusterX, channelHitsX);
            if (hitFlag & 2) appendMatch(clusterY, channelHitsY);
        } else if (hitFlag != 0) {
            appendMatch(selectedCluster, selectedChannelHits);
        }
        const bool hasMatchedHit =
            hasMatchedCluster && !matchedHits.empty();

        const auto detector =
            DetectorFactory::GetInstance().GetDetector(dutID);
        if (!detector || !detector->isDUT() ||
            (!detector->GetPlanarConfig() &&
             !detector->GetPlanarPadConfig())) {
            if (warnedDetectorIDs.insert(dutID).second) {
                cerr << "[TimeResolution] skipping DUT tree entries for "
                     << dutID
                     << ": detector is missing, not a DUT, or unsupported\n";
            }
            continue;
        }

        const auto reference = trackReferences.find(eventID);
        if (reference == trackReferences.end()) continue;
        const bool insideActiveArea = IsInsideDUTActiveArea(
            detector, predictedX, predictedY);
        const TrackKey trackCase{
            eventID, reference->second.trackIndex};
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
        if (rawDetector == rawHits.end()) continue;
        vector<double> stripTimes;
        vector<double> stripTimeErrors;
        vector<double> stripAmplitudes;
        vector<double> stripRiseTimes;
        vector<double> stripWidths;
        vector<int> stripIDs;
        vector<int> stripPlanes;
        size_t validDUTChannels = 0;
        unordered_set<int> fittedRawIndices;
        for (const ChannelHit* selectedHit : matchedHits) {
            if (!selectedHit) continue;
            const int rawIndex = selectedHit->rawIndices;
            if (rawIndex < 0 ||
                rawIndex >= static_cast<int>(rawDetector->second.size()) ||
                !fittedRawIndices.insert(rawIndex).second) {
                continue;
            }
            const ChannelHit fitted =
                timingFitter.ProcessHit(rawDetector->second[rawIndex]);
            ++fitCount;
            if (fitted.isValid && isfinite(fitted.time)) {
                stripTimes.push_back(fitted.time);
                stripTimeErrors.push_back(fitted.timeError);
                stripAmplitudes.push_back(fitted.amp);
                stripRiseTimes.push_back(fitted.riseTime);
                stripWidths.push_back(fitted.width);
                stripIDs.push_back(fitted.id0);
                stripPlanes.push_back(fitted.type);
                ++validDUTChannels;
            }
        }
        if (validDUTChannels == 0) continue;
        const double dutTime = CalculateConfiguredTime(
            stripTimes, stripTimeErrors, stripAmplitudes,
            stripPlanes, timeMethod, amplitudeWeightPower);
        if (!isfinite(dutTime)) continue;
        DUTTimingSample sample;
        sample.eventID = eventID;
        sample.rawEventID = reference->second.rawEventID;
        sample.trackIndex = reference->second.trackIndex;
        sample.detectorID = dutID;
        sample.dutTime = dutTime;
        sample.residual = dutTime - reference->second.time;
        sample.insideActiveArea = insideActiveArea;
        sample.stripTimes = move(stripTimes);
        sample.stripTimeErrors = move(stripTimeErrors);
        sample.stripAmplitudes = move(stripAmplitudes);
        sample.stripRiseTimes = move(stripRiseTimes);
        sample.stripWidths = move(stripWidths);
        sample.stripIDs = move(stripIDs);
        sample.stripPlanes = move(stripPlanes);
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
    gStyle->SetOptStat(1);
    gStyle->SetOptFit(1);
    TH1D histogram(name.c_str(), title.c_str(), bins, low, high);
    histogram.SetStats(true);
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
            times[i] = detector->second.clusterTime;
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

void WriteTimeAmplitudeRelation(
    const vector<double>& amplitudes,
    const vector<double>& residuals,
    TDirectory& directory,
    const string& name,
    const string& title,
    int bins) {
    vector<double> finiteAmplitudes;
    vector<double> finiteResiduals;
    const size_t entries = min(amplitudes.size(), residuals.size());
    finiteAmplitudes.reserve(entries);
    finiteResiduals.reserve(entries);
    for (size_t index = 0; index < entries; ++index) {
        if (!isfinite(amplitudes[index]) ||
            !isfinite(residuals[index]))
            continue;
        finiteAmplitudes.push_back(amplitudes[index]);
        finiteResiduals.push_back(residuals[index]);
    }
    if (finiteAmplitudes.size() < 3) return;

    auto [amplitudeLow, amplitudeHigh] =
        RobustRange(finiteAmplitudes);
    auto [residualLow, residualHigh] =
        RobustRange(finiteResiduals);
    const double amplitudePadding =
        0.05 * max(1.0e-6, amplitudeHigh - amplitudeLow);
    const double residualPadding =
        0.05 * max(1.0e-6, residualHigh - residualLow);
    amplitudeLow -= amplitudePadding;
    amplitudeHigh += amplitudePadding;
    residualLow -= residualPadding;
    residualHigh += residualPadding;
    const int plotBins = max(20, min(bins, 100));

    TDirectory::TContext context(&directory);
    TH2D histogram(
        ("h_" + name).c_str(), title.c_str(),
        plotBins, amplitudeLow, amplitudeHigh,
        plotBins, residualLow, residualHigh);
    TProfile profile(
        ("p_" + name).c_str(),
        (title + " profile").c_str(),
        plotBins, amplitudeLow, amplitudeHigh,
        residualLow, residualHigh);
    histogram.SetStats(false);
    for (size_t index = 0; index < finiteAmplitudes.size(); ++index) {
        histogram.Fill(
            finiteAmplitudes[index], finiteResiduals[index]);
        profile.Fill(
            finiteAmplitudes[index], finiteResiduals[index]);
    }
    profile.SetLineColor(kRed + 1);
    profile.SetMarkerColor(kRed + 1);
    profile.SetMarkerStyle(20);

    TCanvas canvas(
        ("c_" + name).c_str(),
        title.c_str(), 1000, 750);
    canvas.SetRightMargin(0.14);
    histogram.Draw("COLZ");
    profile.Draw("SAME");
    histogram.Write();
    profile.Write();
    canvas.Write();
}

TimeWalkCorrection FitTimeWalkCorrection(
    const vector<double>& amplitudes,
    const vector<double>& residuals) {
    TimeWalkCorrection result;
    vector<double> finiteAmplitudes;
    vector<double> finiteResiduals;
    const size_t entries = min(amplitudes.size(), residuals.size());
    finiteAmplitudes.reserve(entries);
    finiteResiduals.reserve(entries);
    for (size_t index = 0; index < entries; ++index) {
        if (!isfinite(amplitudes[index]) ||
            !isfinite(residuals[index]))
            continue;
        finiteAmplitudes.push_back(amplitudes[index]);
        finiteResiduals.push_back(residuals[index]);
    }
    result.entries = finiteAmplitudes.size();
    if (finiteAmplitudes.size() < 20) return result;

    const auto [amplitudeLow, amplitudeHigh] =
        RobustRange(finiteAmplitudes);
    const auto [residualLow, residualHigh] =
        RobustRange(finiteResiduals);
    if (!(amplitudeHigh > amplitudeLow) ||
        !(residualHigh > residualLow))
        return result;

    static unsigned long fitCounter = 0;
    const string suffix = to_string(fitCounter++);
    TProfile profile(
        ("time_walk_profile_" + suffix).c_str(), "",
        100, amplitudeLow, amplitudeHigh,
        residualLow, residualHigh);
    profile.SetDirectory(nullptr);
    for (size_t index = 0; index < finiteAmplitudes.size(); ++index)
        profile.Fill(
            finiteAmplitudes[index], finiteResiduals[index]);

    TF1 fit(
        ("time_walk_fit_" + suffix).c_str(), "pol3",
        amplitudeLow, amplitudeHigh);
    if (profile.Fit(&fit, "QNR") != 0) return result;

    double fittedMean = 0.0;
    for (double amplitude : finiteAmplitudes)
        fittedMean += fit.Eval(
            clamp(amplitude, amplitudeLow, amplitudeHigh));
    fittedMean /= static_cast<double>(finiteAmplitudes.size());

    result.valid = true;
    result.amplitudeLow = amplitudeLow;
    result.amplitudeHigh = amplitudeHigh;
    for (size_t parameter = 0; parameter < result.parameter.size();
         ++parameter) {
        result.parameter[parameter] =
            fit.GetParameter(static_cast<int>(parameter));
    }
    result.parameter[0] -= fittedMean;
    return result;
}

double CorrectedClusterTime(
    const vector<double>& times,
    const vector<double>& errors,
    const vector<double>& amplitudes,
    const vector<int>& planes,
    const TimeWalkCorrection& correction,
    const string& timeMethod,
    double amplitudeWeightPower) {
    if (!correction.valid) return numeric_limits<double>::quiet_NaN();
    vector<double> correctedTimes;
    vector<double> correctedErrors;
    vector<double> correctedAmplitudes;
    vector<int> correctedPlanes;
    const size_t entries = min(
        {times.size(), amplitudes.size(), planes.size()});
    correctedTimes.reserve(entries);
    correctedErrors.reserve(entries);
    correctedAmplitudes.reserve(entries);
    correctedPlanes.reserve(entries);
    for (size_t index = 0; index < entries; ++index) {
        const double timeCorrection =
            correction.Evaluate(amplitudes[index]);
        if (!isfinite(times[index]) || !isfinite(timeCorrection)) continue;
        correctedTimes.push_back(times[index] - timeCorrection);
        correctedErrors.push_back(
            index < errors.size()
                ? errors[index]
                : numeric_limits<double>::quiet_NaN());
        correctedAmplitudes.push_back(amplitudes[index]);
        correctedPlanes.push_back(planes[index]);
    }
    return CalculateConfiguredTime(
        correctedTimes, correctedErrors, correctedAmplitudes,
        correctedPlanes, timeMethod, amplitudeWeightPower);
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
        "Amplitude-corrected DUT cluster time relative to external T0;"
        "t_{DUT}-T0 [ns];Entries",
        histogramBins);
    TDirectory::TContext context(&directory);
    auto* residualHistogram =
        dynamic_cast<TH1D*>(directory.Get("hTimeResidual"));
    if (residualHistogram) {
        gStyle->SetOptStat(1);
        gStyle->SetOptFit(1);
        residualHistogram->SetStats(true);
        TCanvas residualCanvas(
            "cTimeResidualFit",
            "DUT timing residual, Gaussian fit, and selected window",
            1000, 700);
        residualHistogram->Draw("HIST");
        TF1* gaussian = residualHistogram->GetFunction(
            "hTimeResidual_gaus");
        if (gaussian) {
            gaussian->SetLineColor(kGreen + 2);
            gaussian->SetLineWidth(3);
            gaussian->Draw("SAME");
        }
        const double yMaximum =
            1.08 * residualHistogram->GetMaximum();
        residualHistogram->SetMaximum(yMaximum);
        TLine windowStart(
            result.bestWindowStart, 0.0,
            result.bestWindowStart, yMaximum);
        TLine windowEnd(
            result.bestWindowEnd, 0.0,
            result.bestWindowEnd, yMaximum);
        for (TLine* line : {&windowStart, &windowEnd}) {
            line->SetLineColor(kRed + 1);
            line->SetLineStyle(2);
            line->SetLineWidth(3);
            line->Draw("SAME");
        }
        ostringstream annotation;
        annotation << fixed << setprecision(2)
                   << "selected " << windowWidthNs
                   << " ns window: ["
                   << result.bestWindowStart << ", "
                   << result.bestWindowEnd << "] ns";
        TLatex label;
        label.SetNDC();
        label.SetTextColor(kRed + 1);
        label.SetTextSize(0.035);
        label.DrawLatex(0.14, 0.86, annotation.str().c_str());
        residualCanvas.Modified();
        residualCanvas.Update();
        residualCanvas.Write();
    }
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
    const map<TrackKey, uint64_t>& rawEventIDs,
    const OscilloscopeT0Result& oscilloscopeT0,
    const array<int, 3>& trackerIDs,
    const DUTTimingResult& dutTiming,
    const TrackTimeWeights& trackWeights,
    const string& timeMethod,
    double amplitudeWeightPower,
    double maximumStripAmplitude,
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

    Int_t stripTreeEventID = -1;
    ULong64_t stripTreeRawEventID = 0;
    Int_t stripTreeTrackIndex = -1;
    Int_t stripTreeDetectorID = -1;
    Bool_t stripTreeIsDUT = false;
    Int_t stripTreePlane = -1;
    Int_t stripTreeChannelID = -1;
    Int_t stripTreeClusterStripCount = 0;
    Double_t stripTreeTime = numeric_limits<double>::quiet_NaN();
    Double_t stripTreeAmplitude = numeric_limits<double>::quiet_NaN();
    Double_t stripTreeRiseTime = numeric_limits<double>::quiet_NaN();
    Double_t stripTreeTimeOverThreshold =
        numeric_limits<double>::quiet_NaN();
    Bool_t stripTreeHasExternalT0 = false;
    Double_t stripTreeExternalT0 =
        numeric_limits<double>::quiet_NaN();
    Double_t stripTreeTimeMinusExternalT0 =
        numeric_limits<double>::quiet_NaN();
    TTree stripTimingTree(
        "StripTimingTree",
        "Per-strip timing observables for independent T0 validation");
    stripTimingTree.Branch("eventID", &stripTreeEventID);
    stripTimingTree.Branch("rawEventID", &stripTreeRawEventID);
    stripTimingTree.Branch("trackIndex", &stripTreeTrackIndex);
    stripTimingTree.Branch("detectorID", &stripTreeDetectorID);
    stripTimingTree.Branch("isDUT", &stripTreeIsDUT);
    stripTimingTree.Branch("plane", &stripTreePlane);
    stripTimingTree.Branch("channelID", &stripTreeChannelID);
    stripTimingTree.Branch(
        "clusterStripCount", &stripTreeClusterStripCount);
    stripTimingTree.Branch("time", &stripTreeTime);
    stripTimingTree.Branch("amplitude", &stripTreeAmplitude);
    stripTimingTree.Branch("riseTime", &stripTreeRiseTime);
    stripTimingTree.Branch(
        "timeOverThreshold", &stripTreeTimeOverThreshold);
    stripTimingTree.Branch("hasExternalT0", &stripTreeHasExternalT0);
    stripTimingTree.Branch("externalT0", &stripTreeExternalT0);
    stripTimingTree.Branch(
        "timeMinusExternalT0", &stripTreeTimeMinusExternalT0);
    const auto fillStripTiming =
        [&](int eventID, uint64_t rawEventID, int trackIndex,
            int detectorID, bool isDUT, const auto& timing,
            bool hasExternalT0, double externalT0) {
            const size_t stripCount = min(
                {timing.stripTimes.size(), timing.stripAmplitudes.size(),
                 timing.stripIDs.size(), timing.stripPlanes.size()});
            stripTreeEventID = eventID;
            stripTreeRawEventID = rawEventID;
            stripTreeTrackIndex = trackIndex;
            stripTreeDetectorID = detectorID;
            stripTreeIsDUT = isDUT;
            stripTreeClusterStripCount =
                static_cast<Int_t>(stripCount);
            stripTreeHasExternalT0 = hasExternalT0;
            stripTreeExternalT0 =
                hasExternalT0
                    ? externalT0
                    : numeric_limits<double>::quiet_NaN();
            for (size_t strip = 0; strip < stripCount; ++strip) {
                stripTreePlane = timing.stripPlanes[strip];
                stripTreeChannelID = timing.stripIDs[strip];
                stripTreeTime = timing.stripTimes[strip];
                stripTreeAmplitude = timing.stripAmplitudes[strip];
                stripTreeRiseTime =
                    strip < timing.stripRiseTimes.size()
                        ? timing.stripRiseTimes[strip]
                        : numeric_limits<double>::quiet_NaN();
                stripTreeTimeOverThreshold =
                    strip < timing.stripWidths.size()
                        ? timing.stripWidths[strip]
                        : numeric_limits<double>::quiet_NaN();
                stripTreeTimeMinusExternalT0 =
                    hasExternalT0 && isfinite(stripTreeTime)
                        ? stripTreeTime - externalT0
                        : numeric_limits<double>::quiet_NaN();
                stripTimingTree.Fill();
            }
        };

    map<TrackKey, vector<const DUTTimingSample*>> dutSamplesByTrack;
    for (const auto& [dutID, samples] : dutTiming.samplesByDetector) {
        (void)dutID;
        for (const DUTTimingSample& sample : samples)
            dutSamplesByTrack[
                {sample.eventID, sample.trackIndex}].push_back(&sample);
    }
    const auto exceedsMaximumStripAmplitude =
        [&](const vector<double>& amplitudes) {
            return maximumStripAmplitude > 0.0 &&
                   any_of(
                       amplitudes.begin(), amplitudes.end(),
                       [&](double amplitude) {
                           return isfinite(amplitude) &&
                                  amplitude > maximumStripAmplitude;
                       });
        };

    map<int, vector<double>> rawDetectorTimes;
    map<int, T0ResidualSamples> externalResiduals;
    map<int, vector<double>> externalStripAmplitudes;
    map<int, vector<double>> externalStripResiduals;
    array<vector<double>, 3> trackerPairResiduals;
    vector<double> trackTimeSamples;
    T0ResidualSamples trackExternalResiduals;
    const array<pair<size_t, size_t>, 3> pairIndices = {
        pair<size_t, size_t>{0, 1}, {0, 2}, {1, 2}};

    for (const auto& [key, detectors] : trackTimes) {
        const auto rawEvent = rawEventIDs.find(key);
        if (rawEvent == rawEventIDs.end()) continue;
        const auto external =
            oscilloscopeT0.references.find(rawEvent->second);
        for (const auto& [detectorID, times] : detectors) {
            fillStripTiming(
                key.first, rawEvent->second, key.second, detectorID,
                false, times,
                external != oscilloscopeT0.references.end(),
                external != oscilloscopeT0.references.end()
                    ? external->second.time
                    : numeric_limits<double>::quiet_NaN());
            const double time = times.clusterTime;
            if (!isfinite(time)) continue;
            rawDetectorTimes[detectorID].push_back(time);
            if (external != oscilloscopeT0.references.end()) {
                auto& residualSamples =
                    externalResiduals[detectorID];
                const size_t previousSize =
                    residualSamples.residual.size();
                residualSamples.Add(
                    time - external->second.time,
                    external->second.resolution);
                if (residualSamples.residual.size() > previousSize) {
                    const size_t stripCount = min(
                        times.stripTimes.size(),
                        times.stripAmplitudes.size());
                    for (size_t strip = 0; strip < stripCount; ++strip) {
                        if (!isfinite(times.stripTimes[strip]) ||
                            !isfinite(times.stripAmplitudes[strip]))
                            continue;
                        externalStripAmplitudes[detectorID].push_back(
                            times.stripAmplitudes[strip]);
                        externalStripResiduals[detectorID].push_back(
                            times.stripTimes[strip] -
                            external->second.time);
                    }
                }
            }
        }

        array<double, 3> selectedTimes{};
        bool complete = true;
        for (size_t index = 0; index < trackerIDs.size(); ++index) {
            const auto detector = detectors.find(trackerIDs[index]);
            if (detector == detectors.end() ||
                !isfinite(detector->second.clusterTime)) {
                complete = false;
                break;
            }
            selectedTimes[index] = detector->second.clusterTime;
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
            trackExternalResiduals.Add(
                trackTime - external->second.time,
                external->second.resolution);

        treeRawEventID = rawEvent->second;
        treeTrackIndex = key.second;
        treeEventID = key.first;
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
    map<int, T0ResidualSamples> dutExternalResiduals;
    map<int, T0ResidualSamples> dutEfficiencyExternalResiduals;
    map<int, vector<double>> dutExternalStripAmplitudes;
    map<int, vector<double>> dutExternalStripResiduals;
    for (const auto& [dutID, samples] : dutTiming.samplesByDetector) {
        for (const DUTTimingSample& sample : samples) {
            if (isfinite(sample.dutTime))
                dutTimes[dutID].push_back(sample.dutTime);
            if (isfinite(sample.residual))
                dutMinusTrack[dutID].push_back(sample.residual);
            const auto external =
                oscilloscopeT0.references.find(sample.rawEventID);
            fillStripTiming(
                sample.eventID, sample.rawEventID, sample.trackIndex,
                dutID, true, sample,
                external != oscilloscopeT0.references.end(),
                external != oscilloscopeT0.references.end()
                    ? external->second.time
                    : numeric_limits<double>::quiet_NaN());
            if (external == oscilloscopeT0.references.end() ||
                !isfinite(sample.dutTime))
                continue;
            if (exceedsMaximumStripAmplitude(sample.stripAmplitudes))
                continue;
            const double residual =
                sample.dutTime - external->second.time;
            auto& residualSamples =
                dutExternalResiduals[dutID];
            const size_t previousSize =
                residualSamples.residual.size();
            residualSamples.Add(
                residual, external->second.resolution);
            if (residualSamples.residual.size() > previousSize) {
                const size_t stripCount = min(
                    sample.stripTimes.size(),
                    sample.stripAmplitudes.size());
                for (size_t strip = 0; strip < stripCount; ++strip) {
                    if (!isfinite(sample.stripTimes[strip]) ||
                        !isfinite(sample.stripAmplitudes[strip]))
                        continue;
                    dutExternalStripAmplitudes[dutID].push_back(
                        sample.stripAmplitudes[strip]);
                    dutExternalStripResiduals[dutID].push_back(
                        sample.stripTimes[strip] -
                        external->second.time);
                }
            }
            if (sample.insideActiveArea)
                dutEfficiencyExternalResiduals[dutID].Add(
                    residual, external->second.resolution);
        }
    }

    map<int, TimeWalkCorrection> trackerTimeWalkCorrections;
    for (const auto& [detectorID, amplitudes] :
         externalStripAmplitudes) {
        const TimeWalkCorrection correction =
            FitTimeWalkCorrection(
                amplitudes, externalStripResiduals[detectorID]);
        if (correction.valid)
            trackerTimeWalkCorrections[detectorID] = correction;
    }
    map<int, TimeWalkCorrection> dutTimeWalkCorrections;
    for (const auto& [dutID, amplitudes] :
         dutExternalStripAmplitudes) {
        const TimeWalkCorrection correction =
            FitTimeWalkCorrection(
                amplitudes, dutExternalStripResiduals[dutID]);
        if (correction.valid)
            dutTimeWalkCorrections[dutID] = correction;
    }

    map<int, T0ResidualSamples> correctedExternalResiduals;
    T0ResidualSamples correctedTrackExternalResiduals;
    for (const auto& [key, detectors] : trackTimes) {
        const auto rawEvent = rawEventIDs.find(key);
        if (rawEvent == rawEventIDs.end()) continue;
        const auto external =
            oscilloscopeT0.references.find(rawEvent->second);
        if (external == oscilloscopeT0.references.end()) continue;
        map<int, double> correctedTimes;
        for (const auto& [detectorID, times] : detectors) {
            const auto correction =
                trackerTimeWalkCorrections.find(detectorID);
            if (correction == trackerTimeWalkCorrections.end()) continue;
            const double correctedTime = CorrectedClusterTime(
                times.stripTimes, times.stripTimeErrors,
                times.stripAmplitudes, times.stripPlanes,
                correction->second, timeMethod,
                amplitudeWeightPower);
            if (!isfinite(correctedTime)) continue;
            correctedTimes[detectorID] = correctedTime;
            correctedExternalResiduals[detectorID].Add(
                correctedTime - external->second.time,
                external->second.resolution);
        }
        array<double, 3> selectedTimes{};
        bool complete = true;
        for (size_t tracker = 0; tracker < trackerIDs.size(); ++tracker) {
            const auto time = correctedTimes.find(trackerIDs[tracker]);
            if (time == correctedTimes.end()) {
                complete = false;
                break;
            }
            selectedTimes[tracker] = time->second;
        }
        if (!complete) continue;
        const double correctedTrackTime = inner_product(
            selectedTimes.begin(), selectedTimes.end(),
            trackWeights.value.begin(), 0.0);
        correctedTrackExternalResiduals.Add(
            correctedTrackTime - external->second.time,
            external->second.resolution);
    }

    map<int, T0ResidualSamples> correctedDUTExternalResiduals;
    map<int, T0ResidualSamples> correctedDUTEfficiencyResiduals;
    for (const auto& [dutID, samples] : dutTiming.samplesByDetector) {
        const auto correction = dutTimeWalkCorrections.find(dutID);
        if (correction == dutTimeWalkCorrections.end()) continue;
        for (const DUTTimingSample& sample : samples) {
            const auto external =
                oscilloscopeT0.references.find(sample.rawEventID);
            if (external == oscilloscopeT0.references.end()) continue;
            if (exceedsMaximumStripAmplitude(sample.stripAmplitudes))
                continue;
            const double correctedTime = CorrectedClusterTime(
                sample.stripTimes, sample.stripTimeErrors,
                sample.stripAmplitudes, sample.stripPlanes,
                correction->second, timeMethod,
                amplitudeWeightPower);
            if (!isfinite(correctedTime)) continue;
            const double residual =
                correctedTime - external->second.time;
            correctedDUTExternalResiduals[dutID].Add(
                residual, external->second.resolution);
            if (sample.insideActiveArea)
                correctedDUTEfficiencyResiduals[dutID].Add(
                    residual, external->second.resolution);
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
    map<int, double> externalDetectorResolutions;
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
            const auto corrected =
                correctedExternalResiduals.find(detectorID);
            const T0ResidualSamples& fittedSamples =
                corrected != correctedExternalResiduals.end() &&
                        !corrected->second.residual.empty()
                    ? corrected->second
                    : samples;
            const FitResult fit = WriteAndFit(
                fittedSamples.residual, *detectorDirectory,
                "hClusterTimeMinusT0",
                label +
                    " amplitude-corrected cluster time relative to external T0;"
                        "t-T0 [ns];Entries",
                histogramBins);
            externalDetectorFits[detectorID] = fit;
            externalDetectorResolutions[detectorID] =
                SubtractReferenceResolution(fit, fittedSamples);
            TDirectory* timeAmplitudeDirectory =
                detectorDirectory->mkdir("TimeAmplitude");
            WriteTimeAmplitudeRelation(
                externalStripAmplitudes[detectorID],
                externalStripResiduals[detectorID],
                *timeAmplitudeDirectory,
                "strip_time_minus_t0_vs_amplitude",
                label +
                    " individual-strip time relative to external T0 vs amplitude;"
                    "Fitted strip amplitude [ADC];"
                    "t_{strip}-T0 [ns]",
                histogramBins);
        }
        TDirectory* trackDirectory =
            externalDirectory->mkdir("TrackTime");
        const T0ResidualSamples& fittedTrackSamples =
            !correctedTrackExternalResiduals.residual.empty()
                ? correctedTrackExternalResiduals
                : trackExternalResiduals;
        externalTrackFit = WriteAndFit(
            fittedTrackSamples.residual, *trackDirectory,
            "hTrackTimeMinusT0",
            "Amplitude-corrected track time relative to external T0;"
            "t_{track}-T0 [ns];Entries",
            histogramBins);
        externalTrackResolution = SubtractReferenceResolution(
            externalTrackFit, fittedTrackSamples);

        for (const auto& [dutID, samples] : dutExternalResiduals) {
            TDirectory* detectorDirectory =
                externalDirectory->mkdir(
                    ("DUT" + to_string(dutID)).c_str());
            const auto corrected =
                correctedDUTExternalResiduals.find(dutID);
            const T0ResidualSamples& fittedSamples =
                corrected != correctedDUTExternalResiduals.end() &&
                        !corrected->second.residual.empty()
                    ? corrected->second
                    : samples;
            const FitResult fit = WriteAndFit(
                fittedSamples.residual, *detectorDirectory,
                "hClusterTimeMinusT0",
                "Amplitude-corrected DUT cluster time relative to external T0;"
                "t_{DUT}-T0 [ns];Entries",
                histogramBins);
            externalDUTFits[dutID] = fit;
            externalDUTResolutions[dutID] =
                SubtractReferenceResolution(fit, fittedSamples);
            TDirectory* timeAmplitudeDirectory =
                detectorDirectory->mkdir("TimeAmplitude");
            WriteTimeAmplitudeRelation(
                dutExternalStripAmplitudes[dutID],
                dutExternalStripResiduals[dutID],
                *timeAmplitudeDirectory,
                "strip_time_minus_t0_vs_amplitude",
                "DUT" + to_string(dutID) +
                    " individual-strip time relative to external T0 vs amplitude;"
                    "Fitted strip amplitude [ADC];"
                    "t_{strip}-T0 [ns]",
                histogramBins);
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
                        const auto rawEvent = rawEventIDs.find(key);
                        return rawEvent != rawEventIDs.end() &&
                               oscilloscopeT0.references.count(
                                   rawEvent->second) != 0;
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
                            const auto rawEvent = rawEventIDs.find(key);
                            return rawEvent != rawEventIDs.end() &&
                                   oscilloscopeT0.references.count(
                                       rawEvent->second) != 0;
                        }));
            }
            TDirectory* detectorDirectory =
                efficiencyRoot->mkdir(
                    ("DUT" + to_string(dutID)).c_str());
            const auto corrected =
                correctedDUTEfficiencyResiduals.find(dutID);
            const T0ResidualSamples& efficiencySamples =
                corrected != correctedDUTEfficiencyResiduals.end() &&
                        !corrected->second.residual.empty()
                    ? corrected->second
                    : residuals;
            const TimingEfficiencyResult result =
                WriteCompactTimingEfficiency(
                    efficiencySamples.residual, denominator,
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
    Bool_t resultAmplitudeCorrected = false;
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
    resultsTree.Branch(
        "amplitudeCorrected", &resultAmplitudeCorrected);
    const auto fillResult = [&](const string& analysis,
                                const string& object, int detectorA,
                                int detectorB, const FitResult& fit,
                                double resolution,
                                bool amplitudeCorrected = false) {
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
        resultAmplitudeCorrected = amplitudeCorrected;
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
                   fit, externalDetectorResolutions[detectorID],
                   trackerTimeWalkCorrections.count(detectorID) != 0);
    if (isfinite(externalTrackResolution)) {
        fillResult("with_external_t0", "track_time", -1, -1,
                   externalTrackFit,
                   externalTrackResolution,
                   !correctedTrackExternalResiduals.residual.empty());
    }
    for (const auto& [dutID, fit] : externalDUTFits)
        fillResult("with_external_t0", "dut", dutID, -1, fit,
                   externalDUTResolutions[dutID],
                   dutTimeWalkCorrections.count(dutID) != 0);
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
        resultAmplitudeCorrected =
            correctedDUTEfficiencyResiduals.count(dutID) != 0 &&
            !correctedDUTEfficiencyResiduals.at(dutID).residual.empty();
        resultsTree.Fill();
    }

    string correctionSource;
    Int_t correctionDetectorID = -1;
    Long64_t correctionEntries = 0;
    Double_t correctionAmplitudeLow =
        numeric_limits<double>::quiet_NaN();
    Double_t correctionAmplitudeHigh =
        numeric_limits<double>::quiet_NaN();
    array<Double_t, 4> correctionParameter{};
    TTree correctionTree(
        "TimeAmplitudeCorrections",
        "Applied cubic single-strip time-amplitude correction parameters");
    correctionTree.Branch("source", &correctionSource);
    correctionTree.Branch("detectorID", &correctionDetectorID);
    correctionTree.Branch("entries", &correctionEntries);
    correctionTree.Branch("amplitudeLow", &correctionAmplitudeLow);
    correctionTree.Branch("amplitudeHigh", &correctionAmplitudeHigh);
    correctionTree.Branch(
        "parameter", correctionParameter.data(), "parameter[4]/D");
    const auto fillCorrection = [&](const string& source, int detectorID,
                                    const TimeWalkCorrection& correction) {
        correctionSource = source;
        correctionDetectorID = detectorID;
        correctionEntries =
            static_cast<Long64_t>(correction.entries);
        correctionAmplitudeLow = correction.amplitudeLow;
        correctionAmplitudeHigh = correction.amplitudeHigh;
        copy(correction.parameter.begin(), correction.parameter.end(),
             correctionParameter.begin());
        correctionTree.Fill();
    };
    for (const auto& [detectorID, correction] :
         trackerTimeWalkCorrections)
        fillCorrection("tracker", detectorID, correction);
    for (const auto& [detectorID, correction] : dutTimeWalkCorrections)
        fillCorrection("dut", detectorID, correction);

    output->cd();
    timingTree.Write();
    stripTimingTree.Write();
    resultsTree.Write();
    correctionTree.Write();
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
        row("Associated tracks",
            count(trackExternalResiduals.residual.size()));
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
         << setw(23) << "Self-calibration"
         << setw(24) << "External T0 (A-corr.)"
         << '\n';
    const auto resolutionRow = [&](const string& label,
                                   const string& self,
                                   const string& external) {
        cout << "  " << left << setw(31) << label << right << setw(23)
             << self << setw(24) << external << '\n';
    };
    for (size_t index = 0; index < trackerIDs.size(); ++index) {
        const auto detector =
            DetectorFactory::GetInstance().GetDetector(trackerIDs[index]);
        resolutionRow(
            detector ? detector->GetName()
                     : "Tracker" + to_string(index + 1),
            fixedValue(sqrt(trackWeights.detectorVariance[index]), 2,
                       " ns"),
            externalDetectorResolutions.count(trackerIDs[index])
                ? fixedValue(
                      externalDetectorResolutions.at(trackerIDs[index]),
                      2, " ns")
                : "—");
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
    if (config.contains("timeMethod"))
        m_timeMethod = CanonicalTimeMethod(
            config["timeMethod"].get<string>());
    m_amplitudeWeightPower = config.value(
        "amplitudeWeightPower", m_amplitudeWeightPower);
    m_maximumStripAmplitude = config.value(
        "maximumStripAmplitude", m_maximumStripAmplitude);
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
    if (m_analyzeDUTTiming)
        cout << "  DUT timing file=" << m_dutFile << '\n';
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
           !m_timeMethod.empty() &&
           isfinite(m_amplitudeWeightPower) &&
           m_amplitudeWeightPower >= 0.0 &&
           isfinite(m_maximumStripAmplitude) &&
           m_maximumStripAmplitude != 0.0 &&
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
            m_timeMethod, m_amplitudeWeightPower,
            waveformFits, reconstruction))
        return false;
    if (Terminal::Verbose()) {
        Terminal::Detail(
            "loaded " + Terminal::Count(reconstruction.trackTimes.size()) +
            " tracks · " +
            Terminal::Count(reconstruction.wantedEventIDs.size()) +
            " raw event IDs");
        if (!reconstruction.ambiguousRawEventIDs.empty())
            Terminal::Detail(
                Terminal::Warning(
                    Terminal::Count(
                        reconstruction.ambiguousRawEventIDs.size()) +
                    " repeated raw event IDs excluded from external T0 matching"));
    }
    const TrackTimeWeights trackWeights = CalculateTrackTimeWeights(
        reconstruction.trackTimes, trackerIDs);
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
        auto* dutTree = dynamic_cast<TTree*>(dutFile->Get("DUTTree"));
        if (!dutTree)
            dutTree =
                dynamic_cast<TTree*>(dutFile->Get("PadDUTTree"));
        if (!dutTree) {
            cerr << "[TimeResolution] neither DUTTree nor PadDUTTree is "
                    "present in "
                 << dutPath << '\n';
            return false;
        }
        const auto references = BuildTrackTimeReferences(
            reconstruction.trackTimes, reconstruction.rawEventIDs, trackerIDs,
            trackWeights);
        dutTiming = LoadDUTTiming(
            *dutTree, *parser, references, m_timingWaveformConfig,
            m_timeMethod, m_amplitudeWeightPower, waveformFits);
    }

    const auto& trackTimes = reconstruction.trackTimes;

    const bool wroteOutput =
        WriteCompactTimingOutput(
            outputPath, trackTimes, reconstruction.rawEventIDs,
            oscilloscopeT0, trackerIDs, dutTiming, trackWeights,
            m_timeMethod, m_amplitudeWeightPower,
            m_maximumStripAmplitude,
            m_timingEfficiencyWindowNs,
            m_timingEfficiencyStepNs, m_histogramBins, analysisStarted);
    if (Terminal::Verbose()) {
        Terminal::Detail(Terminal::Muted(
            Terminal::Count(waveformFits) + " waveform fits"));
    }
    return wroteOutput;
}

REGISTER_SCRIPT("TimeResolution", TimeResolutionScript);
