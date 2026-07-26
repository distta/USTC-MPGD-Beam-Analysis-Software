#include "Script/TimeResolutionScript.h"

#include "Algorithm/Analyzer/WaveformProcessor.h"
#include "Detector/DetectorFactory.h"
#include "Event/DataModel.h"
#include "Script/Base/RawDataParser.h"
#include "Script/Base/ScriptFactory.h"

#include <TFile.h>
#include <TDirectory.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TTree.h>

#include <Eigen/Dense>

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
};

struct DUTTimingResult {
    map<int, vector<DUTTimingSample>> samplesByDetector;
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

struct TimeMeasurement {
    double time = numeric_limits<double>::quiet_NaN();
    double error = numeric_limits<double>::quiet_NaN();
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
                                            WaveformProcessor& timingFitter) {
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
                const ChannelHit fitted = timingFitter.ProcessWaveform(detectorRawData[rawIndex]);
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

    WaveformProcessor timingFitter;
    timingFitter.LoadConfig(timingWaveformConfig);
    int loadedEventID = numeric_limits<int>::min();
    unordered_map<int, vector<RawData>> rawHits;
    for (Long64_t entry = 0; entry < validation.GetEntries(); ++entry) {
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
        result.eventIDs[{rawEventID, trackIndex}] = eventID;
        result.wantedEventIDs.insert(rawEventID);
    }
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
    const json& timingWaveformConfig) {
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

    WaveformProcessor timingFitter;
    timingFitter.LoadConfig(timingWaveformConfig);
    int loadedEventID = numeric_limits<int>::min();
    unordered_map<int, vector<RawData>> rawHits;
    for (Long64_t entry = 0; entry < dutTree.GetEntries(); ++entry) {
        dutTree.GetEntry(entry);
        if (hitFlag == 0 || !selectedCluster || !selectedChannelHits ||
            selectedChannelHits->empty()) {
            continue;
        }
        const auto reference = trackReferences.find(eventID);
        if (reference == trackReferences.end()) {
            ++result.unmatchedTrackTimes;
            continue;
        }
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
                timingFitter.ProcessWaveform(rawDetector->second[rawIndex]);
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
        result.samplesByDetector[dutID].push_back(move(sample));
    }
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

int CrossValidationFold(uint64_t eventID) {
    return static_cast<int>(
        (eventID * 11400714819323198485ULL) >> 63);
}

struct QuadraticTimeWalkModel {
    double amplitudeMean = numeric_limits<double>::quiet_NaN();
    double amplitudeScale = numeric_limits<double>::quiet_NaN();
    Eigen::Vector3d coefficients =
        Eigen::Vector3d::Constant(numeric_limits<double>::quiet_NaN());
    bool valid = false;

    double Evaluate(double amplitude) const {
        if (!valid || !isfinite(amplitude)) return 0.0;
        const double normalized =
            (amplitude - amplitudeMean) / amplitudeScale;
        return coefficients[0] + coefficients[1] * normalized +
               coefficients[2] * normalized * normalized;
    }
};

struct DUTTimeWalkResult {
    bool accepted = false;
    array<QuadraticTimeWalkModel, 2> models;
    vector<double> correctedResiduals;
    FitResult rawFit;
    FitResult correctedFit;
    array<FitResult, 2> rawFoldFits;
    array<FitResult, 2> correctedFoldFits;
};

QuadraticTimeWalkModel FitQuadraticTimeWalk(
    const vector<DUTTimingSample>& samples, int heldOutFold) {
    QuadraticTimeWalkModel model;
    double amplitudeSum = 0.0;
    size_t count = 0;
    for (const DUTTimingSample& sample : samples) {
        if (CrossValidationFold(sample.rawEventID) == heldOutFold ||
            !isfinite(sample.amplitude) || !isfinite(sample.residual))
            continue;
        amplitudeSum += sample.amplitude;
        ++count;
    }
    if (count < 100) return model;
    model.amplitudeMean = amplitudeSum / static_cast<double>(count);
    double amplitudeVariance = 0.0;
    for (const DUTTimingSample& sample : samples) {
        if (CrossValidationFold(sample.rawEventID) == heldOutFold ||
            !isfinite(sample.amplitude) || !isfinite(sample.residual))
            continue;
        const double delta = sample.amplitude - model.amplitudeMean;
        amplitudeVariance += delta * delta;
    }
    model.amplitudeScale =
        sqrt(amplitudeVariance / static_cast<double>(count - 1));
    if (!isfinite(model.amplitudeScale) || model.amplitudeScale <= 0.0)
        return model;

    Eigen::Matrix3d normal = Eigen::Matrix3d::Zero();
    Eigen::Vector3d rhs = Eigen::Vector3d::Zero();
    for (const DUTTimingSample& sample : samples) {
        if (CrossValidationFold(sample.rawEventID) == heldOutFold ||
            !isfinite(sample.amplitude) || !isfinite(sample.residual))
            continue;
        const double normalized =
            (sample.amplitude - model.amplitudeMean) / model.amplitudeScale;
        const Eigen::Vector3d features(
            1.0, normalized, normalized * normalized);
        normal.noalias() += features * features.transpose();
        rhs.noalias() += features * sample.residual;
    }
    model.coefficients = normal.ldlt().solve(rhs);
    model.valid = model.coefficients.allFinite();
    return model;
}

DUTTimeWalkResult EvaluateDUTTimeWalk(
    const vector<DUTTimingSample>& samples) {
    DUTTimeWalkResult result;
    if (samples.size() < 1000) return result;
    result.models[0] = FitQuadraticTimeWalk(samples, 0);
    result.models[1] = FitQuadraticTimeWalk(samples, 1);
    if (!result.models[0].valid || !result.models[1].valid) return result;

    vector<double> rawResiduals;
    array<vector<double>, 2> rawByFold;
    array<vector<double>, 2> correctedByFold;
    rawResiduals.reserve(samples.size());
    result.correctedResiduals.reserve(samples.size());
    for (const DUTTimingSample& sample : samples) {
        const int fold = CrossValidationFold(sample.rawEventID);
        const double corrected =
            sample.residual - result.models[fold].Evaluate(sample.amplitude);
        rawResiduals.push_back(sample.residual);
        result.correctedResiduals.push_back(corrected);
        rawByFold[fold].push_back(sample.residual);
        correctedByFold[fold].push_back(corrected);
    }
    result.rawFit = FitGaussianWithoutWriting(rawResiduals);
    result.correctedFit =
        FitGaussianWithoutWriting(result.correctedResiduals);
    for (size_t fold = 0; fold < 2; ++fold) {
        result.rawFoldFits[fold] =
            FitGaussianWithoutWriting(rawByFold[fold]);
        result.correctedFoldFits[fold] =
            FitGaussianWithoutWriting(correctedByFold[fold]);
    }
    result.accepted =
        isfinite(result.rawFit.sigma) &&
        isfinite(result.correctedFit.sigma) &&
        result.correctedFit.sigma < 0.995 * result.rawFit.sigma;
    for (size_t fold = 0; fold < 2; ++fold) {
        result.accepted =
            result.accepted &&
            isfinite(result.rawFoldFits[fold].sigma) &&
            isfinite(result.correctedFoldFits[fold].sigma) &&
            result.correctedFoldFits[fold].sigma <
                result.rawFoldFits[fold].sigma;
    }
    if (!result.accepted) result.correctedResiduals.clear();
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
                       const array<int, 3>& trackerIDs,
                       const DUTTimingResult& dutTiming,
                       const TrackTimeWeights& trackWeights,
                       bool applyDUTTimeWalkCorrection,
                       int histogramBins) {
    unique_ptr<TFile> output(TFile::Open(outputPath.c_str(), "RECREATE"));
    if (!output || output->IsZombie()) {
        cerr << "[TimeResolution] cannot create " << outputPath << '\n';
        return false;
    }
    map<int, DUTTimeWalkResult> dutTimeWalk;
    if (applyDUTTimeWalkCorrection) {
        for (const auto& [detectorID, samples] : dutTiming.samplesByDetector)
            dutTimeWalk.emplace(detectorID, EvaluateDUTTimeWalk(samples));
    }

    TTree eventTree("EventTimes", "Tracker detector times");
    ULong64_t outRawEventID = 0;
    Int_t outTrackIndex = 0, outDetID = 0;
    array<Double_t, kEstimatorCount> estimatorValues{};
    eventTree.Branch("rawEventID", &outRawEventID);
    eventTree.Branch("trackIndex", &outTrackIndex);
    eventTree.Branch("detectorID", &outDetID);
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
    trackTimeTree.Branch("tracker1Weight", &tracker1Weight);
    trackTimeTree.Branch("tracker2Weight", &tracker2Weight);
    trackTimeTree.Branch("tracker3Weight", &tracker3Weight);

    map<pair<int, size_t>, vector<double>> detectorTimeSamples;
    map<tuple<int, int, size_t>, vector<double>> pairResiduals;
    vector<double> trackTimeSamples;
    array<vector<double>, 3> threeTrackerPairResiduals;
    size_t completeThreeTrackerTimes = 0;
    for (const auto& [key, detectors] : trackTimes) {
        for (const auto& [detectorID, times] : detectors) {
            outRawEventID = key.first;
            outTrackIndex = key.second;
            outDetID = detectorID;
            estimatorValues = times.value;
            eventTree.Fill();
            for (size_t estimator = 0; estimator < kEstimatorCount; ++estimator) {
                if (!isfinite(times.value[estimator])) continue;
                detectorTimeSamples[{detectorID, estimator}].push_back(times.value[estimator]);
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
            trackTimeTree.Fill();
            trackTimeSamples.push_back(trackTime);
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
    Double_t dutTimeWalkCorrection = 0.0;
    Double_t dutMinusTrackTimeCorrected =
        numeric_limits<double>::quiet_NaN();
    Bool_t dutTimeWalkApplied = false;
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
    dutTimeTree.Branch("timeWalkCorrection", &dutTimeWalkCorrection);
    dutTimeTree.Branch(
        "dutMinusTrackTimeCorrected", &dutMinusTrackTimeCorrected);
    dutTimeTree.Branch("timeWalkApplied", &dutTimeWalkApplied);
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
    for (const auto& [detectorID, samples] : dutTiming.samplesByDetector) {
        const auto correction = dutTimeWalk.find(detectorID);
        const bool correctionAccepted =
            correction != dutTimeWalk.end() && correction->second.accepted;
        size_t sampleIndex = 0;
        for (const DUTTimingSample& sample : samples) {
            dutEventID = sample.eventID;
            dutRawEventID = sample.rawEventID;
            dutTrackIndex = sample.trackIndex;
            dutDetectorID = detectorID;
            dutTime = sample.dutTime;
            dutTrackTime = sample.trackTime;
            dutMinusTrackTime = sample.residual;
            dutTimeWalkApplied = correctionAccepted;
            dutMinusTrackTimeCorrected =
                correctionAccepted
                    ? correction->second.correctedResiduals[sampleIndex]
                    : sample.residual;
            dutTimeWalkCorrection =
                sample.residual - dutMinusTrackTimeCorrected;
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
            ++sampleIndex;
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

    cout << "[TimeResolution] three-tracker XYMean: IDs=" << trackerIDs[0] << ','
         << trackerIDs[1] << ',' << trackerIDs[2]
         << ", complete tracks=" << completeThreeTrackerTimes;
    if (resolutionValid) {
        cout << ", resolutions=" << tracker1ResolutionNs << ','
             << tracker2ResolutionNs << ',' << tracker3ResolutionNs
             << " ns, weights=" << tracker1Weight << ','
             << tracker2Weight << ',' << tracker3Weight
             << ", weighted track resolution=" << trackResolutionNs << " ns\n";
    } else {
        cout << ", resolution invalid (missing pair fit or negative derived variance)\n";
    }

    if (!dutTiming.samplesByDetector.empty()) {
        TDirectory* dutDirectory = output->mkdir("DUTTimeResolution");
        for (const auto& [dutID, samples] : dutTiming.samplesByDetector) {
            vector<double> dutTimes;
            vector<double> referenceTimes;
            vector<double> residuals;
            vector<double> correctedResiduals;
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
            const auto correction = dutTimeWalk.find(dutID);
            const bool correctionAccepted =
                correction != dutTimeWalk.end() && correction->second.accepted;
            if (correctionAccepted)
                correctedResiduals =
                    correction->second.correctedResiduals;
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
            const FitResult correctedMeasured =
                correctionAccepted
                    ? WriteAndFit(
                          correctedResiduals, *detectorDirectory,
                          "dut_minus_track_time_timewalk_corrected",
                          "Cross-validated amplitude time-walk corrected DUT "
                          "residual;t_{DUT}-t_{track}-#hat{f}(A) [ns];Events",
                          histogramBins)
                    : measured;

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

            if (applyDUTTimeWalkCorrection &&
                correction != dutTimeWalk.end()) {
                const DUTTimeWalkResult& calibration = correction->second;
                TTree timeWalkSummary(
                    "DUTTimeWalkCalibration",
                    "Cross-validated quadratic amplitude time-walk calibration");
                Int_t accepted = calibration.accepted;
                Double_t rawSigma = calibration.rawFit.sigma;
                Double_t correctedSigma = calibration.correctedFit.sigma;
                Double_t rawFoldSigma[2] = {
                    calibration.rawFoldFits[0].sigma,
                    calibration.rawFoldFits[1].sigma};
                Double_t correctedFoldSigma[2] = {
                    calibration.correctedFoldFits[0].sigma,
                    calibration.correctedFoldFits[1].sigma};
                Double_t amplitudeMean[2] = {
                    calibration.models[0].amplitudeMean,
                    calibration.models[1].amplitudeMean};
                Double_t amplitudeScale[2] = {
                    calibration.models[0].amplitudeScale,
                    calibration.models[1].amplitudeScale};
                Double_t coefficients[2][3] = {};
                for (size_t fold = 0; fold < 2; ++fold)
                    for (size_t coefficient = 0; coefficient < 3;
                         ++coefficient)
                        coefficients[fold][coefficient] =
                            calibration.models[fold]
                                .coefficients[coefficient];
                timeWalkSummary.Branch("accepted", &accepted);
                timeWalkSummary.Branch("rawSigmaNs", &rawSigma);
                timeWalkSummary.Branch(
                    "correctedSigmaNs", &correctedSigma);
                timeWalkSummary.Branch(
                    "rawFoldSigmaNs", rawFoldSigma, "rawFoldSigmaNs[2]/D");
                timeWalkSummary.Branch(
                    "correctedFoldSigmaNs", correctedFoldSigma,
                    "correctedFoldSigmaNs[2]/D");
                timeWalkSummary.Branch(
                    "amplitudeMean", amplitudeMean, "amplitudeMean[2]/D");
                timeWalkSummary.Branch(
                    "amplitudeScale", amplitudeScale,
                    "amplitudeScale[2]/D");
                timeWalkSummary.Branch(
                    "coefficients", coefficients, "coefficients[2][3]/D");
                detectorDirectory->cd();
                timeWalkSummary.Fill();
                timeWalkSummary.Write();
            }

            TTree dutSummary(
                "DUTResolution",
                "DUT resolution after subtracting tracker reference in quadrature");
            Int_t valid = 0;
            Int_t summaryDUTID = dutID;
            Long64_t matchedEntries = measured.entries;
            Int_t timeWalkApplied = correctionAccepted;
            Double_t rawResidualMeanNs = measured.mean;
            Double_t rawMeasuredSigmaNs = measured.sigma;
            Double_t rawMeasuredSigmaErrorNs = measured.sigmaError;
            Double_t residualMeanNs = correctedMeasured.mean;
            Double_t measuredSigmaNs = correctedMeasured.sigma;
            Double_t measuredSigmaErrorNs = correctedMeasured.sigmaError;
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
            dutSummary.Branch("timeWalkApplied", &timeWalkApplied);
            dutSummary.Branch("rawResidualMeanNs", &rawResidualMeanNs);
            dutSummary.Branch("rawMeasuredSigmaNs", &rawMeasuredSigmaNs);
            dutSummary.Branch(
                "rawMeasuredSigmaErrorNs", &rawMeasuredSigmaErrorNs);
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

            cout << "[TimeResolution] DUT " << dutID
                 << ": matched=" << samples.size()
                 << ", raw measured sigma=" << rawMeasuredSigmaNs << " ns";
            if (applyDUTTimeWalkCorrection) {
                cout << ", time-walk correction="
                     << (correctionAccepted ? "accepted" : "rejected");
                if (correctionAccepted)
                    cout << ", corrected sigma=" << measuredSigmaNs << " ns";
            }
            if (valid) {
                cout << ", DUT resolution=" << dutResolutionNs << " ns\n";
            } else {
                cout << ", DUT resolution invalid\n";
            }
        }
    }

    output->cd();
    summary.Write();
    cout << "[TimeResolution] tracks=" << trackTimes.size()
         << ", output=" << outputPath << '\n';
    return !trackTimes.empty();
}

}  // namespace

void TimeResolutionScript::LoadConfig(const json& config) {
    m_trackFile = config.value("trackFile", m_trackFile);
    m_dutFile = config.value("dutFile", m_dutFile);
    m_outputFile = config.value("outputFile", m_outputFile);
    m_histogramBins = config.value("histogramBins", m_histogramBins);
    m_trackerIDs = config.value("trackerIDs", m_trackerIDs);
    m_analyzeDUTTiming =
        config.value("analyzeDUTTiming", m_analyzeDUTTiming);
    m_applyDUTTimeWalkCorrection = config.value(
        "applyDUTTimeWalkCorrection", m_applyDUTTimeWalkCorrection);
    if (config.contains("timingWaveform"))
        m_timingWaveformConfig = config["timingWaveform"];
}

void TimeResolutionScript::Print() const {
    cout << "TimeResolutionScript:\n";
    if (!m_trackerIDs.empty()) {
        cout << "  three-tracker IDs=";
        for (size_t i = 0; i < m_trackerIDs.size(); ++i)
            cout << (i == 0 ? "" : ",") << m_trackerIDs[i];
        cout << '\n';
    }
    cout << "  DUT timing=" << (m_analyzeDUTTiming ? "enabled" : "disabled")
         << '\n';
    cout << "  DUT time-walk correction="
         << (m_applyDUTTimeWalkCorrection ? "enabled" : "disabled") << '\n';
    cout << "  timing waveform=" << m_timingWaveformConfig.dump() << '\n';
}

bool TimeResolutionScript::Validate() const {
    return m_histogramBins > 0 &&
           (m_trackerIDs.empty() || m_trackerIDs.size() == 3) &&
           set<int>(m_trackerIDs.begin(), m_trackerIDs.end()).size() ==
               m_trackerIDs.size();
}

bool TimeResolutionScript::Execute() {
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
    if (!LoadReconstructionTimes(
            *validation, *parser, m_timingWaveformConfig, reconstruction))
        return false;
    cout << "[TimeResolution] loaded " << reconstruction.trackTimes.size() << " tracks and "
         << reconstruction.wantedEventIDs.size() << " raw event IDs\n";
    if (reconstruction.missingRawDetectors > 0 || reconstruction.emptyTimingDetectors > 0)
        cout << "[TimeResolution] Fit timing diagnostics: missing raw detector entries="
             << reconstruction.missingRawDetectors
             << ", empty detector times=" << reconstruction.emptyTimingDetectors << '\n';

    const TrackTimeWeights trackWeights =
        CalculateTrackTimeWeights(reconstruction.trackTimes, trackerIDs);
    if (!trackWeights.valid) {
        cerr << "[TimeResolution] cannot derive positive tracker variances "
                "for inverse-variance track-time weights\n";
        return false;
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
            *dutTree, *parser, references, m_timingWaveformConfig);
        cout << "[TimeResolution] DUT timing samples=";
        size_t totalDUTSamples = 0;
        for (const auto& [detectorID, samples] : dutTiming.samplesByDetector) {
            totalDUTSamples += samples.size();
            cout << " DUT" << detectorID << ':' << samples.size();
        }
        cout << ", total=" << totalDUTSamples
             << ", unmatched track references=" << dutTiming.unmatchedTrackTimes
             << ", invalid DUT times=" << dutTiming.invalidDUTTimes << '\n';
    }

    const auto& trackTimes = reconstruction.trackTimes;

    return WriteTimingOutput(outputPath, trackTimes, reconstruction.eventIDs,
                             trackerIDs, dutTiming, trackWeights,
                             m_applyDUTTimeWalkCorrection,
                             m_histogramBins);
}

REGISTER_SCRIPT("TimeResolution", TimeResolutionScript);
