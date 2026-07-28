#include "Algorithm/Oscilloscope/OscilloscopeDataProcessor.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
#include <map>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

namespace {

constexpr int kEventIDBits = 16;
constexpr int kFrameBits = kEventIDBits + 2;
constexpr double kEventBitPeriodSeconds = 25.0e-9;

struct Waveform {
    vector<double> time;
    vector<double> amplitude;
};

enum class TimingStatus {
    NoWindow,
    BelowAmplitude,
    NoCrossing,
    Valid,
};

struct TimingResult {
    double timeNs = numeric_limits<double>::quiet_NaN();
    double amplitude = numeric_limits<double>::quiet_NaN();
    double baseline = numeric_limits<double>::quiet_NaN();
    double threshold = numeric_limits<double>::quiet_NaN();
    TimingStatus status = TimingStatus::NoWindow;
};

struct EventIDResult {
    uint64_t eventID = 0;
    double frameStartNs = numeric_limits<double>::quiet_NaN();
    double threshold = numeric_limits<double>::quiet_NaN();
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

double InterpolatedCrossing(const vector<double>& time,
                            const vector<double>& amplitude,
                            size_t second, double threshold) {
    if (second == 0 || second >= time.size() ||
        second >= amplitude.size())
        return numeric_limits<double>::quiet_NaN();
    const double firstAmplitude = amplitude[second - 1];
    const double secondAmplitude = amplitude[second];
    const double difference = secondAmplitude - firstAmplitude;
    if (difference == 0.0) return time[second];
    const double fraction = (threshold - firstAmplitude) / difference;
    return time[second - 1] +
           fraction * (time[second] - time[second - 1]);
}

double InterpolatedCrossing(const Waveform& waveform, size_t second,
                            double threshold) {
    return InterpolatedCrossing(waveform.time, waveform.amplitude,
                                second, threshold);
}

vector<double> MedianFilter(const vector<double>& values, size_t count,
                            int windowSamples) {
    vector<double> filtered(count);
    vector<double> neighborhood;
    neighborhood.reserve(static_cast<size_t>(windowSamples));
    const size_t radius = static_cast<size_t>(windowSamples / 2);
    for (size_t sample = 0; sample < count; ++sample) {
        const size_t begin = sample > radius ? sample - radius : 0;
        const size_t end = min(count, sample + radius + 1);
        neighborhood.assign(values.begin() + begin, values.begin() + end);
        const size_t middle = neighborhood.size() / 2;
        nth_element(neighborhood.begin(), neighborhood.begin() + middle,
                    neighborhood.end());
        filtered[sample] = neighborhood[middle];
    }
    return filtered;
}

EventIDResult DecodeEventID(const Waveform& waveform,
                            const OscilloscopeProcessingConfig& config) {
    EventIDResult result;
    const size_t count =
        min(waveform.time.size(), waveform.amplitude.size());
    if (count < static_cast<size_t>(kFrameBits)) return result;

    const vector<double> filteredAmplitude = MedianFilter(
        waveform.amplitude, count, config.eventIDMedianFilterSamples);
    const double low = Quantile(filteredAmplitude, 0.05);
    const double high = Quantile(filteredAmplitude, 0.95);
    result.threshold = 0.5 * (low + high);
    if (!isfinite(result.threshold) || high - low < 0.2) return result;

    const double timeWindowMin =
        config.eventIDTimeWindowNs[0] * 1.0e-9;
    const double timeWindowMax =
        config.eventIDTimeWindowNs[1] * 1.0e-9;
    const auto windowBegin =
        lower_bound(waveform.time.begin(), waveform.time.begin() + count,
                    timeWindowMin);
    const auto windowEnd =
        upper_bound(waveform.time.begin(), waveform.time.begin() + count,
                    timeWindowMax);
    size_t firstSample =
        static_cast<size_t>(windowBegin - waveform.time.begin());
    const size_t lastSample =
        static_cast<size_t>(windowEnd - waveform.time.begin());
    firstSample = max<size_t>(1, firstSample);

    size_t firstFalling = count;
    for (size_t sample = firstSample; sample < lastSample; ++sample) {
        if (filteredAmplitude[sample - 1] > result.threshold &&
            filteredAmplitude[sample] <= result.threshold) {
            firstFalling = sample;
            break;
        }
    }
    if (firstFalling == count) return result;

    const double frameStart =
        InterpolatedCrossing(waveform.time, filteredAmplitude,
                             firstFalling, result.threshold);
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
        bits[bit] = filteredAmplitude[sample] < result.threshold ? 1 : 0;
    }

    if (bits.front() != 1 || bits.back() != 0) return result;
    uint64_t eventID = 0;
    for (int bit = 1; bit <= kEventIDBits; ++bit)
        eventID = (eventID << 1U) | static_cast<uint64_t>(bits[bit]);
    result.eventID = eventID;
    result.frameStartNs = frameStart * 1.0e9;
    result.valid = true;
    return result;
}

TimingResult MeasureNegativePulse(
    const Waveform& waveform, double eventIDFrameStartNs,
    const OscilloscopeProcessingConfig& config) {
    TimingResult result;
    const size_t count =
        min(waveform.time.size(), waveform.amplitude.size());
    if (count < 10 || !isfinite(eventIDFrameStartNs)) return result;

    const double windowBeginTime =
        (eventIDFrameStartNs + config.signalTimeWindowNs[0]) * 1.0e-9;
    const double windowEndTime =
        (eventIDFrameStartNs + config.signalTimeWindowNs[1]) * 1.0e-9;
    const auto windowBeginIterator =
        lower_bound(waveform.time.begin(), waveform.time.begin() + count,
                    windowBeginTime);
    const auto windowEndIterator =
        upper_bound(waveform.time.begin(), waveform.time.begin() + count,
                    windowEndTime);
    const size_t windowBegin = static_cast<size_t>(
        windowBeginIterator - waveform.time.begin());
    const size_t windowEnd = static_cast<size_t>(
        windowEndIterator - waveform.time.begin());
    if (windowBegin < 2 || windowEnd <= windowBegin + 1) return result;

    vector<double> baselineValues(
        waveform.amplitude.begin(), waveform.amplitude.begin() + windowBegin);
    result.baseline = Quantile(baselineValues, 0.5);
    if (!isfinite(result.baseline)) return result;

    const auto minimum = min_element(
        waveform.amplitude.begin() + windowBegin,
        waveform.amplitude.begin() + windowEnd);
    const size_t minimumSample = static_cast<size_t>(
        minimum - waveform.amplitude.begin());
    result.amplitude = result.baseline - *minimum;
    if (!isfinite(result.amplitude) ||
        result.amplitude < config.minPulseAmplitude) {
        result.status = TimingStatus::BelowAmplitude;
        return result;
    }

    result.threshold =
        result.baseline - config.cfdFraction * result.amplitude;
    size_t selectedSample = 0;
    for (size_t sample = minimumSample; sample > windowBegin; --sample) {
        if (waveform.amplitude[sample - 1] > result.threshold &&
            waveform.amplitude[sample] <= result.threshold) {
            selectedSample = sample;
            break;
        }
    }
    if (selectedSample == 0) {
        result.status = TimingStatus::NoCrossing;
        return result;
    }

    result.timeNs =
        InterpolatedCrossing(waveform, selectedSample, result.threshold) *
        1.0e9;
    result.status =
        isfinite(result.timeNs) ? TimingStatus::Valid
                                : TimingStatus::NoCrossing;
    return result;
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

void CountTimingResult(
    TimingStatus status, OscilloscopeChannelStatistics& statistics) {
    switch (status) {
        case TimingStatus::Valid:
            ++statistics.valid;
            break;
        case TimingStatus::BelowAmplitude:
            ++statistics.belowAmplitude;
            break;
        case TimingStatus::NoCrossing:
            ++statistics.noCrossing;
            break;
        case TimingStatus::NoWindow:
            ++statistics.noWindow;
            break;
    }
}

}  // namespace

OscilloscopeProcessingResult OscilloscopeDataProcessor::Process(
    const filesystem::path& inputDirectory,
    const OscilloscopeProcessingConfig& config,
    const DiscoveryCallback& onDiscovery,
    const ProgressCallback& onProgress) const {
    OscilloscopeProcessingResult result;
    if (!filesystem::is_directory(inputDirectory)) {
        result.error =
            "CSV directory does not exist: " + inputDirectory.string();
        return result;
    }

    const auto traceFiles = FindTraceFiles(inputDirectory);
    vector<int> traceIndices;
    for (const auto& [index, channels] : traceFiles) {
        bool complete = true;
        for (int channel = 1; channel <= 4; ++channel)
            complete = complete && channels.count(channel) > 0;
        if (complete) traceIndices.push_back(index);
    }
    if (config.maxWaveformFiles > 0 &&
        traceIndices.size() >
            static_cast<size_t>(config.maxWaveformFiles))
        traceIndices.resize(
            static_cast<size_t>(config.maxWaveformFiles));
    result.traceFiles = traceIndices.size();
    if (onDiscovery) onDiscovery(result.traceFiles);
    if (traceIndices.empty()) {
        result.error =
            "no complete C1-C4 trace sets in " + inputDirectory.string();
        return result;
    }

    size_t processedFiles = 0;
    for (int index : traceIndices) {
        array<TraceReader, 4> readers;
        string error;
        for (size_t channel = 0; channel < readers.size(); ++channel) {
            if (!readers[channel].Open(
                    traceFiles.at(index).at(static_cast<int>(channel + 1)),
                    error)) {
                result.error = error;
                return result;
            }
        }

        int segments = readers.front().SegmentCount();
        for (const TraceReader& reader : readers)
            segments = min(segments, reader.SegmentCount());

        array<Waveform, 4> waveforms;
        for (int segment = 0; segment < segments; ++segment) {
            ++result.processedSegments;
            for (size_t channel = 0; channel < readers.size(); ++channel) {
                if (!readers[channel].ReadSegment(waveforms[channel])) {
                    result.error =
                        "incomplete segment " + to_string(segment) +
                        " in trace " + to_string(index);
                    return result;
                }
            }

            const EventIDResult decoded = DecodeEventID(waveforms[3], config);
            if (decoded.valid) {
                ++result.decodedEventIDs;
                result.uniqueEventIDs.insert(decoded.eventID);
            } else {
                ++result.invalidEventIDs;
            }

            array<TimingResult, kChannelCount> timing;
            for (size_t channel = 0; channel < timing.size(); ++channel) {
                timing[channel] = MeasureNegativePulse(
                    waveforms[channel], decoded.frameStartNs, config);
                CountTimingResult(
                    timing[channel].status,
                    result.channelStatistics[channel]);
            }

            const bool complete =
                decoded.valid &&
                all_of(timing.begin(), timing.end(),
                       [](const TimingResult& value) {
                           return value.status == TimingStatus::Valid;
                       });
            if (!complete) continue;

            OscilloscopeEvent event;
            event.eventID = decoded.eventID;
            event.eventIDTime = decoded.frameStartNs;
            for (size_t channel = 0; channel < timing.size(); ++channel) {
                event.time[channel] =
                    timing[channel].timeNs - decoded.frameStartNs;
                event.amplitude[channel] = timing[channel].amplitude;
            }
            result.pairDifferences[0].push_back(
                event.time[0] - event.time[1]);
            result.pairDifferences[1].push_back(
                event.time[0] - event.time[2]);
            result.pairDifferences[2].push_back(
                event.time[1] - event.time[2]);
            result.events.push_back(event);
        }

        ++processedFiles;
        if (onProgress) onProgress(processedFiles, traceIndices.size());
    }

    result.success = true;
    return result;
}
