#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <set>
#include <string>
#include <vector>

struct OscilloscopeProcessingConfig {
    double cfdFraction = 0.5;
    double minPulseAmplitude = 0.1;
    int maxWaveformFiles = -1;
    std::array<double, 2> eventIDTimeWindowNs = {-100.0, 0.0};
    int eventIDMedianFilterSamples = 5;
    std::array<double, 2> signalTimeWindowNs = {-160.0, -80.0};
};

struct OscilloscopeEvent {
    uint64_t eventID = 0;
    double eventIDTime = 0.0;
    std::array<double, 3> time{};
    std::array<double, 3> amplitude{};
};

struct OscilloscopeChannelStatistics {
    size_t valid = 0;
    size_t belowAmplitude = 0;
    size_t noCrossing = 0;
    size_t noWindow = 0;
};

struct OscilloscopeProcessingResult {
    bool success = false;
    std::string error;
    size_t traceFiles = 0;
    size_t processedSegments = 0;
    size_t decodedEventIDs = 0;
    size_t invalidEventIDs = 0;
    std::set<uint64_t> uniqueEventIDs;
    std::array<OscilloscopeChannelStatistics, 3> channelStatistics{};
    std::vector<OscilloscopeEvent> events;
    std::array<std::vector<double>, 3> pairDifferences;
};

class OscilloscopeDataProcessor {
   public:
    static constexpr size_t kChannelCount = 3;

    using DiscoveryCallback = std::function<void(size_t)>;
    using ProgressCallback = std::function<void(size_t, size_t)>;

    OscilloscopeProcessingResult Process(
        const std::filesystem::path& inputDirectory,
        const OscilloscopeProcessingConfig& config,
        const DiscoveryCallback& onDiscovery = {},
        const ProgressCallback& onProgress = {}) const;
};
