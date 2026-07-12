#pragma once

#include "Algorithm/IAlgorithm.h"

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <map>
#include <set>
#include <string>
#include <vector>

class TDirectory;

struct OscilloscopeData {
    uint64_t decodedEventID = 0;
    double triggerTime = 0.0;
    double EvnentIDtime = 0.0;
    double truthT0 = 0.0;
    int fileIndex = -1;
    int segmentIndex = -1;
};

struct OscilloscopeDataResult {
    std::map<uint64_t, OscilloscopeData> dataByEventID;
    size_t decodedEntries = 0;
};

class OscilloscopeDataProcessor : public IAlgorithm {
   public:
    std::string GetName() const override { return "OscilloscopeDataProcessor"; }
    std::string GetVersion() const override { return "1.0.0"; }

    void LoadConfig(const json& config) override;
    void Print() const override;

    bool Initialize();
    OscilloscopeDataResult LoadData(const std::set<uint64_t>& wantedEventIDs) const;

   private:
    struct WaveformSegment {
        std::vector<double> time;
        std::vector<double> amplitude;
    };

    struct StreamTraceReader; 

    double DecodeArriveTime(const WaveformSegment& waveform, double threshold,
                            double startTime, bool risingFlag) const;
    uint64_t DecodeEventID(const WaveformSegment& eventCode, double startTime) const;
    bool DecodeOscilloscopeData(const WaveformSegment& eventCode,
                                const WaveformSegment& trigger,
                                OscilloscopeData& oscilloscopeData) const;

    void WriteWaveformDiagnostic(const WaveformSegment& eventCode,
                                 const WaveformSegment& trigger,
                                 const OscilloscopeData* oscilloscopeData,
                                 int fileIndex, int segmentIndex,
                                 TDirectory& directory) const;
                                 
    bool LoadOscilloscopeDataCache(const std::string& cachePath);
    void WriteOscilloscopeDataCache(const std::string& cachePath, const std::vector<OscilloscopeData>& decodedData) const;

    std::string m_csvDirectory;
    std::string m_outputDir;
    std::string m_dataCacheFile = "OscilloscopeDataCache.root";
    std::string m_waveformDiagnosticFile = "WaveformDiagnostics.root";
    bool m_rebuildTimingCache = false;
    bool m_writeWaveformDiagnostics = false;
    double m_triggerThreshold = 1.0;
    double m_eventThreshold = 0.5;
    int m_eventBits = 16;
    int m_eventCodeChannel = 3;
    int m_triggerChannel = 4;
    int m_maxWaveformFiles = -1;
    bool m_initialized = false;
    std::map<uint64_t, OscilloscopeData> m_dataByEventID;
    size_t m_decodedEntries = 0;
};
