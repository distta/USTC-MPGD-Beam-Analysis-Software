#pragma once

#include "Script/Base/IScript.h"

#include <string>

class TimeResolutionScript : public IScript {
   public:
    std::string GetName() const override { return "TimeResolutionScript"; }
    std::string GetDescription() const override { return "Tracker timing resolution with pair and oscilloscope-T0 references"; }
    void LoadConfig(const json& config) override;
    void Print() const override;
    bool Validate() const override;
    bool Execute() override;

   private:
    std::string m_csvDirectory;
    std::string m_trackFile;
    std::string m_outputFile = "TimeResolution.root";
    std::string m_timingCacheFile = "OscilloscopeTimingCache.root";
    std::string m_waveformDiagnosticFile = "WaveformDiagnostics.root";
    bool m_rebuildTimingCache = false;
    bool m_writeWaveformDiagnostics = false;
    double m_triggerThreshold = 1.0;
    double m_eventThreshold = 0.5;
    int m_eventBits = 16;
    int m_histogramBins = 160;
    int m_maxWaveformFiles = -1;
};
