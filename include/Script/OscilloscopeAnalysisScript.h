#pragma once

#include "Script/Base/IScript.h"

#include <array>
#include <string>

class OscilloscopeAnalysisScript : public IScript {
   public:
    std::string GetName() const override {
        return "OscilloscopeAnalysisScript";
    }
    std::string GetDescription() const override {
        return "Decode oscilloscope event IDs and measure three scintillator times";
    }

    void LoadConfig(const json& config) override;
    void Print() const override;
    bool Validate() const override;
    bool Execute() override;

   private:
    std::string m_csvDirectory;
    std::string m_outputFile = "OscilloscopeAnalysis.root";
    double m_cfdFraction = 0.5;
    double m_minPulseAmplitude = 0.1;
    int m_histogramBins = 160;
    int m_maxWaveformFiles = -1;
    std::array<double, 2> m_eventIDTimeWindowNs = {-100.0, 0.0};
    int m_eventIDMedianFilterSamples = 5;
    std::array<double, 2> m_signalTimeWindowNs = {-160.0, -80.0};
};
