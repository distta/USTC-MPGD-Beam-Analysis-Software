#pragma once

#include "Script/Base/IScript.h"

#include <string>
#include <vector>

class TimeResolutionScript : public IScript {
   public:
    std::string GetName() const override { return "TimeResolutionScript"; }
    std::string GetDescription() const override {
        return "Three-tracker reference and DUT timing resolution";
    }
    void LoadConfig(const json& config) override;
    void Print() const override;
    bool Validate() const override;
    bool Execute() override;

   private:
    std::string m_trackFile;
    std::string m_dutFile = "PadDUTInfo.root";
    std::string m_outputFile = "TimeResolution.root";
    int m_histogramBins = 160;
    std::vector<int> m_trackerIDs;
    bool m_analyzeDUTTiming = false;
    bool m_applyDUTTimeWalkCorrection = false;
    json m_timingWaveformConfig = json{{"mode", "Fit"}};
};
