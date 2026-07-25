#pragma once

#include "Algorithm/Track/TrackReconstruction.h"
#include "Algorithm/Track/TrackerAlignment.h"
#include "Algorithm/Track/TrackPerformance.h"
#include "Script/Base/IScript.h"

class TrackAnalysisScript : public IScript {
   public:
    std::string GetName() const override { return "TrackAnalysisScript"; }
    std::string GetDescription() const override { return "Single-hit calibration followed by multi-hit track reconstruction"; }
    void LoadConfig(const json& config) override;
    void Print() const override;
    bool Execute() override;

   private:
    bool m_runAlignment = false;
    bool m_saveValidationData = true;
    bool m_debug = false;
    bool m_performanceHistograms = true;
    bool m_useEstimatedResolution = true;
    double m_residualHistogramRange = 2.0;
    Tracking::Config m_tracking;
    Tracking::AlignmentConfig m_alignment;
};
