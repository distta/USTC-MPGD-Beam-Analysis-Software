#pragma once

#include "Algorithm/TrackReconstruction.h"
#include "Algorithm/TrackerAlignment.h"
#include "Algorithm/TrackPerformance.h"
#include "Script/Base/IScript.h"

class TrackAnalysisScript : public IScript {
   public:
    std::string GetName() const override { return "TrackAnalysisScript"; }
    std::string GetDescription() const override { return "Multi-track reconstruction and robust tracker alignment"; }
    void LoadConfig(const json& config) override;
    void Print() const override;
    bool Execute() override;

   private:
    bool m_saveValidationData = true;
    bool m_debug = true;
    bool m_performanceHistograms = true;
    double m_residualHistogramRange = 2.0;
    int m_progressInterval = 10000;
    Tracking::Config m_tracking;
    Tracking::AlignmentConfig m_alignment;
};
