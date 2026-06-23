#pragma once

#include "Algorithm/TrackReconstruction.h"
#include "Algorithm/TrackerAlignment.h"
#include <map>
#include <memory>
#include <vector>

class Detector;
class Event;
class TDirectory;
class TH1D;
class TH2D;

namespace Tracking {

class PerformanceAnalyzer {
   public:
    PerformanceAnalyzer(TDirectory* output, std::vector<std::shared_ptr<Detector>> detectors,
                        const Config& config, int totalEvents, double residualRange = 2.0);
    void RecordEvent(const Event& event, const std::vector<Result>& tracks,
                     const ReconstructionStats& stats);
    void RecordAlignment(const std::vector<AlignmentIteration>& history);
    void Write();

   private:
    struct DetectorHistograms {
        TH1D *hitMultiplicity{}, *clusterMultiplicity{}, *clusterSize{};
        TH1D *allHitX{}, *allHitY{}, *selectedHitX{}, *selectedHitY{};
        TH2D *allHitXY{}, *selectedHitXY{};
        TH1D *residualX{}, *residualY{}, *pullX{}, *pullY{}, *summary{};
        TH2D *residualXY{}, *residualXVsHitX{}, *residualYVsHitY{};
        TH2D *residualXVsSlope{}, *residualYVsSlope{}, *residualXVsEvent{}, *residualYVsEvent{};
        std::vector<double> residualValuesX, residualValuesY;
    };

    TDirectory* m_output{};
    std::vector<std::shared_ptr<Detector>> m_detectors;
    Config m_config;
    int m_totalEvents{};
    double m_residualRange{};
    std::map<int, DetectorHistograms> m_detectorHistograms;
    TH1D *m_tracksPerEvent{}, *m_seedCandidates{}, *m_finalCandidates{}, *m_conflictNodes{};
    TH1D *m_chi2{}, *m_layersPerTrack{}, *m_kx{}, *m_ky{}, *m_bx{}, *m_by{}, *m_angle{};
};

}  // namespace Tracking
