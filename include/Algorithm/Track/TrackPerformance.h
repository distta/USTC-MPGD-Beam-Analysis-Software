#pragma once

#include "Algorithm/Track/TrackReconstruction.h"
#include <map>
#include <memory>
#include <utility>
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
                        std::vector<std::shared_ptr<Detector>> referenceDetectors,
                        const Config& config, double residualRange = 2.0);
    void RecordEvent(const Event& event);
    void RecordTrack(const Event& event, const Result& result);
    void Reset();
    std::pair<double, double> EstimateHitResolution();
    void Write();

   private:
    struct DetectorHistograms {
        TH1D *clusterMultiplicity{};
        TH2D *hitMap{};
        TH1D *residualX{}, *residualY{};
    };

    TDirectory* m_output{};
    std::vector<std::shared_ptr<Detector>> m_detectors;
    std::vector<std::shared_ptr<Detector>> m_referenceDetectors;
    Config m_config;
    std::map<int, DetectorHistograms> m_detectorHistograms;
    TH1D *m_commonEquivalentHitX{}, *m_commonEquivalentHitY{};

    void RecordSelection(const Event& event,
                         const std::map<int, int>& hitIndices);
};

}  // namespace Tracking
