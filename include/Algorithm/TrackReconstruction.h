#pragma once

#include "Event/DataModel.h"
#include <map>
#include <memory>
#include <vector>

class Detector;

namespace Tracking {

struct Config {
    double resolutionX = 0.12;       // mm
    double resolutionY = 0.12;       // mm
    double gateSigma = 5.0;
    double maxChi2Ndf = 25.0;
    int maxBranchesPerLayer = 3;
    int maxCandidates = 4000;
    int maxTracks = 32;
    int conflictSearchNodes = 200000;
};

struct Result {
    Track track{};
    std::map<int, int> hitIndices;
    double chi2Ndf = 0.0;
};

struct ReconstructionStats {
    int seedCandidates = 0;
    int finalCandidates = 0;
    int conflictSearchNodes = 0;
};

class Reconstructor {
   public:
    Reconstructor(std::vector<std::shared_ptr<Detector>> detectors, Config config);
    std::vector<Result> Reconstruct(const Event& event, ReconstructionStats* stats = nullptr) const;
    const Config& GetConfig() const { return m_config; }

   private:
    std::vector<std::shared_ptr<Detector>> m_detectors;
    Config m_config;
};

Track FitWeighted(const std::vector<TVector3>& hits, double sigmaX, double sigmaY);

}  // namespace Tracking
