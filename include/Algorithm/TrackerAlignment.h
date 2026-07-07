#pragma once

#include "Algorithm/TrackReconstruction.h"
#include <memory>
#include <map>
#include <vector>

class Detector;

namespace Tracking {

struct AlignmentConfig {
    int maxIterations = 20;
    int minTracks = 20;
    int maxEvents = 10000;          // uniformly sampled events used to build alignment tracks
    int maxTracks = 5000;           // bound the cost of every objective evaluation
    int maxFunctionCalls = 600;
    double huberK = 2.5;
    double maxShiftStep = 0.20;     // mm / iteration
    double maxRotationStep = 0.001; // rad / iteration
    double shiftTolerance = 0.001;
    double rotationTolerance = 1e-6;
    bool debug = true;
};

struct AlignmentIteration {
    int iteration = 0;
    int tracks = 0;
    double lossBefore = 0;
    double lossAfter = 0;
    double maxShift = 0;
    double maxRotation = 0;
    std::map<int, TVector3> alignPosition;
    std::map<int, TVector3> alignRotation;
};

class Aligner {
   public:
    Aligner(std::vector<std::shared_ptr<Detector>> detectors,
            const Reconstructor& reconstructor, AlignmentConfig config);
    bool Run(const std::vector<Event>& events);
    const std::vector<AlignmentIteration>& History() const { return m_history; }

   private:
    std::vector<std::shared_ptr<Detector>> m_detectors;
    const Reconstructor& m_reconstructor;
    AlignmentConfig m_config;
    std::vector<AlignmentIteration> m_history;
};

}  // namespace Tracking
