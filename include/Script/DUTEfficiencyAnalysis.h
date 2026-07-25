#pragma once

#include "Detector/Detector.h"
#include "Event/DataModel.h"

#include <TDirectory.h>

#include <memory>
#include <unordered_set>
#include <vector>

namespace DUTEfficiency {

struct Config {
    int xBins = 8;
    int yBins = 7;
    double xMin = 0.0;
    double xMax = 100.0;
    double yMin = 0.0;
    double yMax = 100.0;
    int minEntriesPerBin = 1;
    std::vector<int> excludedXBins;
    std::vector<int> excludedYBins;

    double margin = 0.6;
    double envelopeScanMin = 0.0;
    double envelopeScanMax = 2.0;
    double envelopeScanStep = 0.05;

    bool enableFakeEfficiency = true;
    unsigned int fakeSeed = 12345;
    int fakePartnersPerEvent = 20;

};

struct Result {
    long long eligibleEvents = 0;
    double eventWeighted2D = 0.0;
    double fake2D = 0.0;
    double nonuniformityX = 0.0;
    double nonuniformityY = 0.0;
    double nonuniformity2D = 0.0;
};

Result Analyze(const std::vector<Event>& events,
               const std::unordered_set<int>& strictSingleHitTrackerEvents,
               const std::shared_ptr<Detector>& detector,
               const Config& config,
               TDirectory* outputDirectory);

}  // namespace DUTEfficiency
