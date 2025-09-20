#pragma once

#include "Config.h"
#include "DataModel.h"
#include <vector>

class ClusterMatcher {
   public:
    explicit ClusterMatcher(const json& config = json{});

    std::vector<LocalHit> MatchClusters(const std::vector<RecCluster>& clusters);
};