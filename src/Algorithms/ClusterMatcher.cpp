#include "Algorithms/ClusterMatcher.h"
#include "DataModel.h"
#include <cmath>

ClusterMatcher::ClusterMatcher(const json& config) {
}

std::vector<LocalHit> ClusterMatcher::MatchClusters(const std::vector<RecCluster>& clusters) {
    std::vector<LocalHit> localHits;

    if (clusters.empty()) return localHits;

    std::unordered_map<int, std::vector<const RecCluster*>> clustersByType;
    for (const auto& aCluster : clusters) {
        clustersByType[aCluster.type].push_back(&aCluster);
    }

    for (const auto& [type, clusterVec] : clustersByType) {
        if (clusterVec.size() != 1)
            return {};
    }

    auto it0 = clustersByType.find(0);
    if (it0 == clustersByType.end() || it0->second.empty()) {
        return {};
    }

    LocalHit aLocalHit;
    aLocalHit.u = it0->second[0]->pos;
    
    auto it1 = clustersByType.find(1);
    if (it1 != clustersByType.end() && !it1->second.empty()) {
        aLocalHit.v = it1->second[0]->pos;
    }

    localHits.push_back(aLocalHit);
    return localHits;
}