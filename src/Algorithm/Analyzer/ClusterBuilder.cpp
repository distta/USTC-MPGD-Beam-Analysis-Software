#include "algorithms/ClusterBuilder.h"
#include "AlgorithmFactory.h"
#include "Detector/Detector.h"
#include "DetectorFrame.h"
#include <algorithm>
#include <limits>

REGISTER_ALGORITHM("ClusterBuilder", ClusterBuilder)

bool ClusterBuilder::Process(DetectorFrame& frame) {
    const auto& stripHits = frame.StripHits();
    if (stripHits.empty()) return false;

    // 调用内部的BuildClusters方法
    auto clusters = BuildClusters(stripHits);

    // 将结果写入frame
    auto& frameClusters = frame.GetMutableClusters();
    frameClusters = std::move(clusters);

    return !frameClusters.empty();
}

std::vector<Cluster> ClusterBuilder::BuildClusters(const std::vector<StripHit>& stripHits) {
    std::vector<Cluster> clusters;
    if (stripHits.empty()) return clusters;

    std::vector<int> currentGroupIndices;
    int currentType = -1;
    bool initialized = false;

    for (size_t i = 0; i < stripHits.size(); ++i) {

        const auto& currentHit = stripHits[i];

        if (!currentHit.isValid) continue;

        if (!initialized) {
            currentType = currentHit.type;
            initialized = true;
        }

        bool shouldEndCluster = false;

        if (!currentGroupIndices.empty()) {

            const auto& prevHit = stripHits[currentGroupIndices.back()];

            if (currentHit.type != currentType) {
                shouldEndCluster = true;
            } else if (currentHit.ID > prevHit.ID + 1 + m_config.maxGap) {
                shouldEndCluster = true;
            }
        }

        if (shouldEndCluster) {

            Cluster cluster;
            cluster.type = currentType;
            cluster.stripHitIndices = currentGroupIndices;

            if (processCluster(cluster, stripHits))
                clusters.push_back(std::move(cluster));

            currentGroupIndices.clear();
            currentType = currentHit.type;
        }

        currentGroupIndices.push_back(i);
    }

    if (!currentGroupIndices.empty()) {

        Cluster cluster;
        cluster.type = currentType;
        cluster.stripHitIndices = currentGroupIndices;

        if (processCluster(cluster, stripHits))
            clusters.push_back(std::move(cluster));
    }

    return clusters;
}

bool ClusterBuilder::processCluster(Cluster& cluster, const std::vector<StripHit>& stripHits) {
    if (cluster.stripHitIndices.empty()) return false;

    // 计算size和range
    cluster.size = cluster.stripHitIndices.size();

    int minID = stripHits[cluster.stripHitIndices.front()].ID;
    int maxID = stripHits[cluster.stripHitIndices.back()].ID;
    cluster.range = maxID - minID + 1;

    // 过滤：检查cluster大小
    if (cluster.size < m_config.minClusterSize || cluster.size > m_config.maxClusterSize) {
        return false;
    }

    // 计算总电荷和最大幅度
    cluster.charge = 0.0;
    cluster.maxAmp = 0.0;
    cluster.time = std::numeric_limits<double>::max();
    double sumPos = 0.0;
    double sumCharge = 0.0;

    for (int idx : cluster.stripHitIndices) {
        const auto& strip = stripHits[idx];

        cluster.charge += strip.amp;
        sumCharge += strip.charge;
        sumPos += strip.ID * strip.charge;  // 电荷加权位置

        if (strip.amp > cluster.maxAmp) {
            cluster.maxAmp = strip.amp;
        }
        if (strip.time < cluster.time) {
            cluster.time = strip.time;
        }
    }

    // 计算质心
    cluster.centroid = (sumCharge > 0) ? sumPos / sumCharge : 0.0;

    // 初始化pos为0（需要由ClusterReconstructor重建）
    cluster.pos = 0.0;

    return true;
}
