#include "algorithms/ClusterBuilder.h"
#include "AlgorithmFactory.h"
#include "Detector/Detector.h"
#include <algorithm>
#include <limits>

REGISTER_ALGORITHM("ClusterBuilder", ClusterBuilder)

std::vector<RecCluster> ClusterBuilder::BuildClusters(const std::map<int, std::vector<StripHit>>& stripHitsByType) {
    // Step 1: 对每个type的StripHit进行聚类
    std::map<int, std::vector<Cluster>> clustersByType;

    for (const auto& [type, stripHits] : stripHitsByType) {
        clustersByType[type] = buildClustersForType(stripHits, type);
    }

    // Step 2: 匹配不同type的Cluster
    return matchClusters(clustersByType);
}

std::vector<Cluster> ClusterBuilder::buildClustersForType(const std::vector<StripHit>& stripHits, int type) {
    std::vector<Cluster> clusters;
    if (stripHits.empty()) return clusters;

    // 按stripID排序
    std::vector<StripHit> sortedHits = stripHits;
    std::sort(sortedHits.begin(), sortedHits.end(),
              [](const StripHit& a, const StripHit& b) { return a.stripID < b.stripID; });

    // 聚类：连续的stripID（允许maxGap间隙）聚成一个cluster
    std::vector<StripHit> currentGroup;
    currentGroup.push_back(sortedHits.front());

    for (size_t i = 1; i < sortedHits.size(); ++i) {
        if (sortedHits[i].stripID <= sortedHits[i - 1].stripID + 1 + m_config.maxGap) {
            if (sortedHits[i].isValid)
                currentGroup.push_back(sortedHits[i]);
        } else {
            // 完成一个cluster
            Cluster cluster;
            cluster.type = type;
            cluster.strips = currentGroup;

            if (processCluster(cluster)) {
                clusters.push_back(std::move(cluster));
            }

            currentGroup.clear();
            currentGroup.push_back(sortedHits[i]);
        }
    }

    // 处理最后一个cluster
    if (!currentGroup.empty()) {
        Cluster cluster;
        cluster.type = type;
        cluster.strips = currentGroup;

        if (processCluster(cluster)) {
            clusters.push_back(std::move(cluster));
        }
    }

    return clusters;
}

bool ClusterBuilder::processCluster(Cluster& cluster) {
    if (cluster.strips.empty()) return false;

    // 计算size和range
    cluster.size = cluster.strips.size();

    int minID = cluster.strips.front().stripID;
    int maxID = cluster.strips.back().stripID;
    cluster.range = maxID - minID + 1;

    // 过滤：检查cluster大小
    if (cluster.size < m_config.minClusterSize || cluster.size > m_config.maxClusterSize) {
        return false;
    }

    // 计算总电荷和最大幅度
    cluster.charge = 0.0;
    cluster.maxAmp = 0.0;
    cluster.time = std::numeric_limits<double>::max();

    for (const auto& strip : cluster.strips) {
        cluster.charge += strip.charge;
        if (strip.amp > cluster.maxAmp) {
            cluster.maxAmp = strip.amp;
        }
        if (strip.time < cluster.time) {
            cluster.time = strip.time;
        }
    }

    // 初始化matchID为-1（未匹配）
    cluster.matchID = -1;

    // 初始化pos为0（需要由ClusterReconstructor重建）
    cluster.pos = 0.0;

    return true;
}

std::vector<RecCluster> ClusterBuilder::matchClusters(std::map<int, std::vector<Cluster>>& clustersByType) {
    std::vector<RecCluster> recClusters;

    // 如果只有一个type，直接转换
    int matchCount = 0;
    if (m_detector->getConfig().readoutPlanePitch.size() == 1) {
        for (auto& cluster : clustersByType.begin()->second) {
            RecCluster rec;
            cluster.matchID = matchCount;
            rec.push_back(cluster);
            recClusters.push_back(rec);
            matchCount++;
        }
        return recClusters;
    }

    // 多type匹配（如X-Y匹配）
    // 简化版本：基于电荷差进行匹配
    std::vector<int> types;
    for (const auto& [type, _] : clustersByType) {
        types.push_back(type);
    }

    if (types.size() == 2) {

        if (m_detector->isTracker()) {
            for (auto& cluster : clustersByType[types[0]]) {
                for (auto& cluster1 : clustersByType[types[1]]) {
                    RecCluster rec;
                    cluster.matchID = matchCount;
                    cluster1.matchID = matchCount;
                    rec.push_back(cluster);
                    rec.push_back(cluster1);
                    recClusters.push_back(rec);
                    matchCount++;
                }
            }
            return recClusters;
        }
        // 双type匹配
        auto& clusters0 = clustersByType[types[0]];
        auto& clusters1 = clustersByType[types[1]];

        std::vector<bool> used0(clusters0.size(), false);
        std::vector<bool> used1(clusters1.size(), false);

        int matchCount = 0;
        for (size_t i = 0; i < clusters0.size(); ++i) {
            if (used0[i]) continue;

            double minChargeDiff = std::numeric_limits<double>::max();
            int bestMatch = -1;

            for (size_t j = 0; j < clusters1.size(); ++j) {
                if (used1[j]) continue;

                double chargeDiff = abs(clusters0[i].charge / clusters1[j].charge - 1.2);
                if (chargeDiff < minChargeDiff && chargeDiff < m_config.MaxChargeDiff) {
                    minChargeDiff = chargeDiff;
                    bestMatch = j;
                }
            }

            if (bestMatch != -1) {
                // 找到匹配
                RecCluster rec;
                clusters0[i].matchID = matchCount;
                clusters1[bestMatch].matchID = matchCount;
                rec.push_back(clusters0[i]);
                rec.push_back(clusters1[bestMatch]);
                recClusters.push_back(rec);
                matchCount++;

                used0[i] = true;
                used1[bestMatch] = true;
            }
        }
    }

    return recClusters;
}
