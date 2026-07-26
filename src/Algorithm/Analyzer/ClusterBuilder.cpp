#include "Algorithm/Analyzer/ClusterBuilder.h"
#include "AlgorithmFactory.h"
#include "Detector/Detector.h"
#include "Detector/PlanarPad.h"
#include "DetectorFrame.h"
#include <algorithm>
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <unordered_map>
#include <utility>

REGISTER_ALGORITHM("ClusterBuilder", ClusterBuilder)

bool ClusterBuilder::Process(DetectorFrame& frame) {
    const auto& channelHits = frame.ChannelHits();
    if (channelHits.empty()) return false;

    auto clusters = m_detector && m_detector->GetPlanarPadConfig()
                        ? BuildPadClusters(channelHits)
                        : BuildClusters(channelHits);

    // 将结果写入frame
    auto& frameClusters = frame.GetMutableClusters();
    frameClusters = std::move(clusters);

    return !frameClusters.empty();
}

std::vector<Cluster> ClusterBuilder::BuildClusters(const std::vector<ChannelHit>& channelHits) {
    std::vector<Cluster> clusters;
    if (channelHits.empty()) return clusters;

    std::map<int, std::vector<int>> indicesByType;
    for (size_t index = 0; index < channelHits.size(); ++index) {
        if (channelHits[index].isValid)
            indicesByType[channelHits[index].type].push_back(
                static_cast<int>(index));
    }

    for (const auto& [type, indices] : indicesByType) {
        std::vector<int> parent(indices.size());
        std::vector<int> rank(indices.size(), 0);
        for (size_t index = 0; index < parent.size(); ++index)
            parent[index] = static_cast<int>(index);

        const auto findRoot = [&](int value, const auto& self) -> int {
            if (parent[value] != value)
                parent[value] = self(parent[value], self);
            return parent[value];
        };
        const auto merge = [&](int left, int right) {
            left = findRoot(left, findRoot);
            right = findRoot(right, findRoot);
            if (left == right) return;
            if (rank[left] < rank[right]) std::swap(left, right);
            parent[right] = left;
            if (rank[left] == rank[right]) ++rank[left];
        };

        for (size_t left = 0; left < indices.size(); ++left) {
            const auto& leftHit = channelHits[indices[left]];
            for (size_t right = left + 1; right < indices.size(); ++right) {
                const auto& rightHit = channelHits[indices[right]];
                if (rightHit.id0 - leftHit.id0 > 1 + m_config.maxGap)
                    break;
                if (std::abs(rightHit.time - leftHit.time) <=
                    m_config.timeWindowNs) {
                    merge(static_cast<int>(left), static_cast<int>(right));
                }
            }
        }

        std::map<int, std::vector<int>> components;
        for (size_t index = 0; index < indices.size(); ++index) {
            components[findRoot(static_cast<int>(index), findRoot)]
                .push_back(indices[index]);
        }
        for (auto& [root, members] : components) {
            (void)root;
            Cluster cluster{};
            cluster.type = type;
            cluster.channelHitIndices = std::move(members);
            if (processCluster(cluster, channelHits))
                clusters.push_back(std::move(cluster));
        }
    }

    std::sort(clusters.begin(), clusters.end(),
              [](const Cluster& left, const Cluster& right) {
                  return left.channelHitIndices.front() <
                         right.channelHitIndices.front();
              });

    return clusters;
}

std::vector<Cluster> ClusterBuilder::BuildPadClusters(
    const std::vector<ChannelHit>& channelHits) {
    std::vector<Cluster> clusters;
    const auto* detector = dynamic_cast<const PlanarPad*>(m_detector);
    if (!detector) return clusters;
    const auto& config = detector->GetPadConfig();

    using PadCoordinate = std::pair<int, int>;
    std::map<PadCoordinate, std::vector<int>> hitsByPad;
    for (size_t index = 0; index < channelHits.size(); ++index) {
        const auto& hit = channelHits[index];
        if (!hit.isValid || hit.type != config.planeType) continue;
        if (!hit.HasID1()) {
            std::cerr << "[ClusterBuilder] WARNING: ignoring PlanarPad hit "
                      << index << " because id1(row) is missing\n";
            continue;
        }
        const int row = hit.id1;
        const int column = hit.id0;
        if (row < 0 || row >= config.rows ||
            column < 0 || column >= config.columns) {
            std::cerr << "[ClusterBuilder] WARNING: ignoring PlanarPad hit "
                      << index << " with out-of-range coordinate (column="
                      << column << ", row=" << row << ") for "
                      << config.columns << 'x' << config.rows << " pads\n";
            continue;
        }
        hitsByPad[{row, column}].push_back(static_cast<int>(index));
    }

    std::set<PadCoordinate> visited;
    for (const auto& [seedPad, seedHits] : hitsByPad) {
        (void)seedHits;
        if (!visited.insert(seedPad).second) continue;

        Cluster cluster{};
        cluster.type = config.planeType;
        std::queue<PadCoordinate> pending;
        pending.push(seedPad);
        int minRow = config.rows, maxRow = -1;
        int minColumn = config.columns, maxColumn = -1;
        int padCount = 0;

        while (!pending.empty()) {
            const PadCoordinate coordinate = pending.front();
            pending.pop();
            ++padCount;
            const auto [row, column] = coordinate;
            minRow = std::min(minRow, row);
            maxRow = std::max(maxRow, row);
            minColumn = std::min(minColumn, column);
            maxColumn = std::max(maxColumn, column);
            const auto& indices = hitsByPad.at(coordinate);
            cluster.channelHitIndices.insert(cluster.channelHitIndices.end(),
                                             indices.begin(), indices.end());

            for (int deltaRow = -1; deltaRow <= 1; ++deltaRow) {
                for (int deltaColumn = -1; deltaColumn <= 1; ++deltaColumn) {
                    if (deltaRow == 0 && deltaColumn == 0) continue;
                    if (m_config.connectivity == 4 &&
                        std::abs(deltaRow) + std::abs(deltaColumn) != 1) continue;
                    const int nextRow = row + deltaRow;
                    const int nextColumn = column + deltaColumn;
                    if (nextRow < 0 || nextRow >= config.rows || nextColumn < 0 ||
                        nextColumn >= config.columns) continue;
                    const PadCoordinate nextPad{nextRow, nextColumn};
                    if (hitsByPad.count(nextPad) && visited.insert(nextPad).second)
                        pending.push(nextPad);
                }
            }
        }

        std::sort(cluster.channelHitIndices.begin(), cluster.channelHitIndices.end());
        if (processCluster(cluster, channelHits, padCount)) {
            cluster.range = std::max(maxRow - minRow, maxColumn - minColumn) + 1;
            clusters.push_back(std::move(cluster));
        }
    }

    std::sort(clusters.begin(), clusters.end(), [](const Cluster& lhs, const Cluster& rhs) {
        return lhs.channelHitIndices.front() < rhs.channelHitIndices.front();
    });
    return clusters;
}

bool ClusterBuilder::processCluster(
    Cluster& cluster,
    const std::vector<ChannelHit>& channelHits,
    int logicalSize) {
    if (cluster.channelHitIndices.empty()) return false;

    // 计算size和range
    if (logicalSize >= 0) {
        cluster.size = logicalSize;
    } else {
        std::set<int> uniqueIDs;
        for (int index : cluster.channelHitIndices)
            uniqueIDs.insert(channelHits[index].id0);
        cluster.size = static_cast<int>(uniqueIDs.size());
    }

    int minID = std::numeric_limits<int>::max();
    int maxID = std::numeric_limits<int>::min();
    for (int index : cluster.channelHitIndices) {
        minID = std::min(minID, channelHits[index].id0);
        maxID = std::max(maxID, channelHits[index].id0);
    }
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

    for (int idx : cluster.channelHitIndices) {
        const auto& strip = channelHits[idx];

        cluster.charge += strip.charge;
        sumCharge += strip.charge;
        sumPos += strip.id0 * strip.charge;  // strip 电荷加权位置

        if (strip.amp > cluster.maxAmp) {
            cluster.maxAmp = strip.amp;
        }
        if (strip.time < cluster.time) {
            cluster.time = strip.time;
        }
    }

    // 计算质心
    cluster.centroid = (sumCharge > 0) ? sumPos / sumCharge
                                      : 0.5 * (minID + maxID);

    // 初始化pos为0（需要由ClusterReconstructor重建）
    cluster.pos = 0.0;

    return true;
}
