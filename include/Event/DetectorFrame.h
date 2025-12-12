#pragma once
#include "Algorithm/algorithms/ClusterBuilder.h"
#include "Algorithm/algorithms/ClusterReconstructor.h"
#include "Algorithm/algorithms/WaveformProcessor.h"
#include "DataModel.h"
#include "Detector.h"
#include <stdexcept>

class DetectorFrame {
   public:
    explicit DetectorFrame(const Detector& det)
        : m_det(det) {}

    // Data Interfaces
    void SetRawData(const std::vector<RawData>& raw) { m_raw = raw; }
    void AddRawData(const RawData& raw) { m_raw.push_back(raw); }

    const std::vector<RawData>& Raw() const { return m_raw; }
    const std::vector<StripHit>& StripHits() const { return m_stripHits; }
    const std::vector<Cluster>& Clusters() const { return m_clusters; }
    const std::vector<Cluster> Clusters(int type) const {
        std::vector<Cluster> result;
        result.reserve(m_clusters.size());
        for (const auto& clu : m_clusters)
            if (clu.type == type)
                result.push_back(clu);
        return result;
    };
    const std::vector<LocalHit>& LocalHits() const { return m_localHits; }
    const std::vector<GlobalHit>& GlobalHits() const { return m_globalHits; }

    // 可修改数据访问接口（供算法使用）
    std::vector<StripHit>& GetMutableStripHits() { return m_stripHits; }
    std::vector<Cluster>& GetMutableClusters() { return m_clusters; }

    // 根据索引获取单个StripHit（带边界检查）
    const StripHit& GetStripHit(int index) const {
        if (index < 0 || index >= static_cast<int>(m_stripHits.size())) {
            throw std::out_of_range("StripHit index out of range");
        }
        return m_stripHits[index];
    }

    std::vector<const StripHit*> GetStripHitsFromCluster(const Cluster& cluster) const {
        std::vector<const StripHit*> result;
        for (int idx : cluster.stripHitIndices) {
            if (idx >= 0 && idx < static_cast<int>(m_stripHits.size()))
                result.push_back(&m_stripHits[idx]);
        }
        return result;
    }

    const RawData* GetRawFromStrip(const StripHit& sh) const {
        if (sh.rawIndices < 0 || sh.rawIndices >= static_cast<int>(m_raw.size())) return nullptr;
        return &m_raw[sh.rawIndices];
    }

    bool Process();

    // Extract StripHit Info from raw Data
    bool AnalyzeRaw();
    // Clustering StripHit
    bool Clustering();
    // Reconstruction: 将 Cluster 匹配并转换为 LocalHit
    bool Reconstruct();

    bool TransformToGlobal();

    void clear() {
        m_raw.clear();
        m_stripHits.clear();
        m_clusters.clear();
        m_localHits.clear();
        m_globalHits.clear();
    }

    // Acess
    const Detector& det() const { return m_det; }

   private:
    //   Raw Level (from DAQ)
    std::vector<RawData> m_raw;

    //   Step 1. Strip Level
    std::vector<StripHit> m_stripHits;

    //   Step 2. Cluster Level
    std::vector<Cluster> m_clusters;

    //   Step 3. Local Reconstruction
    std::vector<LocalHit> m_localHits;

    //   Step 4. Global Coordinates
    std::vector<GlobalHit> m_globalHits;

   private:
    const Detector& m_det;
};
