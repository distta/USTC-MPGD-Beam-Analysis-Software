#pragma once
#include "Algorithm/Analyzer/ClusterBuilder.h"
#include "Algorithm/Analyzer/ClusterReconstructor.h"
#include "Algorithm/Analyzer/HitProcessor.h"
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

    double GetT0() const { return t0; }
    const std::vector<RawData>& Raw() const { return m_raw; }
    const std::vector<ChannelHit>& ChannelHits() const { return m_channelHits; }
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
    std::vector<ChannelHit>& GetMutableChannelHits() { return m_channelHits; }
    std::vector<Cluster>& GetMutableClusters() { return m_clusters; }

    // 根据索引获取单个ChannelHit（带边界检查）
    const ChannelHit& GetChannelHit(int index) const {
        if (index < 0 || index >= static_cast<int>(m_channelHits.size())) {
            throw std::out_of_range("ChannelHit index out of range");
        }
        return m_channelHits[index];
    }

    std::vector<const ChannelHit*> GetChannelHitsFromCluster(const Cluster& cluster) const {
        std::vector<const ChannelHit*> result;
        for (int idx : cluster.channelHitIndices) {
            if (idx >= 0 && idx < static_cast<int>(m_channelHits.size()))
                result.push_back(&m_channelHits[idx]);
        }
        return result;
    }

    const RawData* GetRawFromChannel(const ChannelHit& sh) const {
        if (sh.rawIndices < 0 || sh.rawIndices >= static_cast<int>(m_raw.size())) return nullptr;
        return &m_raw[sh.rawIndices];
    }

    bool Process(double t0);
    bool Process();

    // Extract ChannelHit Info from raw Data
    bool AnalyzeRaw();
    // Clustering ChannelHit
    bool Clustering();
    // Reconstruction: 将 Cluster 匹配并转换为 LocalHit
    bool Reconstruct();

    bool TransformToGlobal();

    void clear() {
        m_raw.clear();
        m_channelHits.clear();
        m_clusters.clear();
        m_localHits.clear();
        m_globalHits.clear();
    }

    // Acess
    const Detector& det() const { return m_det; }

   private:
    double t0;
    //   Raw Level (from DAQ)
    std::vector<RawData> m_raw;

    //   Step 1. Channel Level
    std::vector<ChannelHit> m_channelHits;

    //   Step 2. Cluster Level
    std::vector<Cluster> m_clusters;

    //   Step 3. Local Reconstruction
    std::vector<LocalHit> m_localHits;

    //   Step 4. Global Coordinates
    std::vector<GlobalHit> m_globalHits;

   private:
    const Detector& m_det;
};
