#pragma once

#include "Config.h"
#include "DataModel.h"
#include <vector>

class Clustering {
   public:
    explicit Clustering(const json& config);

    std::vector<std::vector<RawData>> preClustering(const std::vector<RawData>& rawData) const;

    std::vector<RecHit> Reconstruction(const std::vector<RawData>& raw);

    std::vector<RecCluster> BuildClusters(const std::vector<RawData>& raws);

   private:
    WaveformConfig m_waveformConfig;
    ClusterConfig m_clusterConfig;
    ReconstructionConfig m_reconstructionConfig;
    /**
     * @brief 处理波形数据
     */
    bool processWaveform(const RawData& rawData, StripHit& stripData);

    void MatchClusters(std::map<int, std::vector<RecCluster>>&, std::vector<RecHit>& recHits);

    /**
     * @brief 处理聚类
     */
    bool
        processCluster(RecCluster& cluster);

    /**
     * @brief 电荷加权重建
     */
    void reconstructChargeWeighted(RecCluster& cluster);

    /**
     * @brief UTPC重建
     */
    void reconstructUTPC(RecCluster& cluster);
};