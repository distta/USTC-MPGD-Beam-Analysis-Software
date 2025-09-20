#pragma once

#include "Config.h"
#include "DataModel.h"
#include <vector>

class Clustering {
   public:
    explicit Clustering(const json& config);

    std::vector<RecCluster> BuildCluster(const std::vector<RawData>& raw);

   private:
    WaveformConfig m_waveformConfig;
    ClusterConfig m_clusterConfig;
    ReconstructionConfig m_reconstructionConfig;

    /**
     * @brief 处理波形数据
     */
    void processWaveform(const RawData& rawData, StripHit& stripData);

    /**
     * @brief 处理聚类
     */
    void processCluster(RecCluster& cluster, std::vector<RecCluster>& clusters);

    /**
     * @brief 电荷加权重建
     */
    void reconstructChargeWeighted(RecCluster& cluster);

    /**
     * @brief UTPC重建
     */
    void reconstructUTPC(RecCluster& cluster);
};