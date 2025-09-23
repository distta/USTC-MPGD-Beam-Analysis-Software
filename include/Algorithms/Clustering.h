#pragma once

#include "Config.h"
#include "DataModel.h"
#include <vector>

class Clustering {
   public:
    explicit Clustering(const json& config);

    std::vector<std::vector<RawData>> preClustering(const std::vector<RawData>& rawData) const;

    RecCluster BuildCluster(const std::vector<RawData>& raw);

   private:
    WaveformConfig m_waveformConfig;
    ClusterConfig m_clusterConfig;
    ReconstructionConfig m_reconstructionConfig;
    /**
     * @brief 处理波形数据
     */
    void processWaveform(const RawData& rawData, StripHit& stripData);

    /**
     * @brief 处理高质量波形（使用完整拟合算法）
     */
    void processHighQualityWaveform(const std::vector<short>& waveform, int peakAmp, int peakTime, int firstOverTh, StripHit& stripData);

    /**
     * @brief 处理低质量波形（使用简化CFD算法）
     */
    void processSimpleCFD(const std::vector<short>& waveform, int peakAmp, int peakTime, int firstOverTh, int inducedCharge, StripHit& stripData);

    /**
     * @brief 处理聚类
     */
    bool processCluster(RecCluster& cluster);

    /**
     * @brief 电荷加权重建
     */
    void reconstructChargeWeighted(RecCluster& cluster);

    /**
     * @brief UTPC重建
     */
    void reconstructUTPC(RecCluster& cluster);
};