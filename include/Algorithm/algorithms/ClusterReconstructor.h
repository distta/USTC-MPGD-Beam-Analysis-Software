#ifndef CLUSTER_RECONSTRUCTOR_H
#define CLUSTER_RECONSTRUCTOR_H

#include "Config.h"
#include "DataModel.h"
#include "IAlgorithm.h"
#include <vector>

// 前向声明，避免循环依赖
class DetectorFrame;

/**
 * @brief 簇重建算法 - 负责重建cluster的pos值
 *
 * 功能：
 * - 根据cluster中的StripHit信息重建cluster的位置（pos）
 * - 支持多种重建方法（电荷加权、UTPC等）
 */
class ClusterReconstructor : public IAlgorithm {
   public:
    ClusterReconstructor() = default;
    virtual ~ClusterReconstructor() = default;

    // 实现IAlgorithm接口
    std::string GetName() const override { return "ClusterReconstructor"; }
    std::string GetVersion() const override { return "1.0.0"; }

    void LoadConfig(const json& config) override {
        m_config.loadFrom(config);
    }

    void Print() const override {
        std::cout << "[" << GetName() << " v" << GetVersion() << "]" << std::endl;
        m_config.print();
    }

    // ========== 核心接口 ==========

    // 统一接口: 处理DetectorFrame中的Cluster数据
    bool Process(DetectorFrame& frame) override;

   private:
    ReconstructionConfig m_config;

    // 不同的重建方法（以stripID为单位）
    void reconstructChargeWeighted(Cluster& cluster, const std::vector<StripHit>& stripHits);
    void reconstructUTPC(Cluster& cluster, const std::vector<StripHit>& stripHits, double t0);
};

#endif  // CLUSTER_RECONSTRUCTOR_H
