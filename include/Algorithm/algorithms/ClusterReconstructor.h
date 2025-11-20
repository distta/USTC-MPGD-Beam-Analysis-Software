#ifndef CLUSTER_RECONSTRUCTOR_H
#define CLUSTER_RECONSTRUCTOR_H

#include "IAlgorithm.h"
#include "Config.h"
#include "DataModel.h"
#include <vector>

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
    
    /**
     * @brief 重建RecCluster中所有Cluster的pos值（以stripID为单位）
     * @param recClusters RecCluster列表（输入/输出，会修改cluster的pos）
     */
    void ReconstructPositions(std::vector<RecCluster>& recClusters);
    
    /**
     * @brief 重建单个Cluster的pos值（以stripID为单位）
     * @param cluster Cluster（输入/输出，会修改pos）
     */
    void ReconstructPosition(Cluster& cluster);

private:
    ReconstructionConfig m_config;
    
    // 不同的重建方法（以stripID为单位）
    void reconstructChargeWeighted(Cluster& cluster);
    void reconstructUTPC(Cluster& cluster);
};
    

#endif // CLUSTER_RECONSTRUCTOR_H
