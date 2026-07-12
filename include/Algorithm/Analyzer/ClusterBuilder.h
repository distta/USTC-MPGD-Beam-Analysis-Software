#ifndef CLUSTER_BUILDER_H
#define CLUSTER_BUILDER_H

#include "Config.h"
#include "DataModel.h"
#include "IAlgorithm.h"
#include <map>
#include <vector>

// 前向声明，避免循环依赖
class DetectorFrame;

/**
 * @brief 聚类构建算法 - 负责StripHit到Cluster的聚类和匹配
 *
 * 功能：
 * - 将StripHit按照stripID聚类成Cluster
 * - 对不同type的Cluster进行匹配（如X-Y匹配）
 * - 计算cluster的基本物理量（size、range、charge等）
 */
class ClusterBuilder : public IAlgorithm {
   public:
    ClusterBuilder() = default;
    virtual ~ClusterBuilder() = default;

    // 实现IAlgorithm接口
    std::string GetName() const override { return "ClusterBuilder"; }
    std::string GetVersion() const override { return "1.0.0"; }

    void LoadConfig(const json& config) override {
        m_config.loadFrom(config);
    }

    void Print() const override {
        std::cout << "[" << GetName() << " v" << GetVersion() << "]" << std::endl;
        m_config.print();
    }

    // ========== 核心接口 ==========

    // 统一接口: 处理DetectorFrame中的StripHit数据
    bool Process(DetectorFrame& frame) override;

    /**
     * @brief 从 StripHit 构建Cluster（保留旧接口以保证兼容性）
     * @param stripHits 已排序的StripHit列表（按type和stripID有序）
     * @return 所有type的Cluster混合列表
     */
    std::vector<Cluster> BuildClusters(const std::vector<StripHit>& stripHits);

   private:
    ClusterConfig m_config;

    /**
     * @brief 处理单个cluster（计算size、range、charge等）
     * @param cluster 要处理的cluster
     * @param stripHits StripHit列表（用于通过索引访问）
     * @return 是否有效（通过minClusterSize等过滤）
     */
    bool processCluster(Cluster& cluster, const std::vector<StripHit>& stripHits);
};


#endif  // CLUSTER_BUILDER_H
