#ifndef CLUSTER_BUILDER_H
#define CLUSTER_BUILDER_H

#include "Config.h"
#include "DataModel.h"
#include "IAlgorithm.h"
#include <map>
#include <vector>

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

    /**
     * @brief 从StripHit构建Cluster
     * @param stripHitsByType 按type分组的StripHit映射
     * @return 匹配好的RecCluster列表
     */
    std::vector<RecCluster> BuildClusters(const std::map<int, std::vector<StripHit>>& stripHitsByType);

   private:
    ClusterConfig m_config;

    /**
     * @brief 对单个type的StripHit进行聚类
     * @param stripHits 同一type的StripHit列表
     * @param type StripHit的type
     * @return Cluster列表
     */
    std::vector<Cluster> buildClustersForType(const std::vector<StripHit>& stripHits, int type);

    /**
     * @brief 处理单个cluster（计算size、range、charge等）
     * @param cluster 要处理的cluster
     * @return 是否有效（通过minClusterSize等过滤）
     */
    bool processCluster(Cluster& cluster);

    /**
     * @brief 匹配不同type的Cluster
     * @param clustersByType 按type分组的Cluster
     * @return 匹配好的RecCluster列表
     */
    std::vector<RecCluster> matchClusters(std::map<int, std::vector<Cluster>>& clustersByType);
};


#endif  // CLUSTER_BUILDER_H
