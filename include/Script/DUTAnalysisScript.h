#pragma once

#ifndef DUT_ANALYSIS_SCRIPT_H
#define DUT_ANALYSIS_SCRIPT_H

#include "Script/Base/IScript.h"
#include "Event/DataModel.h"
#include <string>
#include <map>
#include <vector>

// DUT分析配置类（从AnalysisEngine迁移）
class DUTAnalysisConfig {
   public:
    // 探测器类型常量
    static constexpr int kTypeX = 0;
    static constexpr int kTypeY = 1;

    // 5sigma残差筛选因子
    static constexpr double kSigmaFactor = 5.0;

    // 无效值标记
    static constexpr double kInvalidValue = -999.0;
    static constexpr int kInvalidSize = -1;

    // 分区配置
    struct BinningConfig {
        double predX_min;
        double predX_max;
        int nBinsX;
        double predY_min;
        double predY_max;
        int nBinsY;

        // 构造函数设置默认值
        BinningConfig()
            : predX_min(0), predX_max(100), nBinsX(5), predY_min(0), predY_max(100), nBinsY(5) {}
    };

    BinningConfig binning;

    // 判断位置是否在有效范围内
    bool IsInValidRange(double predX, double predY) const {
        return predX >= binning.predX_min && predX <= binning.predX_max &&
               predY >= binning.predY_min && predY <= binning.predY_max;
    }
};

// 分区统计数据（从AnalysisEngine迁移）
struct BinData {
    int totalEvents = 0;
    int hitEvents = 0;
    std::vector<double> resX_values;
    std::vector<double> resY_values;
};

// DUT统计辅助类（从AnalysisEngine迁移）
class DUTStatistics {
   public:
    DUTStatistics(const DUTAnalysisConfig& config) : m_config(config) {}

    // 获取bin索引
    std::pair<int, int> GetBinIndices(double predX, double predY) const;

    // 添加统计数据
    void AddBinData(int dutID, int binX, int binY,
                    bool hasValidHit, double resX, double resY);

    // 获取所有分区数据
    const std::map<int, std::map<std::pair<int, int>, BinData>>& GetBinDataMap() const {
        return m_binDataMap;
    }

    // 清空统计数据
    void Clear() { m_binDataMap.clear(); }

   private:
    DUTAnalysisConfig m_config;
    std::map<int, std::map<std::pair<int, int>, BinData>> m_binDataMap;
};

/**
 * @brief DUT Analysis Script
 * 
 * 执行DUT探测器分析，基于Track信息计算残差和效率，
 * 生成DUTInfo.root输出文件
 */
class DUTAnalysisScript : public IScript {
public:
    DUTAnalysisScript() = default;
    ~DUTAnalysisScript() override = default;

    std::string GetName() const override { return "DUTAnalysisScript"; }
    
    std::string GetDescription() const override { 
        return "DUT efficiency and residual analysis"; 
    }

    void LoadConfig(const json& config) override;
    void Print() const override;
    bool Execute() override;

private:
    // 配置参数
    bool m_runAlignment;
    bool m_saveNoiseData;
    int m_progressInterval;
    int m_maxEvents;

    // 私有静态方法
    static Cluster CreateInvalidCluster(int type);
};

#endif // DUT_ANALYSIS_SCRIPT_H
