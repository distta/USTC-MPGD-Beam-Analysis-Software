#pragma once

#ifndef DUT_ANALYSIS_SCRIPT_H
#define DUT_ANALYSIS_SCRIPT_H

#include "Event/DataModel.h"
#include "Script/Base/IScript.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"

#include <map>
#include <string>
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
    bool m_saveEfficiencyMap = true;
    int m_progressInterval;
    int m_maxEvents;
    int m_effXBins = 8;
    int m_effYBins = 7;
    int m_effMinEntriesPerBin = 1;
    std::vector<int> m_effExcludedXBins;
    std::vector<int> m_effExcludedYBins;
    double m_effXMin = 10.0;
    double m_effXMax = 90.0;
    double m_effYMin = 10.0;
    double m_effYMax = 45.0;
    double m_effMatchRadius = 2.0;

    // 私有方法
    static Cluster CreateInvalidCluster(int type);

    // 设置残差直方图格式
    void SetResidualHistogramStyle(TH1D* hist, TF1* fitFunc, const std::string& axis);
};

#endif  // DUT_ANALYSIS_SCRIPT_H
