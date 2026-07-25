#pragma once

#ifndef DUT_ANALYSIS_SCRIPT_H
#define DUT_ANALYSIS_SCRIPT_H

#include "Event/DataModel.h"
#include "Script/DUTEfficiencyAnalysis.h"
#include "Script/Base/IScript.h"
#include <string>
#include <vector>

class DUTAnalysisConfig {
   public:
    static constexpr int kTypeX = 0;
    static constexpr int kTypeY = 1;
    static constexpr double kInvalidValue = -999.0;
    static constexpr int kInvalidSize = -1;
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
    int m_progressInterval;
    int m_maxEvents;
    int m_effXBins = 8;
    int m_effYBins = 7;
    int m_effMinEntriesPerBin = 1;
    std::vector<int> m_effExcludedXBins;
    std::vector<int> m_effExcludedYBins;
    double m_effXMin = 0.0;
    double m_effXMax = 100.0;
    double m_effYMin = 0.0;
    double m_effYMax = 100.0;
    DUTEfficiency::Config m_efficiencyConfig;

    // 私有方法
    static Cluster CreateInvalidCluster(int type);

};

#endif  // DUT_ANALYSIS_SCRIPT_H
