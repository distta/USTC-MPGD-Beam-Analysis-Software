#pragma once


#include <memory>
#include <string>
#include <vector>

#include "Detector/Detector.h"
#include "RawDataParser.h"
#include <nlohmann/json.hpp>

// DUT分析配置类
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
            : predX_min(0), predX_max(100), nBinsX(5),
              predY_min(0), predY_max(100), nBinsY(5) {}
    };
    
    BinningConfig binning;
    
    // 判断位置是否在有效范围内
    bool IsInValidRange(double predX, double predY) const {
        return predX >= binning.predX_min && predX <= binning.predX_max &&
               predY >= binning.predY_min && predY <= binning.predY_max;
    }
};

// 分区统计数据
struct BinData {
    int totalEvents = 0;
    int hitEvents = 0;
    std::vector<double> resX_values;
    std::vector<double> resY_values;
};

// DUT统计辅助类
class DUTStatistics {
public:
    DUTStatistics(const DUTAnalysisConfig& config) : m_config(config) {}
    
    // 计算均值
    static double CalculateMean(const std::vector<double>& values);
    
    // 计算RMS
    static double CalculateRMS(const std::vector<double>& values);
    
    // 获取bin索引
    std::pair<int, int> GetBinIndices(double predX, double predY) const;
    
    // 添加统计数据
    void AddBinData(int dutID, int binX, int binY, 
                    bool hasValidHit, double resX, double resY);
    
    // 获取所有分区数据
    const std::map<int, std::map<std::pair<int,int>, BinData>>& GetBinDataMap() const {
        return m_binDataMap;
    }
    
    // 清空统计数据
    void Clear() { m_binDataMap.clear(); }
    
private:
    DUTAnalysisConfig m_config;
    std::map<int, std::map<std::pair<int,int>, BinData>> m_binDataMap;
};

using json = nlohmann::json;

Track FitTrack(const std::vector<TVector3>& globalHits);

class AnalysisEngine {
   public:
    AnalysisEngine(const std::string& configFile,
                   const std::string& rawDir,
                   const std::string& resultDir,
                   const std::string& runID);

    void Initialize();

    void RunTrackAnalysis();

    void RunDUTAnalysis();
    void RunDUTAlign();

   private:
    std::shared_ptr<RawDataParser> m_parser;

    std::string m_configFile;
    std::string m_rawDir;
    std::string m_resultDir;
    std::string m_runID;
    std::string m_outputDir;  // result/runID/

    json m_config;
    std::vector<Event> m_events;

    // DUT对齐相关私有方法
    double DUTChi2Objective(
        const double* par,
        const std::vector<Event>& events,
        std::shared_ptr<Detector> detector,
        int detID);

    LocalHit CalcuDutResidual(
        std::shared_ptr<Detector> detector,
        const std::vector<Cluster>& clusters,
        const TVector3& predL,
        double& residualX, double& residualY);
    
    // DUT分析辅助函数
    Cluster CreateInvalidCluster(int type);
};
