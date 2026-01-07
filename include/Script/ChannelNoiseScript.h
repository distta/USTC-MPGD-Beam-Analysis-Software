#pragma once

#ifndef CHANNEL_NOISE_SCRIPT_H
#define CHANNEL_NOISE_SCRIPT_H

#include "Script/Base/IScript.h"
#include <string>
#include <map>
#include <vector>

/**
 * @brief 探测器通道噪声分析脚本
 * 
 * 分析每个探测器通道的噪声特性，包括：
 * - 基线（baseline）统计
 * - 噪声（noise）统计
 * - 通道噪声分布
 * - 异常通道识别
 * 
 * 输出ROOT文件包含：
 * - 每个通道的噪声直方图
 * - 每个探测器的通道噪声分布图
 * - 通道噪声统计树
 */
class ChannelNoiseScript : public IScript {
public:
    ChannelNoiseScript() = default;
    ~ChannelNoiseScript() override = default;

    std::string GetName() const override { return "ChannelNoiseScript"; }
    
    std::string GetVersion() const override { return "1.0.0"; }
    
    std::string GetDescription() const override { 
        return "Analyze detector channel noise characteristics and output to ROOT file"; 
    }

    void LoadConfig(const json& config) override;
    void Print() const override;
    bool Execute() override;

private:
    /**
     * @brief 通道噪声统计数据结构
     */
    struct ChannelNoise {
        int detectorID;
        int stripID;
        int type;  // 0=X, 1=Y
        
        double meanBaseline;     // 平均基线
        double rmsBaseline;      // 基线RMS
        double meanNoise;        // 平均噪声（采样点标准差）
        double maxNoise;         // 最大噪声
        int nSamples;            // 采样事例数
        
        bool isNoisy;            // 是否为噪声通道
        bool isDead;             // 是否为死通道
        
        std::vector<double> baselineList;  // 保存原始基线数据用于绘制直方图
    };

    /**
     * @brief 分析所有事例的通道噪声
     */
    void AnalyzeChannelNoise();

    /**
     * @brief 保存结果到ROOT文件
     */
    bool SaveResults();

    // 配置参数
    int m_maxEvents;           // 分析的最大事例数（0表示全部）
    int m_progressInterval;    // 进度显示间隔
    double m_noisyThreshold;   // 噪声通道阈值（ADC）
    double m_deadThreshold;    // 死通道阈值（ADC）

    // 分析结果
    std::map<int, std::map<int, std::map<int, ChannelNoise>>> m_channelNoiseMap;  
    // [detectorID][type][stripID] -> ChannelNoise
};

#endif // CHANNEL_NOISE_SCRIPT_H
