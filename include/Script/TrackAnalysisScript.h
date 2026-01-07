#pragma once

#ifndef TRACK_ANALYSIS_SCRIPT_H
#define TRACK_ANALYSIS_SCRIPT_H

#include "Script/Base/IScript.h"
#include "Event/DataModel.h"
#include <string>
#include <map>
#include <vector>
#include <tuple>

class TFile;

/**
 * @brief Track Analysis Script
 * 
 * 执行径迹分析，处理Tracker探测器数据，
 * 重建轨道并生成TrackInfo.root输出文件
 */
class TrackAnalysisScript : public IScript {
public:
    TrackAnalysisScript() = default;
    ~TrackAnalysisScript() override = default;

    std::string GetName() const override { return "TrackAnalysisScript"; }
    
    std::string GetDescription() const override { 
        return "Track reconstruction and validation analysis"; 
    }

    void LoadConfig(const json& config) override;
    void Print() const override;
    bool Execute() override;

private:
    // 配置参数
    bool m_saveValidationData;
    int m_progressInterval;

    // 从TrackAnalysis迁移的成员
    std::vector<int> m_trackerIDs;                        ///< 所有Tracker探测器的ID列表
    std::map<int, std::pair<double, double>> m_sigmaMap;  ///< 探测器本征分辨率缓存
    std::vector<int> m_seedTrackerIDs;                    ///< 当前事件的种子探测器ID

    // 从TrackAnalysis迁移的方法
    /**
     * @brief 执行径迹对齐流程
     * @param events 事件列表
     * @param file 输出ROOT文件指针
     */
    void RunTrackerAlign(const std::vector<Event>& events, TFile* file);

    /**
     * @brief 计算径迹残差并输出直方图
     * @param events 事件列表
     * @param file 输出ROOT文件指针
     * @return 探测器本征分辨率映射表 (detID -> (sigmaX, sigmaY))
     */
    std::map<int, std::pair<double, double>> ComputeTrackError(
        const std::vector<Event>& events,
        TFile* file);

    /**
     * @brief 执行探测器精对齐
     * @param events 事件列表
     */
    void AlignTrackers(const std::vector<Event>& events);

    /**
     * @brief 为单个事件寻找最佳径迹
     * @param event 单个事件数据
     * @return 元组: (径迹, 击中索引映射, 是否成功)
     */
    std::tuple<Track, std::map<int, int>, bool> FindBestTrack(const Event& event);

    /**
     * @brief 计算预测误差
     * @param targetDetID 目标探测器ID
     * @return X/Y方向预测误差
     */
    std::pair<double, double> ComputePredictionError(int targetDetID);
};

#endif // TRACK_ANALYSIS_SCRIPT_H
