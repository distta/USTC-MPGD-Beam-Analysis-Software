#pragma once

#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "Detector/Detector.h"

class TFile;

/**
 * @brief 径迹分析模块
 *
 * 负责径迹重建、探测器对齐、残差计算等功能
 * 从 AnalysisEngine 中分离出来,专注于径迹分析
 */
class TrackAnalysis {
   public:
    /**
     * @brief 构造函数
     * @param outputDir 输出目录路径
     * @param runID 运行ID
     */
    TrackAnalysis(const std::string& outputDir);

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
     * @brief 获取探测器本征分辨率映射表
     * @return 分辨率映射表
     */
    const std::map<int, std::pair<double, double>>& GetSigmaMap() const {
        return m_sigmaMap;
    }

    /**
     * @brief 设置探测器本征分辨率映射表
     * @param sigmaMap 分辨率映射表
     */
    void SetSigmaMap(const std::map<int, std::pair<double, double>>& sigmaMap) {
        m_sigmaMap = sigmaMap;
    }

   private:
    /**
     * @brief 计算预测误差
     * @param targetDetID 目标探测器ID
     * @return X/Y方向预测误差
     */
    std::pair<double, double> ComputePredictionError(int targetDetID);

    std::string m_outputDir;  ///< 结果输出目录

    std::vector<int> m_trackerIDs;                        ///< 所有Tracker探测器的ID列表
    std::map<int, std::pair<double, double>> m_sigmaMap;  ///< 探测器本征分辨率缓存
    std::vector<int> m_seedTrackerIDs;                    ///< 当前事件的种子探测器ID
};
