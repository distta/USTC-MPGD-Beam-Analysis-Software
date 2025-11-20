#pragma once

#include "DataModel.h"
#include "RawDataParser.h"
#include "Detector/DetectorFactory.h"
#include "TCanvas.h"

#include <string>
#include <vector>
#include <map>

class EventDisplayManager {
public:
    // rawDir: 原始 root 所在目录（会找到 rawDir + "/run" + runID + ".root"）
    // resultDir: 结果目录（TrackInfo.root 在 resultDir/runID/TrackInfo.root）
    EventDisplayManager(const std::string& rawDir,
                        const std::string& resultDir,
                        const std::string& runID);

    ~EventDisplayManager();

    // 初始化 RawDataParser & trackInfo 路径等
    bool Initialize();

    // 交互式运行：列出 track entries，选择 entry -> 绘制并保存所有内容
    void RunInteractive();

private:
    struct TrackEntry {
        int eventID;
        Track track;
    };

    // 读取 TrackInfo.root 中所有 Tracks
    bool LoadTrackEntries();

    // 对单个 entry 的完整流程：LoadEvent -> Reconstruct DUTs -> 绘图 -> 保存
    bool ProcessEntry(const TrackEntry& te);

    // 显示 DUT 概览（按类型画幅度/时间，inset local XY，标出 nearest cluster 与 predicted hit）
    void DrawDUTOverview(int eventID, std::shared_ptr<Detector> det, const Track& track);

    // 在同一画布中为一个 cluster 绘制所有 strip 的波形（所有 ADC 采样点）
    void DrawClusterWaveforms(int eventID, std::shared_ptr<Detector> det, const Cluster& cluster, TDirectory* dir);

    // 寻找与 predicted local hit 最近的 cluster（返回指向 cluster 的指针或 nullptr）
    const Cluster* FindNearestCluster(const std::vector<RecCluster>& recClusters,
                                      const LocalHit& predLocal,
                                      double& outDist,
                                      int& outClusterType);

    // 帮助：把 stripID -> RawData 对应起来（det->GetRawData() 中查找）
    const RawData* FindRawForStrip(const std::vector<RawData>& raw, int stripID, int type) const;

private:
    std::string m_rawDir;
    std::string m_resultDir;
    std::string m_runID;

    std::string m_rawFilePath;    // rawDir + "/run" + runID + ".root"
    std::string m_trackFilePath;  // resultDir + "/" + runID + "/TrackInfo.root"

    std::unique_ptr<RawDataParser> m_parser;

    std::vector<TrackEntry> m_trackEntries;
    std::vector<TCanvas*> m_canvases; // 管理 ROOT canvas

    // 输出文件基名（每个 event 会写入一个 root）
    std::string m_outBaseDir;
};
