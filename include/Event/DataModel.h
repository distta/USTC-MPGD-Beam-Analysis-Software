#pragma once

#include "TVector3.h"
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <vector>

class DetectorFrame;

struct RawData {
    int id0;                 // strip ID，或 pad column
    int id1{-1};             // pad row；strip 读出为 -1
    int type;                // 读出面编号
    std::vector<short> adc;  // 波形采样点 (ADC counts)
    int chip{-1};            // 前端芯片编号
    int channel{-1};         // 芯片通道编号
    unsigned int rawHitID{0};
    double adcValue{std::numeric_limits<double>::quiet_NaN()};
    double hitTimeNs{std::numeric_limits<double>::quiet_NaN()};

    bool HasID1() const { return id1 >= 0; }
    bool HasDirectHit() const {
        return std::isfinite(adcValue) && std::isfinite(hitTimeNs);
    }
};

struct ChannelHit {
    int id0;                 // strip ID，或 pad column
    int id1{-1};             // pad row；strip 读出为 -1
    int type;                // 读出面/通道类型
    int rawIndices;

    bool HasID1() const { return id1 >= 0; }

    // ---- 提取的关键量 ----
    double amp;       // 峰值 - baseline
    double charge;    // 积分电荷
    int peakTime;     // 峰值时间 (采样点索引)
    double time;      // 信号时间 (前沿拟合+CFD)
    double riseTime;  // 上升时间 (10%-90%)
    double width;     // 信号宽度

    // ---- 误差信息 ----
    double timeError;  // 时间误差

    // ---- 标志位 ----
    bool isSaturated;  // 是否饱和
    bool isValid;      // 是否有效
};

struct Cluster {
    int type;                            // 所属读出面/通道类型
    std::vector<int> channelHitIndices;  // ChannelHit在DetectorFrame::m_channelHits中的全局索引

    // ---- 聚类整体量 ----
    int size;       // 聚类通道数
    int range;      // 一维通道范围；二维 pad 聚类算法不应依赖该字段
    double charge;  // 聚类总电荷
    double maxAmp;
    double time;      // 聚类时间 (最早)
    double centroid;  // 聚类重心
    double pos;       // cluster重建位置
    TVector3 localPosition;  // 二维读出重建位置
    bool hasLocalPosition{false};
};

typedef TVector3 GlobalHit;

struct LocalHit {
    TVector3 localPos;
    std::vector<int> clusterIndices;
};

// 径迹数据结构
struct Track {
    double kx;    // x方向斜率
    double ky;    // y方向斜率
    double bx;    // x截距
    double by;    // y截距
    double chi2;  // 拟合质量
};

struct Event {
    int eventID;
    std::map<int, std::shared_ptr<DetectorFrame>> detectorFramesMap;
    Track track;
    double t0;
};
