#pragma once

#include "Algorithm/Config.h"
#include "TVector3.h"
#include <iostream>
#include <map>
#include <nlohmann/json.hpp>

using json = nlohmann::json;


struct RawData {
    int stripID;             // 读出板编号
    int type;                // 通道编号
    std::vector<short> adc;  // 波形采样点 (ADC counts)
};

struct StripHit {
    int stripID;  // 条号
    int type;     // 读出条类型 (X, Y, U/V…)

    // ---- 提取的关键量 ----
    double amp;       // 峰值 - baseline
    double charge;    // 积分电荷
    int peakTime;     // 峰值时间 (采样点索引)
    double time;      // 信号时间 (前沿拟合+CFD)
    double riseTime;  // 上升时间 (10%-90%)

    // ---- 误差信息 ----
    double timeError;  // 时间误差

    // ---- 标志位 ----
    bool isSaturated;  // 是否饱和
    bool isValid;      // 是否有效
};

struct Cluster {
    int type;  // 属于X/Y/U/V 哪个方向
    int matchID;
    std::vector<StripHit> strips;  // 参与聚类的条

    // ---- 聚类整体量 ----
    int size;       // 聚类条数
    int range;      // 聚类范围 (最大ID - 最小ID)
    double charge;  // 聚类总电荷
    double maxAmp;
    double time;  // 聚类时间 (最早)
    double pos;   // cluster重建位置
};

typedef TVector3 GlobalHit;
typedef TVector3 LocalHit;
typedef std::vector<Cluster> RecCluster;

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
    std::map<int, std::vector<LocalHit>> recLocalHits;     // detID -> LocalHits
    std::map<int, std::vector<RecCluster>> recClusters;    // detID -> RecClusters
    Track track;
};
 