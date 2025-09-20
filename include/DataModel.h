#pragma once

#include <map>
#include <nlohmann/json.hpp>
#include <optional>
#include "Rtypes.h"
#include <vector>

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
    double amplitude;  // 峰值 - baseline
    double charge;     // 积分电荷
    int peakTime;      // 峰值时间 (采样点索引)
    double time;       // 信号时间 (前沿拟合+CFD)
    double riseTime;   // 上升时间 (10%-90%)

    // ---- 误差信息 ----
    double timeError;  // 时间误差

    // ---- 标志位 ----
    bool isSaturated;  // 是否饱和
    bool isValid;      // 波形是否有效 (fit失败/噪声事件标记)
};

struct RecCluster {
    int type;                      // 属于X/Y/U/V 哪个方向
    std::vector<StripHit> strips;  // 参与聚类的条

    // ---- 聚类整体量 ----
    int size;       // 聚类条数
    int range;      // 聚类范围 (最大ID - 最小ID)
    double charge;  // 聚类总电荷
    double maxAmplitude;
    double time;  // 聚类时间 (最早)
    double pos;   // cluster重建位置
};

struct LocalHit {
    double u = 0.0;                          // 一维量 (始终存在)
    std::optional<double> v = std::nullopt;  // 二维探测器时有效，1D 探测器为空
};

struct GlobalHit {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

// 径迹数据结构
struct Track {

    double slope_x;      // x方向斜率
    double slope_y;      // y方向斜率
    double intercept_x;  // x截距
    double intercept_y;  // y截距
    double chi2;         // 拟合质量

    std::pair<double, double> getPositionAtZ(double z) const {
        double x = intercept_x + slope_x * z;
        double y = intercept_y + slope_y * z;
        return {x, y};
    }
};

struct Event {
    int eventID;

    std::map<int, std::vector<RecCluster>> recClusters;
    std::map<int, std::vector<LocalHit>> recLocalHits;
    std::map<int, std::vector<GlobalHit>> recGlobalHits;
    Track track;
};
