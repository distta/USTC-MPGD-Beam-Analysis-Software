#pragma once

#include <nlohmann/json.hpp>
#include <string>

using json = nlohmann::json;

// 波形处理配置
struct WaveformConfig {
    double cfdFraction = 0.1;
    double noiseThreshold = 60.0;
    double saturationLevel = 1800.0;
    double timePitch = 25.0;

    void loadFrom(const json& config, const std::string& section = "waveform") {
        if (!config.contains(section)) return;

        const auto& cfg = config[section];
        if (cfg.contains("cfdFraction")) cfdFraction = cfg["cfdFraction"];
        if (cfg.contains("noiseThreshold")) noiseThreshold = cfg["noiseThreshold"];
        if (cfg.contains("saturationLevel")) saturationLevel = cfg["saturationLevel"];
        if (cfg.contains("timePitch")) timePitch = cfg["timePitch"];
    }
};

// 聚类配置
struct ClusterConfig {
    int maxGap = 1;           // 最大间隙
    int minClusterSize = 2;   // 最小聚类大小
    int maxClusterSize = 20;  // 最大聚类大小

    void loadFrom(const json& config, const std::string& section = "cluster") {
        if (!config.contains(section)) return;

        const auto& cfg = config[section];
        if (cfg.contains("maxGap")) maxGap = cfg["maxGap"];
        if (cfg.contains("minClusterSize")) minClusterSize = cfg["minClusterSize"];
        if (cfg.contains("maxClusterSize")) maxClusterSize = cfg["maxClusterSize"];
    }
};

enum class ReconstructionMethod {
    ChargeWeighted,  // 电荷加权
    UTPC,            // UTPC算法
};

// 重建配置
struct ReconstructionConfig {
    ReconstructionMethod method = ReconstructionMethod::ChargeWeighted;

    void loadFrom(const json& config, const std::string& section = "reconstruction") {
        if (!config.contains(section)) return;

        const auto& cfg = config[section];
        if (cfg.contains("method")) {
            std::string methodStr = cfg["method"];
            if (methodStr == "UTPC") {
                method = ReconstructionMethod::UTPC;
            }
        } else {
            method = ReconstructionMethod::ChargeWeighted;
        }
    }
};

struct planarConfig {
    std::map<int, double> readoutPlaneAngle = {{0, 0}, {1, 90}};
    std::map<int, double> readoutPlanePitch = {{0, 0.4}, {1, 0.4}};
};

struct cylinderConfig {
    double radius = 65;
    std::map<int, double> readoutPlaneAngle = {{0, 0}};
    std::map<int, double> readoutPlanePitch = {{0, 0.4}};
};
