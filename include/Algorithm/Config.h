#pragma once

#include <iostream>
#include <map>
#include <nlohmann/json.hpp>
#include <string>

using json = nlohmann::json;

class AlgorithmConfig {
   public:
    virtual ~AlgorithmConfig() = default;

    virtual void loadFrom(const json& config) = 0;

    virtual void print() const = 0;
};

// ========== 波形处理配置 ==========
class WaveformConfig : public AlgorithmConfig {
   public:
    std::string mode = "Default";
    double cfdFraction = 0.1;
    double noiseThreshold = 60.0;
    double saturationLevel = 1800.0;
    double timePitch = 25.0;

    // 实现AlgorithmConfig接口
    void loadFrom(const json& config) override {
        const json* cfg = &config;
        if (cfg->contains("mode")) mode = (*cfg)["mode"];
        if (cfg->contains("cfdFraction")) cfdFraction = (*cfg)["cfdFraction"];
        if (cfg->contains("noiseThreshold")) noiseThreshold = (*cfg)["noiseThreshold"];
        if (cfg->contains("saturationLevel")) saturationLevel = (*cfg)["saturationLevel"];
        if (cfg->contains("timePitch")) timePitch = (*cfg)["timePitch"];
    }

    void print() const override {
        std::cout << "WaveformConfig:" << std::endl;
        std::cout << "  Mode:" << mode << std::endl;
        std::cout << "  CFD Fraction:" << cfdFraction << std::endl;
        std::cout << "  Noise Threshold: " << noiseThreshold << std::endl;
        std::cout << "  Saturation Level: " << saturationLevel << std::endl;
        std::cout << "  Time Pitch: " << timePitch << std::endl;
    }
};

// ========== 聚类配置 ==========
class ClusterConfig : public AlgorithmConfig {
   public:
    int maxGap = 0;              // 最大间隙
    int minClusterSize = 2;      // 最小聚类大小
    int maxClusterSize = 6;     // 最大聚类大小
    double MaxChargeDiff = 0.4;  // 对于不同Cluster匹配最大电荷差

    void loadFrom(const json& config) override {
        const json* cfg = &config;
        if (cfg->contains("maxGap")) maxGap = (*cfg)["maxGap"];
        if (cfg->contains("minClusterSize")) minClusterSize = (*cfg)["minClusterSize"];
        if (cfg->contains("maxClusterSize")) maxClusterSize = (*cfg)["maxClusterSize"];
        if (cfg->contains("MaxChargeDiff")) MaxChargeDiff = (*cfg)["MaxChargeDiff"];
    }

    void print() const override {
        std::cout << "ClusterConfig:" << std::endl;
        std::cout << "  Max Gap: " << maxGap << std::endl;
        std::cout << "  Min Cluster Size: " << minClusterSize << std::endl;
        std::cout << "  Max Cluster Size: " << maxClusterSize << std::endl;
        std::cout << "  Max Charge Diff: " << MaxChargeDiff << std::endl;
    }
};

enum class ReconstructionMethod {
    ChargeWeighted,  // 电荷加权
    UTPC,            // UTPC算法
};

// ========== 重建配置 ==========
class ReconstructionConfig : public AlgorithmConfig {
   public:
    ReconstructionMethod method = ReconstructionMethod::ChargeWeighted;

    void loadFrom(const json& config) override {
        const json* cfg = &config;

        if (cfg->contains("method")) {
            std::string methodStr = (*cfg)["method"];
            if (methodStr == "UTPC") {
                method = ReconstructionMethod::UTPC;
            } else {
                method = ReconstructionMethod::ChargeWeighted;
            }
        }
    }

    void print() const override {
        std::cout << "ReconstructionConfig:" << std::endl;
        std::cout << "  Method: " << (method == ReconstructionMethod::ChargeWeighted ? "ChargeWeighted" : "UTPC") << std::endl;
    }
};

struct planarConfig {
    std::map<int, double> readoutPlaneAngle = {{0, 0}, {1, 90}};
    std::map<int, double> readoutPlanePitch = {{0, 0.4}, {1, 0.4}};
    std::map<int, int> readoutPlaneStripNumber = {{0, 256}, {1, 256}};
};

struct cylinderConfig {
    double radius = 65;
    std::map<int, double> readoutPlaneAngle = {{0, 0}};
    std::map<int, double> readoutPlanePitch = {{0, 0.4}};
};
