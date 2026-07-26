#pragma once

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <nlohmann/json.hpp>
#include <stdexcept>
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
    double cfdFraction = 0.5;
    double noiseThreshold = 100.0;
    double saturationLevel = 2000.0;
    double timePitch = 25.0;
    double timeWindowStart = -50;
    double timeWindowEnd = 500.0;

    // 实现AlgorithmConfig接口
    void loadFrom(const json& config) override {
        const json* cfg = &config;
        if (cfg->contains("mode")) mode = (*cfg)["mode"];
        if (cfg->contains("cfdFraction")) cfdFraction = (*cfg)["cfdFraction"];
        if (cfg->contains("noiseThreshold")) noiseThreshold = (*cfg)["noiseThreshold"];
        if (cfg->contains("saturationLevel")) saturationLevel = (*cfg)["saturationLevel"];
        if (cfg->contains("timePitch")) timePitch = (*cfg)["timePitch"];
        if (cfg->contains("timeWindowStart")) timeWindowStart = (*cfg)["timeWindowStart"];
        if (cfg->contains("timeWindowEnd")) timeWindowEnd = (*cfg)["timeWindowEnd"];
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
    int minClusterSize = 1;      // 最小聚类大小
    int maxClusterSize = 10;     // 最大聚类大小
    double MaxChargeDiff = 0.4;  // 对于不同Cluster匹配最大电荷差
    int connectivity = 4;        // pad 聚类邻接方式：4 或 8

    void loadFrom(const json& config) override {
        const json* cfg = &config;
        if (cfg->contains("maxGap")) maxGap = (*cfg)["maxGap"];
        if (cfg->contains("minClusterSize")) minClusterSize = (*cfg)["minClusterSize"];
        if (cfg->contains("maxClusterSize")) maxClusterSize = (*cfg)["maxClusterSize"];
        if (cfg->contains("MaxChargeDiff")) MaxChargeDiff = (*cfg)["MaxChargeDiff"];
        if (cfg->contains("connectivity")) connectivity = (*cfg)["connectivity"];
        if (connectivity != 4 && connectivity != 8) {
            throw std::runtime_error("ClusterConfig.connectivity must be 4 or 8");
        }
    }

    void print() const override {
        std::cout << "ClusterConfig:" << std::endl;
        std::cout << "  Max Gap: " << maxGap << std::endl;
        std::cout << "  Min Cluster Size: " << minClusterSize << std::endl;
        std::cout << "  Max Cluster Size: " << maxClusterSize << std::endl;
        std::cout << "  Max Charge Diff: " << MaxChargeDiff << std::endl;
        std::cout << "  Pad Connectivity: " << connectivity << std::endl;
    }
};

enum class ReconstructionMethod {
    ChargeWeighted,  // 电荷加权
    PadChargeWeighted,  // pad 二维电荷加权
    UTPC,            // UTPC算法
    RawUTPC          // raw UTPC算法
};

// ========== 重建配置 ==========
class ReconstructionConfig : public AlgorithmConfig {
   public:
    ReconstructionMethod method = ReconstructionMethod::ChargeWeighted;
    std::string weightSource = "amp";
    double weightPower = 1.0;

    void loadFrom(const json& config) override {
        const json* cfg = &config;

        if (cfg->contains("method")) {
            std::string methodStr = (*cfg)["method"];
            if (methodStr == "UTPC") {
                method = ReconstructionMethod::UTPC;
            } else if (methodStr == "rawUTPC") {
                method = ReconstructionMethod::RawUTPC;
            } else if (methodStr == "PadChargeWeighted") {
                method = ReconstructionMethod::PadChargeWeighted;
            } else {
                method = ReconstructionMethod::ChargeWeighted;
            }
        }
        if (cfg->contains("weightSource"))
            weightSource = (*cfg)["weightSource"].get<std::string>();
        if (cfg->contains("weightPower"))
            weightPower = (*cfg)["weightPower"].get<double>();
        if (weightSource != "charge" && weightSource != "amp")
            throw std::runtime_error(
                "ReconstructionConfig.weightSource must be charge or amp");
        if (!std::isfinite(weightPower) || weightPower < 0.0)
            throw std::runtime_error(
                "ReconstructionConfig.weightPower must be finite and non-negative");
    }

    void print() const override {
        std::cout << "ReconstructionConfig:" << std::endl;
        const char* name = "ChargeWeighted";
        if (method == ReconstructionMethod::PadChargeWeighted) name = "PadChargeWeighted";
        else if (method == ReconstructionMethod::UTPC) name = "UTPC";
        else if (method == ReconstructionMethod::RawUTPC) name = "rawUTPC";
        std::cout << "  Method: " << name << std::endl;
        if (method == ReconstructionMethod::ChargeWeighted ||
            method == ReconstructionMethod::PadChargeWeighted) {
            std::cout << "  Weight Source: " << weightSource << std::endl;
            std::cout << "  Weight Power: " << weightPower << std::endl;
        }
    }
};

struct planarConfig {
    std::vector<int> readoutPlaneType = {0, 1};
    std::map<int, double> readoutPlaneAngle = {{0, 0}, {1, 90}};
    std::map<int, double> readoutPlanePitch = {{0, 0.4}, {1, 0.4}};
    std::map<int, int> readoutPlaneStripNumber = {{0, 256}, {1, 256}};
};

struct planarPadConfig {
    int rows = 0;
    int columns = 0;
    double pitchX = 0.0;
    double pitchY = 0.0;
    double sizeX = 0.0;
    double sizeY = 0.0;
    int planeType = 0;
};

struct cylinderConfig {
    double radius = 65;
    std::map<int, double> readoutPlaneAngle = {{0, 0}};
    std::map<int, double> readoutPlanePitch = {{0, 0.4}};
};
