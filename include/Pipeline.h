#pragma once

#include "Detector.h"
#include "DataModel.h"
#include <TChain.h>
#include <memory>
#include <map>
#include <string>
#include <vector>

class Pipeline {
public:
    Pipeline() = default;
    ~Pipeline() = default;

    // 配置初始化
    void Initialize(const std::string& configFile);

    // 添加探测器
    void AddDetector(std::shared_ptr<Detector> det) {
        detectors_[det->GetID()] = det;
    }

    // 直接运行主循环
    void Run(const std::string& inputFile);

private:
    std::tuple<int, int, int> ElectronicMap(int boardID, int channelID);

private:
    std::map<int, std::shared_ptr<Detector>> detectors_;

    // ROOT 变量
    TChain* rawChain_ = nullptr;
    unsigned int apv_evt_ = 0;
    std::vector<unsigned int>* apv_id_ = nullptr;
    std::vector<unsigned int>* apv_ch_ = nullptr;
    std::vector<unsigned int>* mm_strip_ = nullptr;
    std::vector<std::vector<double>>* apv_q_ = nullptr;
};