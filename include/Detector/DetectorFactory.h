#pragma once

#include <map>
#include <memory>
#include <nlohmann/json.hpp>
#include <string>
#include <vector>

#include "Detector.h"

using json = nlohmann::json;

class DetectorFactory {
   public:
    /**
     * @brief 获取工厂单例实例
     * @return DetectorFactory引用
     */
    static DetectorFactory& GetInstance();

    bool Initialize(const json& config);

    std::shared_ptr<Detector> GetDetector(int id) const;

    const std::map<int, std::shared_ptr<Detector>>& GetAllDetectors() const;

    std::vector<std::shared_ptr<Detector>> GetDetectorsByRole(Detector::Role role) const;

    std::vector<int> GetDetectorIDsByRole(Detector::Role role) const;

    void Clear();

    // 禁用拷贝构造和赋值
    DetectorFactory(const DetectorFactory&) = delete;
    DetectorFactory& operator=(const DetectorFactory&) = delete;

   private:
    DetectorFactory() = default;
    ~DetectorFactory() = default;

    std::shared_ptr<Detector> CreateDetector(const json& detConfig);

    std::map<int, std::shared_ptr<Detector>> m_detectors;  ///< 存储所有探测器实例
};
