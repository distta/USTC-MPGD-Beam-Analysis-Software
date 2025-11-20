#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>

#include "Detector.h"

using json = nlohmann::json;

class DetectorFactory {
public:
    /**
     * @brief 获取工厂单例实例
     * @return DetectorFactory引用
     */
    static DetectorFactory& GetInstance();

    /**
     * @brief 从JSON配置初始化所有探测器
     * @param config JSON配置对象,包含detectors数组
     * @return 是否初始化成功
     */
    bool Initialize(const json& config);

    /**
     * @brief 根据ID获取探测器
     * @param id 探测器ID
     * @return 探测器共享指针,如果不存在则返回nullptr
     */
    std::shared_ptr<Detector> GetDetector(int id) const;

    /**
     * @brief 获取所有探测器
     * @return ID到探测器的映射表
     */
    const std::map<int, std::shared_ptr<Detector>>& GetAllDetectors() const;

    /**
     * @brief 按角色筛选探测器
     * @param role 探测器角色(Tracker/DUT/Ignored)
     * @return 符合角色的探测器共享指针列表
     */
    std::vector<std::shared_ptr<Detector>> GetDetectorsByRole(Detector::Role role) const;

    /**
     * @brief 清空所有探测器实例
     */
    void Clear();

    /**
     * @brief 检查是否已初始化
     * @return 是否已初始化
     */
    bool IsInitialized() const { return m_initialized; }

    // 禁用拷贝构造和赋值
    DetectorFactory(const DetectorFactory&) = delete;
    DetectorFactory& operator=(const DetectorFactory&) = delete;

private:
    DetectorFactory() = default;
    ~DetectorFactory() = default;

    /**
     * @brief 创建单个探测器实例
     * @param detConfig 探测器配置JSON对象
     * @return 探测器共享指针
     * @throws std::runtime_error 配置错误或不支持的类型
     */
    std::shared_ptr<Detector> CreateDetector(const json& detConfig);

    std::map<int, std::shared_ptr<Detector>> m_detectors;  ///< 存储所有探测器实例
    bool m_initialized{false};  ///< 标记是否已初始化
};
