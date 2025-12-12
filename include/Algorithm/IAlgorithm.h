#pragma once

#ifndef IALGORITHM_H
#define IALGORITHM_H

#include <memory>
#include <nlohmann/json.hpp>
#include <string>

using json = nlohmann::json;

// 前向声明，避免循环依赖
class Detector;
class DetectorFrame;

/**
 * @brief 算法基础接口
 *
 * 所有算法必须继承此接口，提供基本的元信息和配置加载能力
 * 不强制具体的处理接口，各算法根据需要自行定义功能接口
 */
class IAlgorithm {
   public:
    virtual ~IAlgorithm() = default;

    /**
     * @brief 获取算法名称
     * @return 算法名称字符串
     */
    virtual std::string GetName() const = 0;

    /**
     * @brief 获取算法版本
     * @return 版本字符串
     */
    virtual std::string GetVersion() const { return "1.0.0"; }

    /**
     * @brief 加载配置
     * @param config JSON配置对象，内部结构由具体算法自行定义和解析
     *
     * 统一的配置加载接口，实现概念层面的统一
     * 各算法内部可以自由定义config的结构
     */
    virtual void LoadConfig(const json& config) = 0;

    /**
     * @brief 打印算法配置信息
     */
    virtual void Print() const = 0;

    /**
     * @brief 设置拥有此算法的探测器
     * @param detector 探测器指针
     */
    void SetDetector(Detector* detector) { m_detector = detector; }

    /**
     * @brief 获取拥有此算法的探测器
     * @return 探测器指针，如果未设置则返回nullptr
     */
    Detector* GetDetector() const { return m_detector; }

    /**
     * @brief 统一的算法处理接口
     * @param frame DetectorFrame的引用，算法可通过此对象访问和修改数据
     * @return 处理是否成功
     * 
     * 各算法子类根据自身职责实现该方法:
     * - WaveformProcessor: 读取Raw数据，生成StripHit
     * - ClusterBuilder: 读取StripHit，生成Cluster
     * - ClusterReconstructor: 读取StripHit和Cluster，重建Cluster位置
     */
    virtual bool Process(DetectorFrame& frame) { return true; }

   protected:
    Detector* m_detector = nullptr;  ///< 拥有此算法的探测器指针
};

/**
 * @brief 算法创建函数类型
 */
using AlgorithmCreateFunc = std::shared_ptr<IAlgorithm> (*)();

#endif  // IALGORITHM_H
