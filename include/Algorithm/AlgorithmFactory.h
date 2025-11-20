#pragma once

#ifndef ALGORITHM_FACTORY_H
#define ALGORITHM_FACTORY_H

#include "IAlgorithm.h"
#include <map>
#include <string>
#include <memory>
#include <functional>
#include <mutex>
#include <vector>

/**
 * @brief 算法工厂类（单例模式）
 * 
 * 负责算法的注册和创建
 * 线程安全，支持在程序启动时通过静态对象自动注册算法
 */
class AlgorithmFactory {
public:
    /**
     * @brief 获取工厂单例实例
     */
    static AlgorithmFactory& Instance();

    /**
     * @brief 注册算法
     * @param name 算法类型名称
     * @param createFunc 创建函数
     */
    void RegisterAlgorithm(const std::string& name, AlgorithmCreateFunc createFunc);

    /**
     * @brief 创建算法实例
     * @param name 算法类型名称
     * @param config 配置JSON对象
     * @return 算法实例智能指针
     * @throws std::runtime_error 如果算法未注册
     */
    std::shared_ptr<IAlgorithm> CreateAlgorithm(const std::string& name, const json& config);

    /**
     * @brief 检查算法是否已注册
     * @param name 算法类型名称
     * @return true如果已注册，否则false
     */
    bool IsRegistered(const std::string& name) const;

    /**
     * @brief 获取所有已注册的算法名称
     * @return 算法名称列表
     */
    std::vector<std::string> GetRegisteredNames() const;

private:
    AlgorithmFactory() = default;
    ~AlgorithmFactory() = default;
    
    // 禁止拷贝和赋值
    AlgorithmFactory(const AlgorithmFactory&) = delete;
    AlgorithmFactory& operator=(const AlgorithmFactory&) = delete;

    std::unordered_map<std::string, AlgorithmCreateFunc> m_registry;
    mutable std::mutex m_mutex;  // 线程安全
};

/**
 * @brief 算法注册辅助类
 * 
 * 通过静态对象的构造函数自动注册算法
 * 使用方法：在算法实现文件中定义静态AlgorithmRegistrar对象
 */
class AlgorithmRegistrar {
public:
    AlgorithmRegistrar(const std::string& name, AlgorithmCreateFunc createFunc) {
        AlgorithmFactory::Instance().RegisterAlgorithm(name, createFunc);
    }
};

/**
 * @brief 算法注册宏
 * 
 * 简化算法注册过程
 * 使用方法：REGISTER_ALGORITHM("algorithmName", AlgorithmClass)
 */
#define REGISTER_ALGORITHM(name, AlgorithmClass) \
    namespace { \
        std::shared_ptr<IAlgorithm> Create##AlgorithmClass() { \
            return std::make_shared<AlgorithmClass>(); \
        } \
        static AlgorithmRegistrar registrar_##AlgorithmClass(name, Create##AlgorithmClass); \
    }

#endif // ALGORITHM_FACTORY_H
