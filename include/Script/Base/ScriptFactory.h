#pragma once

#ifndef SCRIPT_FACTORY_H
#define SCRIPT_FACTORY_H

#include "Script/Base/IScript.h"
#include <map>
#include <string>
#include <memory>
#include <functional>
#include <mutex>
#include <vector>
#include <unordered_map>

/**
 * @brief 脚本工厂类（单例模式）
 * 
 * 负责脚本的注册和创建
 * 线程安全，支持在程序启动时通过静态对象自动注册脚本
 */
class ScriptFactory {
public:
    /**
     * @brief 获取工厂单例实例
     */
    static ScriptFactory& Instance();

    /**
     * @brief 注册脚本
     * @param name 脚本类型名称
     * @param createFunc 创建函数
     */
    void RegisterScript(const std::string& name, ScriptCreateFunc createFunc);

    /**
     * @brief 创建脚本实例
     * @param name 脚本类型名称
     * @param config 配置JSON对象
     * @return 脚本实例智能指针
     * @throws std::runtime_error 如果脚本未注册
     */
    std::shared_ptr<IScript> CreateScript(const std::string& name, const json& config);

    /**
     * @brief 检查脚本是否已注册
     * @param name 脚本类型名称
     * @return true如果已注册，否则false
     */
    bool IsRegistered(const std::string& name) const;

    /**
     * @brief 获取所有已注册的脚本名称
     * @return 脚本名称列表
     */
    std::vector<std::string> GetRegisteredNames() const;

private:
    ScriptFactory() = default;
    ~ScriptFactory() = default;
    
    // 禁止拷贝和赋值
    ScriptFactory(const ScriptFactory&) = delete;
    ScriptFactory& operator=(const ScriptFactory&) = delete;

    std::unordered_map<std::string, ScriptCreateFunc> m_registry;
    mutable std::mutex m_mutex;  // 线程安全
};

/**
 * @brief 脚本注册辅助类
 * 
 * 通过静态对象的构造函数自动注册脚本
 * 使用方法：在脚本实现文件中定义静态ScriptRegistrar对象
 */
class ScriptRegistrar {
public:
    ScriptRegistrar(const std::string& name, ScriptCreateFunc createFunc) {
        ScriptFactory::Instance().RegisterScript(name, createFunc);
    }
};

/**
 * @brief 脚本注册宏
 * 
 * 简化脚本注册过程
 * 使用方法：REGISTER_SCRIPT("scriptName", ScriptClass)
 */
#define REGISTER_SCRIPT(name, ScriptClass) \
    namespace { \
        std::shared_ptr<IScript> Create##ScriptClass() { \
            return std::make_shared<ScriptClass>(); \
        } \
        static ScriptRegistrar registrar_##ScriptClass(name, Create##ScriptClass); \
    }

#endif // SCRIPT_FACTORY_H
