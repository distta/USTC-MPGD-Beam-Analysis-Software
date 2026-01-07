#pragma once

#ifndef SCRIPT_MANAGER_H
#define SCRIPT_MANAGER_H

#include "Script/Base/IScript.h"
#include "Script/Base/ScriptFactory.h"
#include <memory>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

// 前向声明
class RawDataParser;

/**
 * @brief 脚本信息结构
 */
struct ScriptInfo {
    std::string name;                  ///< 脚本显示名称
    std::string type;                  ///< 脚本类型（对应注册名称）
    bool enabled;                      ///< 是否启用
    json config;                       ///< 脚本配置
    std::shared_ptr<IScript> instance; ///< 脚本实例
};

/**
 * @brief 脚本管理器
 * 
 * 负责脚本的加载、配置和执行
 * 管理所有脚本实例的生命周期
 * 同时负责资源初始化和注入
 */
class ScriptManager {
public:
    /**
     * @brief 构造函数
     * @param configFile 配置文件路径
     * @param rawDir 原始数据目录
     * @param resultDir 结果输出目录
     * @param runID 运行ID
     */
    ScriptManager(const std::string& configFile, const std::string& rawDir,
                  const std::string& resultDir, const std::string& runID);
    ~ScriptManager() = default;

    /**
     * @brief 初始化资源
     * 读取配置、初始化DetectorFactory、创建Parser、创建输出目录
     * @return 初始化是否成功
     */
    bool Initialize();

    /**
     * @brief 从配置加载脚本定义
     * @param scriptsConfig scripts配置节点的JSON数组
     */
    void LoadScripts(const json& scriptsConfig);

    /**
     * @brief 获取所有已启用的脚本信息
     * @return 已启用的脚本信息列表
     */
    std::vector<ScriptInfo> GetEnabledScripts() const;

    /**
     * @brief 根据显示名称执行脚本
     * @param name 脚本显示名称
     * @return 执行是否成功
     */
    bool ExecuteScript(const std::string& name);

    /**
     * @brief 执行所有已启用的脚本
     * @return 成功执行的脚本数量
     */
    int ExecuteAllEnabled();

    /**
     * @brief 获取脚本数量
     * @return 脚本总数
     */
    size_t GetScriptCount() const { return m_scripts.size(); }

    /**
     * @brief 获取已启用脚本数量
     * @return 已启用的脚本数量
     */
    size_t GetEnabledScriptCount() const;

    /**
     * @brief 清空所有脚本
     */
    void Clear() { m_scripts.clear(); }

private:
    std::vector<ScriptInfo> m_scripts;  ///< 脚本列表

    // 资源管理
    std::shared_ptr<RawDataParser> m_parser;  ///< 数据解析器实例
    json m_config;                            ///< 全局配置对象
    std::string m_configFile;                 ///< 配置文件路径
    std::string m_rawDir;                     ///< 原始数据目录
    std::string m_resultDir;                  ///< 结果输出目录
    std::string m_runID;                      ///< 运行ID
    std::string m_outputDir;                  ///< 完整输出路径
};

#endif // SCRIPT_MANAGER_H
