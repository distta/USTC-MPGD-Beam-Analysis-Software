#pragma once

#ifndef ISCRIPT_H
#define ISCRIPT_H

#include <memory>
#include <nlohmann/json.hpp>
#include <string>

using json = nlohmann::json;

// 前向声明，避免循环依赖
class RawDataParser;

/**
 * @brief 脚本基础接口
 *
 * 所有自定义分析脚本必须继承此接口
 * 提供基本的元信息、配置加载和执行能力
 */
class IScript {
   public:
    virtual ~IScript() = default;

    /**
     * @brief 获取脚本名称
     * @return 脚本类型标识字符串
     */
    virtual std::string GetName() const = 0;

    /**
     * @brief 获取脚本版本
     * @return 版本字符串，默认 1.0.0
     */
    virtual std::string GetVersion() const { return "1.0.0"; }

    /**
     * @brief 获取脚本描述
     * @return 脚本功能描述字符串
     */
    virtual std::string GetDescription() const { return ""; }

    /**
     * @brief 加载配置
     * @param config JSON配置对象，内部结构由具体脚本自行定义和解析
     */
    virtual void LoadConfig(const json& config) = 0;

    /**
     * @brief 打印脚本配置信息
     */
    virtual void Print() const = 0;

    /**
     * @brief 执行脚本主逻辑
     * @return 执行是否成功
     */
    virtual bool Execute() = 0;

    /**
     * @brief 设置数据解析器
     * @param parser 数据解析器指针
     */
    void SetParser(std::shared_ptr<RawDataParser> parser) { m_parser = parser; }

    /**
     * @brief 获取数据解析器
     * @return 数据解析器指针，如果未设置则返回nullptr
     */
    std::shared_ptr<RawDataParser> GetParser() const { return m_parser; }

    /**
     * @brief 设置全局配置
     * @param config 全局配置JSON对象
     */
    void SetConfig(const json& config) { m_config = config; }

    /**
     * @brief 获取全局配置
     * @return 全局配置JSON对象引用
     */
    const json& GetConfig() const { return m_config; }

    /**
     * @brief 设置输出目录
     * @param outputDir 输出目录路径
     */
    void SetOutputDir(const std::string& outputDir) { m_outputDir = outputDir; }

    /**
     * @brief 获取输出目录
     * @return 输出目录路径
     */
    std::string GetOutputDir() const { return m_outputDir; }

    /**
     * @brief 初始化钩子（预留扩展接口）
     * @return 初始化是否成功
     */
    virtual bool Initialize() { return true; }

    /**
     * @brief 清理钩子（预留扩展接口）
     */
    virtual void Finalize() {}

    /**
     * @brief 配置验证接口（预留扩展接口）
     * @return 配置是否有效
     */
    virtual bool Validate() const { return true; }

   protected:
    std::shared_ptr<RawDataParser> m_parser;  ///< 数据解析器指针
    json m_config;                            ///< 全局配置对象
    std::string m_outputDir;                  ///< 输出目录路径
};

/**
 * @brief 脚本创建函数类型
 */
using ScriptCreateFunc = std::shared_ptr<IScript> (*)();

#endif  // ISCRIPT_H
