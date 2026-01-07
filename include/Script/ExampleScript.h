#pragma once

#ifndef EXAMPLE_SCRIPT_H
#define EXAMPLE_SCRIPT_H

#include "Script/Base/IScript.h"
#include <string>

/**
 * @brief 示例脚本
 * 
 * 演示如何创建自定义分析脚本
 * 展示基本的数据访问和输出功能
 */
class ExampleScript : public IScript {
public:
    ExampleScript() = default;
    ~ExampleScript() override = default;

    std::string GetName() const override { return "ExampleScript"; }
    
    std::string GetDescription() const override { 
        return "Example script demonstrating the script framework"; 
    }

    void LoadConfig(const json& config) override;
    void Print() const override;
    bool Execute() override;

private:
    std::string m_message;
    int m_eventLimit;
};

#endif // EXAMPLE_SCRIPT_H
