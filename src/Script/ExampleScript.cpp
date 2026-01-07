#include "Script/ExampleScript.h"
#include "Script/Base/ScriptFactory.h"
#include <fstream>
#include <iostream>

void ExampleScript::LoadConfig(const json& config) {
    m_message = config.value("message", "Hello from ExampleScript!");
    m_eventLimit = config.value("eventLimit", 10);
}

void ExampleScript::Print() const {
    std::cout << "ExampleScript Configuration:" << std::endl;
    std::cout << "  Message: " << m_message << std::endl;
    std::cout << "  Event Limit: " << m_eventLimit << std::endl;
}

bool ExampleScript::Execute() {
    std::cout << "\n"
              << m_message << std::endl;

    // 获取配置信息
    const auto& config = GetConfig();
    std::cout << "Config loaded successfully" << std::endl;

    // 创建输出文件
    std::string outputFile = GetOutputDir() + "example_output.txt";
    std::ofstream outFile(outputFile);
    if (!outFile.is_open()) {
        std::cerr << "Error: Cannot create output file: " << outputFile << std::endl;
        return false;
    }

    outFile << "Example Script Output" << std::endl;
    outFile << "=====================" << std::endl;
    outFile << "Message: " << m_message << std::endl;
    outFile << "Output Directory: " << GetOutputDir() << std::endl;
    outFile << std::endl;

    outFile.close();
    std::cout << "Output saved to: " << outputFile << std::endl;

    return true;
}

// 注册脚本
REGISTER_SCRIPT("ExampleScript", ExampleScript);
