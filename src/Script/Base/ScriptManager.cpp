#include "Script/Base/ScriptManager.h"
#include "Detector/DetectorFactory.h"
#include "Script/Base/RawDataParser.h"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

ScriptManager::ScriptManager(const std::string& configFile, const std::string& rawDir,
                             const std::string& resultDir, const std::string& runID)
    : m_configFile(configFile), m_rawDir(rawDir), m_resultDir(resultDir), m_runID(runID) {
    // 构造输出目录
    m_outputDir = m_resultDir + "/" + m_runID + "/";
}

bool ScriptManager::Initialize() {
    std::cout << "\n========================================" << std::endl;
    std::cout << "Initializing BeamAnalysis" << std::endl;
    std::cout << "========================================" << std::endl;

    // 读取配置
    std::ifstream in(m_configFile);
    if (!in.is_open()) {
        std::cerr << "Cannot open config: " << m_configFile << std::endl;
        return false;
    }
    in >> m_config;

    if (!m_config.contains("detectors")) {
        std::cerr << "No detectors in config" << std::endl;
        return false;
    }

    // 使用DetectorFactory创建探测器
    auto& factory = DetectorFactory::GetInstance();
    if (!factory.Initialize(m_config)) {
        std::cerr << "Failed to initialize DetectorFactory" << std::endl;
        return false;
    }

    // A canonical input can be selected explicitly. Keep the historical
    // run<ID>.root convention as the fallback.
    std::string rawFile = m_rawDir + "/run" + m_runID + ".root";
    if (m_config.contains("input") && m_config["input"].contains("file")) {
        rawFile = m_config["input"]["file"].get<std::string>();
        const std::string token = "{run_id}";
        if (const auto position = rawFile.find(token); position != std::string::npos)
            rawFile.replace(position, token.size(), m_runID);
    }

    std::cout << "Config file: " << m_configFile << std::endl;
    std::cout << "Raw file   : " << rawFile << std::endl;
    std::cout << "Output dir : " << m_outputDir << std::endl;
    std::cout << "Run ID     : " << m_runID << std::endl;

    // 创建输出目录
    std::filesystem::create_directories(m_outputDir);

    // 初始化parser
    m_parser = std::make_shared<RawDataParser>(rawFile);
    if (!m_parser->Initialize()) {
        std::cerr << "Failed to initialize parser" << std::endl;
        return false;
    }

    if (m_config.contains("scripts")) {
        std::cout << "Loading custom scripts..." << std::endl;
        LoadScripts(m_config["scripts"]);
    }

    std::cout << "Initialization complete" << std::endl;

    return true;
}

void ScriptManager::LoadScripts(const json& scriptsConfig) {
    if (scriptsConfig.is_null() || !scriptsConfig.is_array()) {
        std::cout << "No scripts configuration found or invalid format." << std::endl;
        return;
    }

    m_scripts.clear();

    for (const auto& scriptConfig : scriptsConfig) {
        try {
            ScriptInfo info;

            // 读取必填字段
            if (!scriptConfig.contains("name") || !scriptConfig.contains("type")) {
                std::cerr << "Script configuration missing required fields (name or type). Skipping." << std::endl;
                continue;
            }

            info.name = scriptConfig["name"].get<std::string>();
            info.type = scriptConfig["type"].get<std::string>();

            // 读取可选字段
            info.enabled = scriptConfig.value("enabled", true);
            info.config = scriptConfig.value("config", json::object());

            // 检查脚本类型是否已注册
            if (!ScriptFactory::Instance().IsRegistered(info.type)) {
                std::cerr << "Script type '" << info.type << "' (name: " << info.name
                          << ") is not registered. Skipping." << std::endl;
                continue;
            }

            // 创建脚本实例
            info.instance = ScriptFactory::Instance().CreateScript(info.type, info.config);
            if (info.instance) {
                // 注入资源
                info.instance->SetParser(m_parser);
                info.instance->SetConfig(m_config);
                info.instance->SetOutputDir(m_outputDir);

                // 验证配置
                if (!info.instance->Validate()) {
                    std::cerr << "Script '" << info.name << "' configuration validation failed. Skipping." << std::endl;
                    continue;
                }
            }

            std::cout << "Script loaded: " << info.name << " (type: " << info.type
                      << ", enabled: " << (info.enabled ? "yes" : "no") << ")" << std::endl;
            m_scripts.push_back(std::move(info));

        } catch (const std::exception& e) {
            std::cerr << "Error loading script: " << e.what() << std::endl;
        }
    }

    std::cout << "Total scripts loaded: " << m_scripts.size()
              << " (enabled: " << GetEnabledScriptCount() << ")" << std::endl;
}

std::vector<ScriptInfo> ScriptManager::GetEnabledScripts() const {
    std::vector<ScriptInfo> enabledScripts;

    for (const auto& script : m_scripts) {
        if (script.enabled) {
            enabledScripts.push_back(script);
        }
    }

    return enabledScripts;
}

bool ScriptManager::ExecuteScript(const std::string& name) {
    for (auto& script : m_scripts) {
        if (script.name == name) {
            if (!script.enabled) {
                std::cerr << "Script '" << name << "' is disabled." << std::endl;
                return false;
            }

            if (!script.instance) {
                std::cerr << "Script '" << name << "' instance is null." << std::endl;
                return false;
            }

            try {
                std::cout << "\n=== Executing Script: " << script.name << " ===" << std::endl;
                script.instance->Print();

                bool success = script.instance->Execute();

                if (success) {
                    std::cout << "=== Script '" << script.name << "' completed successfully ===" << std::endl;
                } else {
                    std::cerr << "=== Script '" << script.name << "' failed ===" << std::endl;
                }

                return success;

            } catch (const std::exception& e) {
                std::cerr << "Exception while executing script '" << name << "': "
                          << e.what() << std::endl;
                return false;
            }
        }
    }

    std::cerr << "Script '" << name << "' not found." << std::endl;
    return false;
}

int ScriptManager::ExecuteAllEnabled() {
    int successCount = 0;
    auto enabledScripts = GetEnabledScripts();

    std::cout << "\nExecuting " << enabledScripts.size() << " enabled script(s)..." << std::endl;

    for (const auto& script : enabledScripts) {
        if (ExecuteScript(script.name)) {
            successCount++;
        }
    }

    std::cout << "\nScript execution summary: " << successCount << "/"
              << enabledScripts.size() << " succeeded." << std::endl;

    return successCount;
}

size_t ScriptManager::GetEnabledScriptCount() const {
    size_t count = 0;
    for (const auto& script : m_scripts) {
        if (script.enabled) {
            count++;
        }
    }
    return count;
}
