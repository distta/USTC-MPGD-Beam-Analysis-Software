#include "Event/EventDisplayManager.h"
#include "Input/ConverterFactory.h"
#include "Script/Base/ScriptManager.h"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <string>

using json = nlohmann::json;

void printUsage(const char* prog) {
    std::cout << "Usage: " << prog << " <run_id> [config_file]" << std::endl;
    std::cout << "       " << prog << " convert <conversion_config> [output.root]" << std::endl;
    std::cout << "  run_id     : Run ID (e.g., 1813)" << std::endl;
    std::cout << "  config_file: Config path (default: config/config.json)" << std::endl;
    std::cout << "\nExamples:" << std::endl;
    std::cout << "  " << prog << " 1813" << std::endl;
    std::cout << "  " << prog << " 1813 config/config1813.json" << std::endl;
    std::cout << "\nInput/Output:" << std::endl;
    std::cout << "  Raw data : raw/run<run_id>.root" << std::endl;
    std::cout << "  Results  : result/<run_id>/" << std::endl;
}

int runConversion(int argc, char* argv[]) {
    if (argc < 3) {
        printUsage(argv[0]);
        return 1;
    }

    std::ifstream input(argv[2]);
    if (!input) {
        std::cerr << "Cannot open conversion config: " << argv[2] << std::endl;
        return 1;
    }
    json config;
    input >> config;
    const json& conversion = config.contains("conversion") ? config["conversion"] : config;
    const std::string type = conversion.value("type", "");
    const std::string output = argc > 3 ? argv[3] : conversion.value("output", "");
    if (type.empty() || output.empty()) {
        std::cerr << "Conversion config requires type and output" << std::endl;
        return 1;
    }

    auto converter = ConverterFactory::Create(type);
    if (!converter) {
        std::cerr << "Unknown converter type: " << type << std::endl;
        return 1;
    }
    return converter->Convert(conversion, output) ? 0 : 1;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        printUsage(argv[0]);
        return 1;
    }

    if (std::string(argv[1]) == "convert") return runConversion(argc, argv);

    std::string runID = argv[1];
    std::string configFile = argc > 2 ? argv[2] : "config.json";
    if (!std::filesystem::exists(configFile))
        configFile = "/home/qxhuang/fs/ustcfs/workarea/BeamResult_202511/config/" + configFile;
    // std::string rawDir = "/home/qxhuang/fs/ustcfs/workarea/BeamResult_202511/raw";
    std::string rawDir = "/home/qxhuang/fs/ustcfs/workarea/BeamResult_LHCb/converted";
    std::string resultDir = "/home/qxhuang/fs/ustcfs/workarea/BeamResult_LHCb/LHCB/result";

    try {
        // 创建ScriptManager，传入配置参数
        ScriptManager scriptManager(configFile, rawDir, resultDir, runID);

        // 初始化资源
        if (!scriptManager.Initialize()) {
            std::cerr << "Failed to initialize ScriptManager" << std::endl;
            return 1;
        }

        while (true) {
            std::cout << "\nSelect mode:" << std::endl;
            std::cout << "  1) Event Display Mode (DUT only)" << std::endl;
            std::cout << "  2) Run Custom Scripts" << std::endl;
            std::cout << "  0) Exit" << std::endl;
            std::cout << "Choice: ";

            int choice;
            std::cin >> choice;

            if (std::cin.fail()) {
                std::cin.clear();                                                    // 清除错误标志
                std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');  // 忽略错误输入
                std::cerr << "Invalid input! Please enter a number." << std::endl;
                continue;
            }

            if (choice == 1) {
                std::cout << "Event Display Mode not implemented yet." << std::endl;
                EventDisplayManager edm(rawDir, resultDir, runID);
                if (!edm.Initialize()) {
                    std::cerr << "Failed to init EventDisplayManager\n";
                    continue;
                }
                edm.RunInteractive();
            } else if (choice == 2) {
                // Run Custom Scripts
                std::cout << "\n========================================" << std::endl;
                std::cout << " Custom Scripts Execution" << std::endl;
                std::cout << "========================================" << std::endl;

                auto enabledScripts = scriptManager.GetEnabledScripts();

                if (enabledScripts.empty()) {
                    std::cout << "No enabled scripts found." << std::endl;
                    continue;
                }

                // 显示可用脚本列表
                std::cout << "\nAvailable Scripts:" << std::endl;
                for (size_t i = 0; i < enabledScripts.size(); ++i) {
                    const auto& script = enabledScripts[i];
                    std::cout << "  [" << (i + 1) << "] " << script.name;
                    if (!script.instance->GetDescription().empty()) {
                        std::cout << " - " << script.instance->GetDescription();
                    }
                    std::cout << std::endl;
                }
                std::cout << "  [0] Return to main menu" << std::endl;

                // 用户选择
                std::cout << "\nSelect script to execute: ";
                int scriptChoice;
                std::cin >> scriptChoice;

                if (std::cin.fail()) {
                    std::cin.clear();
                    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    std::cerr << "Invalid input!" << std::endl;
                    continue;
                }

                if (scriptChoice == 0) {
                    std::cout << "Returning to main menu..." << std::endl;
                    continue;
                }

                if (scriptChoice < 1 || scriptChoice > static_cast<int>(enabledScripts.size())) {
                    std::cerr << "Invalid choice!" << std::endl;
                    continue;
                }

                // 执行选中的脚本
                const auto& selectedScript = enabledScripts[scriptChoice - 1];
                bool success = scriptManager.ExecuteScript(selectedScript.name);

                if (success) {
                    std::cout << "\nScript execution completed successfully." << std::endl;
                } else {
                    std::cerr << "\nScript execution failed." << std::endl;
                }
            } else if (choice == 0) {
                std::cout << "Exiting program..." << std::endl;
                break;
            } else {
                std::cerr << "Invalid choice! Please try again." << std::endl;
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "\nError: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
