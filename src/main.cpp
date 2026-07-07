#include "Event/EventDisplayManager.h"
#include "Input/ConverterFactory.h"
#include "Script/Base/ScriptManager.h"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <string>

using json = nlohmann::json;
namespace fs = std::filesystem;

namespace {

struct RunPaths {
    std::string runID;
    std::string profile = "config/defaults.json";
    std::string configFile;
    std::string baseDir;
    std::string configDir;
    std::string rawDir;
    std::string convertedDir;
    std::string resultDir;
    bool dryRun = false;
    bool help = false;
};

const json kBuiltInDefaults = {
    {"baseDir", "/home/qxhuang/fs/ustcfs/workarea/BeamResult_202511"},
    {"configDir", "{baseDir}/config"},
    {"rawDir", "{baseDir}/raw"},
    {"convertedDir", "{baseDir}/converted"},
    {"resultDir", "{baseDir}/result"},
    {"defaultConfig", "config.json"}};

std::string Expand(std::string value, const RunPaths& paths) {
    const std::pair<std::string, std::string> tokens[] = {
        {"{baseDir}", paths.baseDir},
        {"{run_id}", paths.runID}};
    for (const auto& [token, replacement] : tokens) {
        size_t pos = 0;
        while ((pos = value.find(token, pos)) != std::string::npos) {
            value.replace(pos, token.size(), replacement);
            pos += replacement.size();
        }
    }
    return value;
}

std::string JoinPath(const std::string& directory, const std::string& name) {
    return (fs::path(directory) / name).lexically_normal().string();
}

std::string DisplayPath(const std::string& path) {
    if (path.empty()) return "(empty)";
    return fs::absolute(fs::path(path)).lexically_normal().string();
}

std::string NextArg(int& i, int argc, char* argv[], const std::string& option) {
    if (i + 1 >= argc) throw std::runtime_error("Missing value for " + option);
    return argv[++i];
}

void LoadProfileDefaults(RunPaths& paths, json& defaults) {
    if (!fs::exists(paths.profile)) return;
    std::ifstream input(paths.profile);
    if (!input) throw std::runtime_error("Cannot open profile: " + paths.profile);
    json profileDefaults;
    input >> profileDefaults;
    defaults.update(profileDefaults);
}

RunPaths ParseRunOptions(int argc, char* argv[]) {
    RunPaths paths;
    paths.runID = argv[1];

    bool positionalConfigUsed = false;
    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        const auto equal = arg.find('=');
        std::string value;
        if (equal != std::string::npos) {
            value = arg.substr(equal + 1);
            arg = arg.substr(0, equal);
        }

        auto optionValue = [&] { return value.empty() ? NextArg(i, argc, argv, arg) : value; };

        if (arg == "-h" || arg == "--help") {
            paths.help = true;
        } else if (arg == "--dry-run") {
            paths.dryRun = true;
        } else if (arg == "--profile") {
            paths.profile = optionValue();
        } else if (arg == "-c" || arg == "--config") {
            paths.configFile = optionValue();
        } else if (arg == "--base-dir") {
            paths.baseDir = optionValue();
        } else if (arg == "--config-dir") {
            paths.configDir = optionValue();
        } else if (arg == "--raw-dir") {
            paths.rawDir = optionValue();
        } else if (arg == "--converted-dir") {
            paths.convertedDir = optionValue();
        } else if (arg == "--result-dir") {
            paths.resultDir = optionValue();
        } else if (!arg.empty() && arg[0] != '-' && !positionalConfigUsed && paths.configFile.empty()) {
            paths.configFile = arg;
            positionalConfigUsed = true;
        } else {
            throw std::runtime_error("Unknown argument: " + arg);
        }
    }

    json defaults = kBuiltInDefaults;
    LoadProfileDefaults(paths, defaults);

    auto pick = [&](const std::string& value, const char* key) {
        return Expand(value.empty() ? defaults.value(key, std::string()) : value, paths);
    };

    paths.baseDir = pick(paths.baseDir, "baseDir");
    paths.configDir = pick(paths.configDir, "configDir");
    paths.rawDir = pick(paths.rawDir, "rawDir");
    paths.convertedDir = pick(paths.convertedDir, "convertedDir");
    paths.resultDir = pick(paths.resultDir, "resultDir");
    paths.configFile = pick(paths.configFile, "defaultConfig");

    fs::path configPath(paths.configFile);
    if (!configPath.is_absolute() && !fs::exists(configPath)) {
        configPath = fs::path(paths.configDir) / configPath;
    }
    paths.configFile = configPath.lexically_normal().string();

    return paths;
}

void PrintRuntimePaths(const RunPaths& paths) {
    std::cout << "\nResolved run setup" << std::endl;
    std::cout << "  Run ID      : " << paths.runID << std::endl;
    std::cout << "  Config file : " << DisplayPath(paths.configFile) << std::endl;
    std::cout << "  Raw dir     : " << DisplayPath(paths.rawDir) << std::endl;
    std::cout << "  Converted   : " << DisplayPath(paths.convertedDir) << std::endl;
    std::cout << "  Result dir  : " << DisplayPath(paths.resultDir) << std::endl;
    std::cout << "  Mode        : "
              << (paths.dryRun ? "dry-run" : "interactive") << std::endl;
}

bool CheckFileExists(const std::string& label, const std::string& path) {
    const bool ok = fs::exists(path);
    std::cout << "  " << label << ": " << (ok ? "OK " : "MISSING ")
              << DisplayPath(path) << std::endl;
    return ok;
}

bool CheckInputIfNeeded(const std::string& input, const std::string& output) {
    const bool inputExists = CheckFileExists("raw input", input);
    return fs::exists(output) || inputExists;
}

bool RunDryRun(const RunPaths& paths) {
    std::cout << "\nDry-run checks" << std::endl;
    bool ok = CheckFileExists("config", paths.configFile);
    if (!ok) return false;

    std::ifstream input(paths.configFile);
    json config;
    input >> config;

    if (!config.contains("detectors") || !config["detectors"].is_array()) {
        std::cerr << "  detectors: MISSING or invalid" << std::endl;
        return false;
    }
    std::cout << "  detectors: " << config["detectors"].size() << std::endl;

    const std::string rawFile = JoinPath(paths.rawDir, "run" + paths.runID + ".root");
    std::string canonicalFile = JoinPath(paths.convertedDir, "run" + paths.runID + ".root");
    bool conversionCanCreateCanonical = false;
    if (config.contains("input") && config["input"].is_object() &&
        config["input"].contains("file")) {
        canonicalFile = Expand(config["input"]["file"].get<std::string>(), paths);
    }

    if (config.contains("conversion") && config["conversion"].is_object()) {
        json conversion = config["conversion"];
        if (!conversion.contains("input")) conversion["input"] = rawFile;
        if (!conversion.contains("output")) conversion["output"] = canonicalFile;
        const std::string type = conversion.value("type", "SRS");
        const std::string conversionInput = Expand(conversion.value("input", rawFile), paths);
        canonicalFile = Expand(conversion.value("output", canonicalFile), paths);
        const bool conversionEnabled = conversion.value("enabled", true);
        std::cout << "  conversion : "
                  << (conversionEnabled ? type : "disabled")
                  << std::endl;
        if (conversionEnabled) {
            if (!ConverterFactory::Create(type)) {
                std::cerr << "  converter : UNKNOWN " << type << std::endl;
                ok = false;
            }
            conversionCanCreateCanonical = true;
            ok = CheckInputIfNeeded(conversionInput, canonicalFile) && ok;
        }
    } else if (!config.contains("input")) {
        std::cout << "  conversion : automatic SRS fallback" << std::endl;
        conversionCanCreateCanonical = true;
        ok = CheckInputIfNeeded(rawFile, canonicalFile) && ok;
    }

    const bool canonicalExists = CheckFileExists("canonical", canonicalFile);
    if (!conversionCanCreateCanonical && !canonicalExists) {
        ok = false;
    }

    if (config.contains("scripts") && config["scripts"].is_array()) {
        size_t enabled = 0;
        for (const auto& script : config["scripts"]) {
            if (script.value("enabled", true)) ++enabled;
        }
        std::cout << "  scripts   : " << enabled << " enabled / "
                  << config["scripts"].size() << " configured" << std::endl;
    } else {
        std::cout << "  scripts   : none configured" << std::endl;
    }

    std::cout << "\nDry-run complete. No analysis was executed." << std::endl;
    return ok;
}

}  // namespace

void printUsage(const char* prog) {
    std::cout << "Usage: " << prog << " <run_id> [config_file] [options]" << std::endl;
    std::cout << "       " << prog << " convert <conversion_config> [output.root]" << std::endl;
    std::cout << "  run_id     : Run ID (e.g., 1813)" << std::endl;
    std::cout << "  config_file: Optional config name/path (default from config/defaults.json)" << std::endl;
    std::cout << "\nOptions:" << std::endl;
    std::cout << "  --profile <file>       Defaults profile (default: config/defaults.json)" << std::endl;
    std::cout << "  -c, --config <file>    Config file or name under configDir" << std::endl;
    std::cout << "  --base-dir <dir>       Dataset base directory" << std::endl;
    std::cout << "  --config-dir <dir>     Directory for config names" << std::endl;
    std::cout << "  --raw-dir <dir>        Raw input directory" << std::endl;
    std::cout << "  --converted-dir <dir>  Canonical ROOT directory" << std::endl;
    std::cout << "  --result-dir <dir>     Result output directory" << std::endl;
    std::cout << "  --dry-run              Resolve pa ths and validate config only" << std::endl;
    std::cout << "  -h, --help             Show this help" << std::endl;
    std::cout << "\nExamples:" << std::endl;
    std::cout << "  " << prog << " 1813" << std::endl;
    std::cout << "  " << prog << " 1813 config1813.json" << std::endl;
    std::cout << "  " << prog << " 1813 --base-dir /path/to/BeamResult" << std::endl;
    std::cout << "  " << prog << " 1813 --dry-run" << std::endl;
    std::cout << "\nInput/Output:" << std::endl;
    std::cout << "  Raw data : {baseDir}/raw/run<run_id>.root" << std::endl;
    std::cout << "  Converted: {baseDir}/converted/run<run_id>.root" << std::endl;
    std::cout << "  Results  : {baseDir}/result/<run_id>/" << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc < 2 || std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help") {
        printUsage(argv[0]);
        return argc < 2 ? 1 : 0;
    }

    try {
        const RunPaths paths = ParseRunOptions(argc, argv);
        if (paths.help) {
            printUsage(argv[0]);
            return 0;
        }
        PrintRuntimePaths(paths);

        if (paths.dryRun) {
            return RunDryRun(paths) ? 0 : 1;
        }

        // 创建ScriptManager，传入配置参数
        ScriptManager scriptManager(paths.configFile, paths.rawDir,
                                    paths.convertedDir, paths.resultDir,
                                    paths.runID);

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
                EventDisplayManager edm(paths.convertedDir, paths.resultDir, paths.runID);
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
