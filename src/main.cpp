/**
 * @file main.cpp
 * @brief BeamAnalysis 项目主程序入口
 * @author Huang Qixuan
 */

#include "Pipeline.h"
#include <iostream>
#include <string>
#include <filesystem>
#include <sstream>

void printUsage(const char* programName) {
    std::cout << "用法: " << programName << " <run_number> [config_file]" << std::endl;
    std::cout << "参数说明:" << std::endl;
    std::cout << "  run_number  : 实验运行编号 (例如: 1537, 1538, ...)" << std::endl;
    std::cout << "  config_file : 配置文件路径 (可选, 默认: config/config.json)" << std::endl;
    std::cout << "示例:" << std::endl;
    std::cout << "  " << programName << " 1537" << std::endl;
    std::cout << "  " << programName << " 1538 config/custom.json" << std::endl;
    std::cout << "\n输入输出规则:" << std::endl;
    std::cout << "  输入: raw/run<run_number>.root" << std::endl;
    std::cout << "  输出: result/<run_number>/" << std::endl;
}

std::string formatDataFileName(const std::string& runNumber) {
    // 构造数据文件名: raw/run<runNumber>.root
    std::stringstream ss;
    ss << "raw/run" << runNumber << ".root";
    return ss.str();
}

std::string formatOutputDir(const std::string& runNumber) {
  
    std::stringstream ss;
    ss << "result/" << runNumber;
    return ss.str();
}

int main(int argc, char* argv[]) {
    try {

        if (argc < 2) {
            std::cerr << "错误: 缺少必要参数" << std::endl;
            printUsage(argv[0]);
            return 1;
        }
        
        std::string runNumber = argv[1];
        std::string configFile = "config/config.json";

        if (argc > 2) {
            configFile = argv[2];
        }
        
        std::string dataFile = formatDataFileName(runNumber);
        std::string outputDir = formatOutputDir(runNumber);
        
        if (!std::filesystem::exists(dataFile)) {
            std::cerr << "错误: 数据文件不存在: " << dataFile << std::endl;
            std::cerr << "请确保 raw/ 目录下存在对应的数据文件" << std::endl;
            return 1;
        }
        
        if (!std::filesystem::exists(configFile)) {
            std::cerr << "错误: 配置文件不存在: " << configFile << std::endl;
            return 1;
        }
        
        std::cout << "BeamAnalysis 启动中..." << std::endl;
        std::cout << "配置文件: " << configFile << std::endl;
        std::cout << "数据文件: " << dataFile << std::endl;
        std::cout << "运行编号: " << runNumber << std::endl;
        std::cout << "输出目录: " << outputDir << std::endl;
        
        // 创建输出目录
        std::filesystem::create_directories(outputDir);
        // 创建并运行分析管道
        Pipeline pipeline;
        pipeline.Initialize(configFile);
        pipeline.SetRawDataFile(dataFile);  // 设置要处理的数据文件
        pipeline.SetOutputDirectory(outputDir);  // 设置输出目录
        pipeline.Run();
        pipeline.Finalize();
        
        std::cout << "分析完成！" << std::endl;
        std::cout << "结果已保存到 " << outputDir << " 目录" << std::endl;
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "错误: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "未知错误发生" << std::endl;
        return 1;
    }
}