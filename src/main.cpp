#include "AnalysisEngine.h"
#include "EventDisplay/EventDisplayManager.h"
#include <iostream>
#include <string>

void printUsage(const char* prog) {
    std::cout << "Usage: " << prog << " <run_id> [config_file]" << std::endl;
    std::cout << "  run_id     : Run ID (e.g., 1813)" << std::endl;
    std::cout << "  config_file: Config path (default: config/config.json)" << std::endl;
    std::cout << "\nExamples:" << std::endl;
    std::cout << "  " << prog << " 1813" << std::endl;
    std::cout << "  " << prog << " 1813 config/config1813.json" << std::endl;
    std::cout << "\nInput/Output:" << std::endl;
    std::cout << "  Raw data : raw/run<run_id>.root" << std::endl;
    std::cout << "  Results  : result/<run_id>/" << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        printUsage(argv[0]);
        return 1;
    }

    std::string runID = argv[1];
    std::string configFile = argc > 2 ? argv[2] : "config.json";
    configFile = "/home/qxhuang/fs/ustcfs/workarea/BeamResult_202511/config/" + configFile;
    std::string rawDir = "/home/qxhuang/fs/ustcfs/workarea/BeamResult_202511/raw";
    std::string resultDir = "/home/qxhuang/fs/ustcfs/workarea/BeamResult_202511/result";

    try {
        AnalysisEngine engine(configFile, rawDir, resultDir, runID);
        engine.Initialize();

        while (true) {
            std::cout << "\nSelect mode:" << std::endl;
            std::cout << "  1) Track Analysis" << std::endl;
            std::cout << "  2) DUT Analysis" << std::endl;
            std::cout << "  3) Event Display Mode (DUT only)" << std::endl;
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
                engine.RunTrackAnalysis();
            } else if (choice == 2) {
                engine.RunDUTAnalysis();
            } else if (choice == 3) {
                EventDisplayManager edm(rawDir, resultDir, runID);
                if (!edm.Initialize()) {
                    std::cerr << "Failed to init EventDisplayManager\n";
                    continue;  // 不退出，回到菜单
                }
                edm.RunInteractive();
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