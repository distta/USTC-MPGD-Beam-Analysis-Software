// #include "ClusterFinder.h"
// #include "DetectorTypes.h"
// #include "PreProcessor.h"
// #include "TrackReconstructor.h"
// #include <filesystem>
// #include <iostream>

// int main(int argc, char* argv[]) {
//    if (argc < 2) {
//       std::cerr << "Usage: " << argv[0] << " <run_number>" << std::endl;
//       return 1;
//    }

//    int runNumber = std::stoi(argv[1]);
//    ConfigManager::Instance().LoadConfig("config/config.yaml");

//    // 设置当前运行号
//    ConfigManager::Instance().SetCurrentRun(runNumber);

//    // 获取文件路径
//    std::string rawPath = ConfigManager::Instance().GetRawDataPath();
//    std::string preprocessedPath = ConfigManager::Instance().GetPreprocessedPath();
//    std::string clustersPath = ConfigManager::Instance().GetClustersPath();
//    std::string tracksPath = ConfigManager::Instance().GetTracksPath();

//    // 检查并执行预处理
//    if (!std::filesystem::exists(preprocessedPath)) {
//       std::cout << "Running preprocessing for run " << runNumber << std::endl;
//       PreProcessor preprocessor;
//       preprocessor.ProcessAndSave(rawPath, preprocessedPath);
//    }

//    // 检查并执行聚类查找
//    if (!std::filesystem::exists(clustersPath)) {
//       std::cout << "Finding clusters for run " << runNumber << std::endl;
//       ClusterFinder clusterFinder;
//       clusterFinder.FindAndSaveClusters(preprocessedPath, clustersPath);
//    }

//    // 检查并执行径迹重建
//    if (!std::filesystem::exists(tracksPath)) {
//       std::cout << "Reconstructing tracks for run " << runNumber << std::endl;
//       TrackReconstructor reconstructor;
//       reconstructor.ReconstructAndSave(clustersPath, tracksPath);
//    }

//    std::cout << "Analysis completed for run " << runNumber << std::endl;
//    return 0;
// }

#include "Analysis.h"
int main(int argc, char* argv[]) {
   Analysis analysis("/home/qxhuang/workarea/Beam/config/config.json", 1591);
   analysis.initialize();
}