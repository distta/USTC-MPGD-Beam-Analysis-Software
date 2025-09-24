

#pragma once

#include <map>
#include <memory>
#include <string>

#include "DataModel.h"
#include "Detector/Detector.h"
#include "Detector/Planar.h"
#include <nlohmann/json.hpp>

using json = nlohmann::json;

class Pipeline {
   public:
    explicit Pipeline(const std::string& configFile);
    ~Pipeline() = default;

    /**
     * Run 完整流程
     *  rawFile   : 输入的原始 ROOT 文件（含 raw tree）
     *  cacheFile : heavy 阶段缓存文件路径（若存在则跳过 heavy 阶段）
     *  outFile   : 最终输出文件（包含 Local/Global/Track）
     */
    void Run(const std::string& rawFile, const std::string& cacheFile, const std::string& outFile);

   private:
    void InitializeDetectors();

    void RunClustering(const std::string& rawFile, const std::string& cacheFile);

    void RunTracking(const std::string& cacheFile, const std::string& outFile);

    // 把硬件 (boardID, channelID) 映射到 (detID, stripID, type)
    std::tuple<int, int, int> MapBoardChannel(unsigned int boardID, unsigned int channelID) const;

    bool EventFilter(const std::vector<RecHit>&);

    // 简单轨迹拟合（直线拟合 x(z), y(z)）
    Track FitTrack(const std::map<int, std::vector<GlobalHit>>& globalHits) const;

   private:
    std::string m_configFile;
    std::string m_rawDataFileName;

    json m_config;

    // detID -> shared_ptr<Detector>
    std::map<int, std::shared_ptr<Detector>> m_dets;
};
