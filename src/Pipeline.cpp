#include "Pipeline.h"
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

void Pipeline::Initialize(const std::string& configFile) {
   std::ifstream in(configFile);
   if (!in.is_open()) throw std::runtime_error("Failed to open config file!");

   json config;
   in >> config;

   for (auto& detConfig : config["detectors"]) {
      int detID = detConfig["id"];
      std::string name = detConfig["name"];
      auto detector = std::make_shared<Detector>(detID, name);
      detector->initialize(detConfig);  // 初始化几何+算法
      AddDetector(detector);
   }
}

void Pipeline::Run(const std::string& inputFile) {
   rawChain_ = new TChain("raw");
   rawChain_->Add(inputFile.c_str());

   // 设置分支
   apv_id_ = new std::vector<unsigned int>;
   apv_ch_ = new std::vector<unsigned int>;
   mm_strip_ = new std::vector<unsigned int>;
   apv_q_ = new std::vector<std::vector<double>>;

   rawChain_->SetBranchAddress("apv_evt", &apv_evt_);
   rawChain_->SetBranchAddress("apv_id", &apv_id_);
   rawChain_->SetBranchAddress("apv_ch", &apv_ch_);
   rawChain_->SetBranchAddress("mm_strip", &mm_strip_);
   rawChain_->SetBranchAddress("apv_q", &apv_q_);

   int totalEntries = rawChain_->GetEntries();
   std::cout << "Processing file: " << inputFile
             << ", total events: " << totalEntries << std::endl;

   // ---- 核心事件循环 ----
   for (int i = 0; i < totalEntries; ++i) {
      rawChain_->GetEntry(i);
      if (i % 1000 == 0) {
         std::cout << "Event " << i << " / " << totalEntries << std::endl;
      }

      // 1. 构建RawData
      std::map<int, std::vector<RawData>> rawDataBuffer;
      for (size_t j = 0; j < apv_id_->size(); j++) {
         int boardID = apv_id_->at(j);
         int channelID = apv_ch_->at(j);

         auto [detID, stripID, type] = ElectronicMap(boardID, channelID);

         RawData raw;
         raw.StripID = stripID;
         raw.type = type;
         raw.adc = apv_q_->at(j);

         rawDataBuffer[detID].push_back(raw);
      }

      // 2. 逐探测器重建
      Event event;
      event.eventID = apv_evt_;

      for (auto& [detID, det] : detectors_) {
         auto& rawHits = rawDataBuffer[detID];

         // Raw -> StripHit -> Cluster
         auto stripHits = det->ProcessRawData(rawHits);

         auto recClusters = det->BuildClusters(stripHits);
         event.clusterMap[detID] = recClusters;

         // Cluster -> LocalHit
         auto localHits = det->RecLocalHit(recClusters);
         event.localHits[detID] = localHits;

         // Local -> Global
         for (auto& lh : localHits) {
            event.globalHits[detID].push_back(det->LocalToGlobal(lh));
         }
      }

      // 3. 简单Track拟合（可直接内联，不单独函数）
      if (event.globalHits.size() >= 2) {
      }

      // 4. 输出或保存
      std::cout << "Event " << event.eventID
                << " hits: " << event.globalHits.size() << std::endl;
   }
}

std::tuple<int, int, int> Pipeline::ElectronicMap(int boardID, int channelID) {
   // 映射规则示例
   int detID = boardID / 10;  // 例如：10x = detectorID
   int stripID = channelID % 128;
   int type = (channelID < 128) ? 0 : 1;  // 0 = X, 1 = Y
   return {detID, stripID, type};
}