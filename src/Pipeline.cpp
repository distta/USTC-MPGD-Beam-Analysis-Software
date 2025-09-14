#include "Pipeline.h"
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

Pipeline::~Pipeline() {
   delete rawChain_;
   delete apv_id_;
   delete apv_ch_;
   delete apv_q_;
}

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

std::tuple<int, int, int> Pipeline::ElectronicMap(int boardID, int channelID) {
   // 映射规则示例
   int detID = boardID / 10;  // 例如：10x = detectorID
   int stripID = channelID % 128;
   int type = (channelID < 128) ? 0 : 1;  // 0 = X, 1 = Y
   return {detID, stripID, type};
}

void Pipeline::Run(const std::string& inputFile) {
   rawChain_ = new TChain("raw");
   rawChain_->Add(inputFile.c_str());

   apv_id_ = new std::vector<unsigned int>;
   apv_ch_ = new std::vector<unsigned int>;
   apv_q_ = new std::vector<std::vector<double>>;

   rawChain_->SetBranchAddress("apv_evt", &apv_evt_);
   rawChain_->SetBranchAddress("apv_id", &apv_id_);
   rawChain_->SetBranchAddress("apv_ch", &apv_ch_);
   rawChain_->SetBranchAddress("apv_q", &apv_q_);

   int totalEntries = rawChain_->GetEntries();
   std::cout << "Processing file: " << inputFile
             << ", total events: " << totalEntries << std::endl;

   // ---- 核心事件循环 ----
   for (int i = 0; i < totalEntries; ++i) {
      rawChain_->GetEntry(i);

      // 1. 构建RawData
      std::map<int, std::vector<RawData>> rawDataBuffer;
      for (size_t j = 0; j < apv_id_->size(); j++) {
         int boardID = apv_id_->at(j);
         int channelID = apv_ch_->at(j);

         auto [detID, stripID, type] = ElectronicMap(boardID, channelID);
         RawData raw{stripID, type, apv_q_->at(j)};
         rawDataBuffer[detID].push_back(raw);
      }

      // 2. algorithms: rawData->StripHit->RecCluster->LocalHit->GlobleHit
      Event event;
      event.eventID = apv_evt_;

      for (auto& [detID, det] : m_dets) {
         auto& rawHits = rawDataBuffer[detID];
         if (rawHits.empty()) continue;

         // Raw -> StripHit
         auto stripHits = det->BuildStripHit(rawHits);

         // StripHit -> Cluster
         auto recClusters = det->BuildClusters(stripHits);
         event.clusterMap[detID] = recClusters;

         // Cluster -> LocalHit
         auto localHits = det->RecLocalHit(recClusters);
         event.localHits[detID] = localHits;

         // LocalHit -> GlobalHit
         for (auto& lh : localHits) {
            GlobalHit gh = det->LocalToGlobal(lh);
            event.globalHits[detID].push_back(gh);
         }
      }

      // 3. Find Track And Pick Event you want
      if (!EventFilter(event)) continue;

      std::cout << "Event " << event.eventID
                << " hits: " << event.globalHits.size() << std::endl;
   }
}

bool Pipeline::EventFilter(Event& event) {
   std::vector<double> z_vals;
   std::vector<double> x_vals;
   std::vector<double> y_vals;

   // 1. 检查每个Tracker是否恰好只有一个GlobalHit
   for (auto& [detID, det] : m_dets) {
      if (!det->isTracker()) continue;

      auto it = event.globalHits.find(detID);
      if (it == event.globalHits.end() || it->second.size() != 1) {
         return false;
      }

      const GlobalHit& hit = it->second[0];
      z_vals.push_back(hit.z);
      x_vals.push_back(hit.x);
      y_vals.push_back(hit.y);
   }

   if (z_vals.size() < 2) {
      return false;
   }

   // 计算z、x、y的均值
   auto mean = [](const std::vector<double>& v) {
      double sum = 0.0;
      for (auto& val : v) sum += val;
      return sum / v.size();
   };

   double z_mean = mean(z_vals);
   double x_mean = mean(x_vals);
   double y_mean = mean(y_vals);

   // 计算斜率和截距
   auto calcSlopeIntercept = [&](const std::vector<double>& z, const std::vector<double>& val, double mean_z, double mean_val) {
      double num = 0.0;
      double den = 0.0;
      for (size_t i = 0; i < z.size(); ++i) {
         double dz = z[i] - mean_z;
         num += dz * (val[i] - mean_val);
         den += dz * dz;
      }
      double slope = num / den;
      double intercept = mean_val - slope * mean_z;
      return std::make_pair(slope, intercept);
   };

   auto [slope_x, intercept_x] = calcSlopeIntercept(z_vals, x_vals, z_mean, x_mean);
   auto [slope_y, intercept_y] = calcSlopeIntercept(z_vals, y_vals, z_mean, y_mean);

   // 计算chi2和ndf
   double chi2 = 0.0;
   for (size_t i = 0; i < z_vals.size(); ++i) {
      double dx = x_vals[i] - (intercept_x + slope_x * z_vals[i]);
      double dy = y_vals[i] - (intercept_y + slope_y * z_vals[i]);
      chi2 += dx * dx + dy * dy;
   }
   int ndf = static_cast<int>(2 * z_vals.size() - 4);

   // 5. 保存结果到event.track
   event.track.slope_x = slope_x;
   event.track.intercept_x = intercept_x;
   event.track.slope_y = slope_y;
   event.track.intercept_y = intercept_y;
   event.track.chi2 = chi2 / ndf;

   return true;
}