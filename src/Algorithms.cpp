#include "Algorithms.h"
#include "DataModel.h"
#include <nlohmann/json.hpp>

using json = nlohmann::json;

class DefaultStripHitBuilder : public IStripHitBuilder {
  public:
   DefaultStripHitBuilder(const json& params) {
      // 初始化参数，如果有的话
      if (params.contains("threshold")) {
         threshold_ = params["threshold"];
      }
   }

   const std::vector<StripHit> BuildStripHit(const std::vector<RawData>& raw) override {
      std::vector<StripHit> hits;
      // 示例：构建 StripHit（根据需求填充）
      for (const auto& rawData : raw) {
         StripHit hit;
         hit.stripID = rawData.StripID;
         hit.type = rawData.type;
         hit.raw = rawData;
         // 根据具体需求进行波形分析
         hits.push_back(hit);
      }
      return hits;
   }

  private:
   double threshold_ = 0.0;  // 例如：自定义的参数
};

class ChargeWeightedCluster : public IClusterBuilder {
  public:
   ChargeWeightedCluster(const json& params) {
      // 初始化参数
      if (params.contains("maxClusterDistance")) {
         maxClusterDistance_ = params["maxClusterDistance"];
      }
   }

   const std::vector<RecCluster> BuildCluster(const std::vector<StripHit>& hits) override {
      std::vector<RecCluster> clusters;
      RecCluster currentCluster;

      for (const auto& hit : hits) {
         // 示例：按一定逻辑聚类
         if (currentCluster.strips.empty() ||
             std::abs(hit.stripID - currentCluster.strips.back().stripID) <= maxClusterDistance_) {
            currentCluster.strips.push_back(hit);
         } else {
            clusters.push_back(currentCluster);
            currentCluster.strips.clear();
            currentCluster.strips.push_back(hit);
         }
      }

      // 最后一组聚类也需要保存
      if (!currentCluster.strips.empty()) {
         clusters.push_back(currentCluster);
      }
      return clusters;
   }

  private:
   int maxClusterDistance_ = 5;  // 聚类的最大距离
};
