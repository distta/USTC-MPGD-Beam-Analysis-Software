#pragma once

#include "DataModel.h"
#include <vector>

class IStripHitBuilder {
  public:
   virtual ~IStripHitBuilder() = default;
   virtual const std::vector<StripHit> BuildStripHit(const std::vector<RawData>& raw) = 0;
};

class IClusterBuilder {
  public:
   virtual ~IClusterBuilder() = default;
   virtual const std::vector<RecCluster> BuildCluster(const std::vector<StripHit>& hits) = 0;
};