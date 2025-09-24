#pragma once

#include "Config.h"
#include "DataModel.h"
#include "Detector.h"

class Planar : public Detector {
   public:
    Planar(int id, const std::string& name, const json& config);

    GlobalHit LocalToGlobal(const RecHit& local) const override;

   protected:
    planarConfig m_config;
};