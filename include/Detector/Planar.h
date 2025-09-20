#pragma once

#include "Config.h"
#include "DataModel.h"
#include "Detector.h"

class Planar : public Detector {
   public:
    Planar(int id, const std::string& name, const json& config);

    LocalHit GlobalToLocal(const GlobalHit& global) const override;
    GlobalHit LocalToGlobal(const LocalHit& local) const override;

   protected:
    planarConfig m_config;
};