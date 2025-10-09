#pragma once

#include "Config.h"
#include "DataModel.h"
#include "Detector.h"

class Cylinder : public Detector {
   public:
    Cylinder(int id, const std::string& name, const json& config);

    double GetLocalHit(const Track& track, int type) const override;

   protected:
    cylinderConfig m_config;
};