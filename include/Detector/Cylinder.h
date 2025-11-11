#pragma once

#include "Config.h"
#include "DataModel.h"
#include "Detector.h"

class Cylinder : public Detector {
   public:
    Cylinder(int id, const std::string& name, const json& config);

    GlobalHit GetHitFromTrack(const Track& track) const override;

   protected:
    cylinderConfig m_config;
};