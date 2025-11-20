#pragma once

#include "Detector.h"

class Planar : public Detector {
   public:
    Planar(int id, const std::string& name, const json& config);

    GlobalHit GetHitFromTrack(const Track& track) const override;

    TVector3 GetPlaneNormal() const;

    LocalHit GetLocalHitFromCluster(const RecCluster& cluster) const override;
};