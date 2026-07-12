#pragma once

#include "Detector.h"

class Planar : public Detector {
   public:
    Planar(int id, const std::string& name, const json& config);

    GlobalHit CalcHitFromTrack(const Track& track) const override;

    TVector3 GetPlaneNormal() const;

    std::vector<LocalHit> CalcLocalHitsFromClusters(const std::vector<Cluster>& clusters) const override;

    const planarConfig* GetPlanarConfig() const override { return &m_config; }

   private:
    planarConfig m_config;
};
