#pragma once

#include "Detector.h"

class PlanarPad : public Detector {
   public:
    PlanarPad(int id, const std::string& name, const json& config);

    GlobalHit CalcHitFromTrack(const Track& track) const override;
    std::vector<LocalHit> CalcLocalHitsFromClusters(
        const std::vector<Cluster>& clusters) const override;

    const planarPadConfig& GetPadConfig() const { return m_padConfig; }
    const planarPadConfig* GetPlanarPadConfig() const override {
        return &m_padConfig;
    }

    TVector3 PadPosition(int row, int column) const;
    bool ContainsLocal(const TVector3& localPosition) const;

   private:
    planarPadConfig m_padConfig;
};
