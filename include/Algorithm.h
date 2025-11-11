#pragma once

#include "DataModel.h"

class Algorithm {
   public:
    Algorithm() = default;
    explicit Algorithm(const json& config);

    std::vector<RecCluster> Reconstruct(const std::vector<RawData>& raw);

    bool processCluster(Cluster& aCluster);

   private:
    WaveformConfig m_waveformConfig;
    ClusterConfig m_clusterConfig;
    ReconstructionConfig m_reconstructionConfig;

    StripHit processWaveform(const RawData& rawData);

    StripHit processFastWaveform(const RawData& rawData);

    std::vector<Cluster> BuildClusters(const std::vector<StripHit>& stripHits);

    std::vector<RecCluster> MatchClusters(std::map<int, std::vector<Cluster>>& clusters);

    // assisting functions

    void reconstructChargeWeighted(Cluster& cluster);

    void reconstructUTPC(Cluster& cluster);
};