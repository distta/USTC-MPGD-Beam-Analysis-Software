#include "DetectorFrame.h"

bool DetectorFrame::AnalyzeRaw() {
    if (m_raw.empty()) return false;

    auto hitProcessor = m_det.GetAlgorithm<HitProcessor>("HitProcessor");
    return hitProcessor->Process(*this);
}

// Clustering ChannelHit
bool DetectorFrame::Clustering() {
    if (m_channelHits.empty()) return false;

    auto clusterBuilder = m_det.GetAlgorithm<ClusterBuilder>("ClusterBuilder");
    return clusterBuilder->Process(*this);
}

// Reconstruction: 重建 Cluster 位置并生成 LocalHit
bool DetectorFrame::Reconstruct() {
    if (m_clusters.empty()) return false;

    // Step 1: 调用 ClusterReconstructor 计算每个 Cluster 的 pos
    auto clusterRecon = m_det.GetAlgorithm<ClusterReconstructor>("ClusterReconstructor");
    clusterRecon->Process(*this);

    // Step 2: 调用Detector几何逻辑生成LocalHits
    m_localHits = m_det.CalcLocalHitsFromClusters(m_clusters);

    return !m_localHits.empty();
}

bool DetectorFrame::TransformToGlobal() {
    if (m_localHits.empty()) return false;

    // 预分配内存优化
    m_globalHits.reserve(m_localHits.size());

    // LocalHits -> GlobalHits
    for (const auto& lh : m_localHits) {
        m_globalHits.push_back(m_det.LocalToGlobal(lh.localPos));
    }

    return true;
}

bool DetectorFrame::Process(double t0) {
    this->t0 = t0;
    return AnalyzeRaw() && Clustering() && Reconstruct() && TransformToGlobal();
}

bool DetectorFrame::Process() {
    t0 = std::numeric_limits<double>::quiet_NaN();
    return AnalyzeRaw() && Clustering() && Reconstruct() && TransformToGlobal();
}
