#include "Detector/PlanarPad.h"

#include <cmath>
#include <stdexcept>

namespace {

planarPadConfig ParsePadConfig(const json& detectorConfig) {
    if (!detectorConfig.contains("pads") ||
        !detectorConfig["pads"].is_object()) {
        throw std::runtime_error(
            "[PlanarPad] Missing or invalid 'pads' configuration");
    }

    const auto& pads = detectorConfig["pads"];
    planarPadConfig config;
    config.rows = pads.value("rows", 0);
    config.columns = pads.value("columns", 0);
    config.pitchX = pads.value("pitchX", 0.0);
    config.pitchY = pads.value("pitchY", 0.0);
    config.sizeX = pads.value("sizeX", config.pitchX);
    config.sizeY = pads.value("sizeY", config.pitchY);
    config.planeType = pads.value("planeType", 0);

    if (config.rows <= 0 || config.columns <= 0) {
        throw std::runtime_error(
            "[PlanarPad] 'rows' and 'columns' must be positive");
    }
    if (config.pitchX <= 0.0 || config.pitchY <= 0.0) {
        throw std::runtime_error(
            "[PlanarPad] 'pitchX' and 'pitchY' must be positive");
    }
    if (config.sizeX <= 0.0 || config.sizeY <= 0.0 ||
        config.sizeX > config.pitchX || config.sizeY > config.pitchY) {
        throw std::runtime_error(
            "[PlanarPad] pad size must be positive and not exceed its pitch");
    }
    return config;
}

}  // namespace

PlanarPad::PlanarPad(int id, const std::string& name, const json& config)
    : Detector(id, name, config), m_padConfig(ParsePadConfig(config)) {}

GlobalHit PlanarPad::CalcHitFromTrack(const Track& track) const {
    const TVector3 rotation = GetRot();
    const double cx = std::cos(rotation.X());
    const double sx = std::sin(rotation.X());
    const double cy = std::cos(rotation.Y());
    const double sy = std::sin(rotation.Y());
    const double cz = std::cos(rotation.Z());
    const double sz = std::sin(rotation.Z());

    TVector3 normal(cz * sy * cx + sz * sx,
                    sz * sy * cx - cz * sx,
                    cy * cx);
    normal = normal.Unit();

    const TVector3 origin(track.bx, track.by, 0.0);
    const TVector3 direction(track.kx, track.ky, 1.0);
    const double denominator = direction.Dot(normal);
    if (std::abs(denominator) < 1e-12) return {};

    const double parameter = (GetPos() - origin).Dot(normal) / denominator;
    if (parameter < 0.0) return {};

    return origin + parameter * direction;
}

TVector3 PlanarPad::PadPosition(int row, int column) const {
    if (row < 0 || row >= m_padConfig.rows ||
        column < 0 || column >= m_padConfig.columns) {
        throw std::out_of_range("[PlanarPad] pad row or column is out of range");
    }
    return {column * m_padConfig.pitchX,
            row * m_padConfig.pitchY, 0.0};
}

bool PlanarPad::ContainsLocal(const TVector3& localPosition) const {
    const double x = localPosition.X();
    const double y = localPosition.Y();
    const double xMinimum = -0.5 * m_padConfig.pitchX;
    const double xMaximum =
        (m_padConfig.columns - 0.5) * m_padConfig.pitchX;
    const double yMinimum = -0.5 * m_padConfig.pitchY;
    const double yMaximum =
        (m_padConfig.rows - 0.5) * m_padConfig.pitchY;
    if (x < xMinimum || x > xMaximum ||
        y < yMinimum || y > yMaximum) {
        return false;
    }

    const double withinPadX = std::fmod(
        x - xMinimum, m_padConfig.pitchX);
    const double withinPadY = std::fmod(
        y - yMinimum, m_padConfig.pitchY);
    const double marginX = 0.5 * (m_padConfig.pitchX - m_padConfig.sizeX);
    const double marginY = 0.5 * (m_padConfig.pitchY - m_padConfig.sizeY);
    return withinPadX >= marginX &&
           withinPadX <= m_padConfig.pitchX - marginX &&
           withinPadY >= marginY &&
           withinPadY <= m_padConfig.pitchY - marginY;
}

std::vector<LocalHit> PlanarPad::CalcLocalHitsFromClusters(
    const std::vector<Cluster>& clusters) const {
    std::vector<LocalHit> localHits;
    localHits.reserve(clusters.size());

    for (size_t index = 0; index < clusters.size(); ++index) {
        const auto& cluster = clusters[index];
        if (cluster.type != m_padConfig.planeType) continue;

        if (!cluster.hasLocalPosition) continue;

        LocalHit hit;
        hit.localPos = cluster.localPosition;
        hit.clusterIndices.push_back(static_cast<int>(index));
        localHits.push_back(std::move(hit));
    }

    return localHits;
}
