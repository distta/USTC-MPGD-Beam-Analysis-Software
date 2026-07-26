#include "Detector/Planar.h"
#include "DataModel.h"
#include "TMath.h"
#include "iostream"
#include <algorithm>
#include <cmath>
#include <map>

Planar::Planar(int id, const std::string& name, const json& config) : Detector(id, name, config) {

    if (config.contains("readoutPlanes")) {
        m_config.readoutPlaneType.clear();
        m_config.readoutPlaneAngle.clear();
        m_config.readoutPlanePitch.clear();
        m_config.readoutPlaneStripNumber.clear();
        for (auto& [planeIDStr, planeData] : config["readoutPlanes"].items()) {
            try {
                int type = std::stoi(planeIDStr);
                m_config.readoutPlaneType.push_back(type);
                m_config.readoutPlaneAngle.emplace(type, planeData["angle"]);
                m_config.readoutPlanePitch.emplace(type, planeData["pitch"]);
                m_config.readoutPlaneStripNumber.emplace(type, int(planeData["num"]));
            } catch (std::exception& e) {
                std::cout << "Error: " << e.what() << " in readout plane configuration" << std::endl;
            }
        }
        std::sort(m_config.readoutPlaneType.begin(),
                  m_config.readoutPlaneType.end());
    }
}

GlobalHit Planar::CalcHitFromTrack(const Track& track) const {
    GlobalHit globalHit;

    TVector3 rot = GetRot();
    double rx = rot.X(), ry = rot.Y(), rz = rot.Z();

    double cx = cos(rx), sx = sin(rx);
    double cy = cos(ry), sy = sin(ry);
    double cz = cos(rz), sz = sin(rz);

    // 正确的法向量计算（旋转矩阵的第三列）
    TVector3 normal(cz * sy * cx + sz * sx,
                    sz * sy * cx - cz * sx,
                    cy * cx);
    normal = normal.Unit();

    // 更合理的轨迹参数化
    TVector3 trackOrigin(track.bx, track.by, 0);

    // 假设track包含完整的方向向量
    TVector3 trackDirection(track.kx, track.ky, 1);
    trackDirection = trackDirection.Unit();

    TVector3 planePoint = GetPos();
    TVector3 diff = planePoint - trackOrigin;
    double denominator = trackDirection.Dot(normal);

    if (fabs(denominator) < 1e-12) {
        return globalHit;  // 平行，无交点
    }

    double t = diff.Dot(normal) / denominator;
    if (t < 0) {
        return globalHit;  // 反向，无交点
    }

    globalHit = trackOrigin + t * trackDirection;
    return globalHit;
}

std::vector<LocalHit> Planar::CalcLocalHitsFromClusters(const std::vector<Cluster>& clusters) const {
    std::vector<LocalHit> localHits;

    if (clusters.empty()) return localHits;

    // 按type分组clusters
    std::map<int, std::vector<size_t>> clustersByType;  // type -> cluster indices
    for (size_t i = 0; i < clusters.size(); ++i) {
        clustersByType[clusters[i].type].push_back(i);
    }

    // 获取所有type
    std::vector<int> types = m_config.readoutPlaneType;
    const auto angle = [&](int type) {
        return m_config.readoutPlaneAngle.at(type) * TMath::DegToRad();
    };
    const auto coordinate = [&](const Cluster& cluster) {
        return cluster.pos * m_config.readoutPlanePitch.at(cluster.type);
    };

    // 根据type数量决定匹配策略
    if (types.size() == 1) {
        // 单type：每个cluster直接生成LocalHit
        int type = types[0];
        for (size_t idx : clustersByType[type]) {
            LocalHit lh;
            const double phi = angle(type);
            const double value = coordinate(clusters[idx]);
            lh.localPos.SetXYZ(value * std::cos(phi),
                               value * std::sin(phi), 0);
            lh.clusterIndices.push_back(idx);
            localHits.push_back(lh);
        }
    } else if (types.size() == 2) {
        // 双type：X-Y匹配
        int type0 = types[0];
        int type1 = types[1];

        const auto& indices0 = clustersByType[type0];
        const auto& indices1 = clustersByType[type1];

        // 简单匹配：所有组合（可根据时间、电荷等优化）
        for (size_t idx0 : indices0) {
            for (size_t idx1 : indices1) {
                LocalHit lh;
                const double phi0 = angle(type0);
                const double phi1 = angle(type1);
                const double c0 = std::cos(phi0);
                const double s0 = std::sin(phi0);
                const double c1 = std::cos(phi1);
                const double s1 = std::sin(phi1);
                const double determinant = c0 * s1 - s0 * c1;
                if (std::abs(determinant) < 1e-12) continue;
                const double u0 = coordinate(clusters[idx0]);
                const double u1 = coordinate(clusters[idx1]);
                const double x = (u0 * s1 - s0 * u1) / determinant;
                const double y = (c0 * u1 - u0 * c1) / determinant;

                lh.localPos.SetXYZ(x, y, 0);
                lh.clusterIndices.push_back(idx0);
                lh.clusterIndices.push_back(idx1);
                localHits.push_back(lh);
            }
        }
    } else {
        // 多type：复杂匹配逻辑（暂不实现）
        std::cerr << "[Planar] Warning: More than 2 types not supported yet" << std::endl;
    }

    return localHits;
}
