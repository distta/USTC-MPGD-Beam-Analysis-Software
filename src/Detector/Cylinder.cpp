#include "Detector/Cylinder.h"
#include "DataModel.h"
#include "iostream"
#include <cmath>
#include <map>

Cylinder::Cylinder(int id, const std::string& name, const json& config) : Detector(id, name, config) {

    if (config.contains("readoutPlanes")) {
        m_config.readoutPlaneAngle.clear();
        m_config.readoutPlanePitch.clear();
        for (auto& [planeIDStr, planeData] : config["readoutPlanes"].items()) {
            try {
                int type = std::stoi(planeIDStr);
                m_config.readoutPlaneAngle.emplace(type, planeData["angle"]);
                m_config.readoutPlanePitch.emplace(type, planeData["pitch"]);
            } catch (std::exception& e) {
                std::cout << "Error: " << e.what() << " in readout plane configuration" << std::endl;
            }
        }
    }

    if (config.contains("radius")) {
        m_config.radius = config["radius"];
    }
}

GlobalHit Cylinder::CalcHitFromTrack(const Track& track) const {
    // track parameters in global: x(z) = kx * z + bx, y(z) = ky * z + by
    const double kx = track.kx;
    const double ky = track.ky;
    const double bx = track.bx;
    const double by = track.by;

    // inverse rotate (R(-rot.z())) and translate to detector-local XY
    const double cos_r = std::cos(m_rot.Z());
    const double sin_r = std::sin(m_rot.Z());
    const double R = m_config.radius;

    // intercept vector relative to detector center, then rotate
    const double dx0 = bx - m_pos.X();
    const double dy0 = by - m_pos.Y();
    const double bx_loc = cos_r * dx0 + sin_r * dy0;  // R(-θ) * (bx - m_pos.X())
    const double by_loc = -sin_r * dx0 + cos_r * dy0;

    // slopes rotated the same way
    const double kx_loc = cos_r * kx + sin_r * ky;
    const double ky_loc = -sin_r * kx + cos_r * ky;

    // Quadratic in z: A z^2 + B z + C = 0
    const double A = kx_loc * kx_loc + ky_loc * ky_loc;
    const double B = 2.0 * (kx_loc * bx_loc + ky_loc * by_loc);
    const double C = bx_loc * bx_loc + by_loc * by_loc - R * R;

    const double eps = 1e-12;
    std::vector<double> sols;

    if (std::abs(A) < eps) {
        // 降为一次方程 B z + C = 0
        if (std::abs(B) > eps) {
            sols.push_back(-C / B);
        }
    } else {
        const double disc = B * B - 4.0 * A * C;
        if (disc < 0.0) {
            return GlobalHit();  // 无实根，返回默认构造的GlobalHit
        }
        const double sqrtD = std::sqrt(disc);
        sols.push_back((-B + sqrtD) / (2.0 * A));
        sols.push_back((-B - sqrtD) / (2.0 * A));
    }

    if (sols.empty()) {
        return GlobalHit();  // 无解，返回默认构造的GlobalHit
    }

    // 选择与探测器中心 z (m_pos.Z()) 最接近的解（一个常用的策略）
    double best_z = sols[0];
    double best_diff = std::abs(best_z - m_pos.Z());
    for (double zcand : sols) {
        double d = std::abs(zcand - m_pos.Z());
        if (d < best_diff) {
            best_diff = d;
            best_z = zcand;
        }
    }

    // 计算交点在局域系的坐标
    const double x_hit = kx_loc * best_z + bx_loc;
    const double y_hit = ky_loc * best_z + by_loc;
    const double z_local = best_z - m_pos.Z();
    const double r_hit = std::sqrt(x_hit * x_hit + y_hit * y_hit);
    double phi = std::atan2(y_hit, x_hit);
    if (phi < 0.0) phi += 2.0 * std::acos(-1.0);
    const double arc_len = r_hit * phi;  // φ方向取弧长

    // 返回全局坐标系中的击中点
    return GlobalHit(m_pos.X() + x_hit, m_pos.Y() + y_hit, m_pos.Z() + z_local);
}

std::vector<LocalHit> Cylinder::CalcLocalHitsFromClusters(const std::vector<Cluster>& clusters) const {
    std::vector<LocalHit> localHits;

    if (clusters.empty()) return localHits;

    // 按type分组clusters
    std::map<int, std::vector<size_t>> clustersByType;
    for (size_t i = 0; i < clusters.size(); ++i) {
        clustersByType[clusters[i].type].push_back(i);
    }

    // 获取所有type
    std::vector<int> types;
    for (const auto& [type, _] : clustersByType) {
        types.push_back(type);
    }

    // 根据type数量决定匹配策略
    if (types.size() == 1) {
        // 单type：每个cluster直接生成LocalHit
        int type = types[0];
        for (size_t idx : clustersByType[type]) {
            LocalHit lh;
            if (type == 0) {
                lh.localPos.SetXYZ(clusters[idx].pos * m_config.readoutPlanePitch.at(type), 0, 0);
            } else if (type == 1) {
                lh.localPos.SetXYZ(0, clusters[idx].pos * m_config.readoutPlanePitch.at(type), 0);
            }
            lh.clusterIndices.push_back(idx);
            localHits.push_back(lh);
        }
    } else if (types.size() == 2) {
        // 双type：组合匹配
        int type0 = types[0];
        int type1 = types[1];

        const auto& indices0 = clustersByType[type0];
        const auto& indices1 = clustersByType[type1];

        // 对于圆柱探测器，使用类似平面的匹配策略
        for (size_t idx0 : indices0) {
            for (size_t idx1 : indices1) {
                LocalHit lh;
                double x = 0, y = 0;

                if (type0 == 0) {
                    x = clusters[idx0].pos * m_config.readoutPlanePitch.at(type0);
                    y = clusters[idx1].pos * m_config.readoutPlanePitch.at(type1);
                } else {
                    x = clusters[idx1].pos * m_config.readoutPlanePitch.at(type1);
                    y = clusters[idx0].pos * m_config.readoutPlanePitch.at(type0);
                }

                lh.localPos.SetXYZ(x, y, 0);
                lh.clusterIndices.push_back(idx0);
                lh.clusterIndices.push_back(idx1);
                localHits.push_back(lh);
            }
        }
    } else {
        // 多type：复杂匹配逻辑（暂不实现）
        std::cerr << "[Cylinder] Warning: More than 2 types not supported yet" << std::endl;
    }

    return localHits;
}
