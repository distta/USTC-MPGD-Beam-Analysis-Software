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

double Cylinder::GetLocalHit(const Track& track, int type) const {

    // track parameters in global: x(z) = sx * z + ix, y(z) = sy * z + iy
    const double sx = track.slope_x;
    const double sy = track.slope_y;
    const double ix = track.intercept_x;
    const double iy = track.intercept_y;

    // inverse rotate (R(-rot_z)) and translate to detector-local XY
    const double cos_r = std::cos(m_rotZ);
    const double sin_r = std::sin(m_rotZ);
    const double R = m_config.radius;

    // intercept vector relative to detector center, then rotate
    const double dx0 = ix - m_posX;
    const double dy0 = iy - m_posY;
    const double ix_loc = cos_r * dx0 + sin_r * dy0;  // R(-θ) * (ix - m_posX)
    const double iy_loc = -sin_r * dx0 + cos_r * dy0;

    // slopes rotated the same way
    const double sx_loc = cos_r * sx + sin_r * sy;
    const double sy_loc = -sin_r * sx + cos_r * sy;

    // Quadratic in z: A z^2 + B z + C = 0
    const double A = sx_loc * sx_loc + sy_loc * sy_loc;
    const double B = 2.0 * (sx_loc * ix_loc + sy_loc * iy_loc);
    const double C = ix_loc * ix_loc + iy_loc * iy_loc - R * R;

    const double eps = 1e-12;
    std::vector<double> sols;

    if (std::abs(A) < eps) {
        // 降为一次方程 B z + C = 0
        if (std::abs(B) > eps) {
            sols.push_back(-C / B);
        }
    } else {
        const double disc = B * B - 4.0 * A * C;
        if (disc < 0.0) return std::numeric_limits<double>::quiet_NaN();  // 无实根
        const double sqrtD = std::sqrt(disc);
        sols.push_back((-B + sqrtD) / (2.0 * A));
        sols.push_back((-B - sqrtD) / (2.0 * A));
    }

    // 选择与探测器中心 z (cz) 最接近的解（一个常用的策略）
    double best_z = sols[0];
    double best_diff = std::abs(best_z - m_posZ);
    for (double zcand : sols) {
        double d = std::abs(zcand - m_posZ);
        if (d < best_diff) {
            best_diff = d;
            best_z = zcand;
        }
    }

    // 计算交点在局域系的坐标
    const double x_hit = sx_loc * best_z + ix_loc;
    const double y_hit = sy_loc * best_z + iy_loc;
    const double z_local = best_z - m_posZ;
    const double r_hit = std::sqrt(x_hit * x_hit + y_hit * y_hit);
    double phi = std::atan2(y_hit, x_hit);
    if (phi < 0.0) phi += 2.0 * std::acos(-1.0);
    const double arc_len = r_hit * phi;  // φ方向取弧长

    switch (type) {
        case 0:
            return arc_len;  // φ方向（弧长）
        case 1:
            return z_local;  // Z（局域）
        case 2:
            return r_hit;  // 径向
        default:
            return std::numeric_limits<double>::quiet_NaN();
    }
}