#include "Detector/Planar.h"
#include "DataModel.h"
#include "iostream"
#include <cmath>
#include <map>

Planar::Planar(int id, const std::string& name, const json& config) : Detector(id, name, config) {

    if (config.contains("readoutPlanes")) {
        m_config.readoutPlaneAngle.clear();
        m_config.readoutPlanePitch.clear();
        m_config.readoutPlaneStripNumber.clear();
        for (auto& [planeIDStr, planeData] : config["readoutPlanes"].items()) {
            try {
                int type = std::stoi(planeIDStr);
                m_config.readoutPlaneAngle.emplace(type, planeData["angle"]);
                m_config.readoutPlanePitch.emplace(type, planeData["pitch"]);
                m_config.readoutPlaneStripNumber.emplace(type, int(planeData["num"]));
            } catch (std::exception& e) {
                std::cout << "Error: " << e.what() << " in readout plane configuration" << std::endl;
            }
        }
    }
}

GlobalHit Planar::GetHitFromTrack(const Track& track) const {
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

LocalHit Planar::GetLocalHitFromCluster(const RecCluster& recCluster) const {

    LocalHit localHit = TVector3(-999, -999, 0);
    for (int i = 0; i < recCluster.size(); i++) {
        if (recCluster[i].type == 0) {
            localHit.SetX(recCluster[i].pos * m_config.readoutPlanePitch.at(0));
        } else if (recCluster[i].type == 1) {
            localHit.SetY(recCluster[i].pos * m_config.readoutPlanePitch.at(1));
        }
    }

    return localHit;
}
