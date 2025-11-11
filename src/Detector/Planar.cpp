#include "Detector/Planar.h"
#include "DataModel.h"
#include "iostream"
#include <cmath>
#include <map>

Planar::Planar(int id, const std::string& name, const json& config) : Detector(id, name, config) {

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
}

GlobalHit Planar::GetHitFromTrack(const Track& track) const {
    GlobalHit globalHit;
    double rx = m_rot.X(), ry = m_rot.Y(), rz = m_rot.Z();

    double cosRx = cos(rx), sinRx = sin(rx);
    double cosRy = cos(ry), sinRy = sin(ry);
    double cosRz = cos(rz), sinRz = sin(rz);

    double nx = sinRy;
    double ny = -sinRx * cosRy;
    double nz = cosRx * cosRy;

    TVector3 normal(nx, ny, nz);
    normal = normal.Unit();

    TVector3 trackOrigin(track.bx, track.by, 0);
    TVector3 trackDirection(track.kx, track.ky, 1);
    trackDirection = trackDirection.Unit();

    TVector3 planePoint = m_pos;
    TVector3 diff = planePoint - trackOrigin;
    double denominator = trackDirection.Dot(normal);

    if (fabs(denominator) < 1e-12) {
        return globalHit;
    }

    double t = diff.Dot(normal) / denominator;

    if (t < 0) {
        return globalHit;
    }

    globalHit = trackOrigin + t * trackDirection;

    return globalHit;
}

LocalHit Planar::GetLocalHitFromCluster(const RecCluster& recCluster) const {
    double pitchX = m_config.readoutPlanePitch.at(0);
    double pitchY = m_config.readoutPlanePitch.at(1);

    LocalHit localHit = TVector3(0, 0, 0);
    for (int i = 0; i < recCluster.size(); i++) {
        if (recCluster[i].type == 0) {
            localHit.SetX(recCluster[i].pos * pitchX);
        } else if (recCluster[i].type == 1) {
            localHit.SetY(recCluster[i].pos * pitchY);
        }
    }

    return localHit;
}
