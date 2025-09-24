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


GlobalHit Planar::LocalToGlobal(const RecHit& recHit) const {
    GlobalHit global;
    double u = recHit.cluster[0].pos;
    double v = recHit.cluster[1].pos;

    // 简单旋转+平移
    global.x = m_posX + u * cos(m_rotZ) - v * sin(m_rotZ);
    global.y = m_posY + u * sin(m_rotZ) + v * cos(m_rotZ);
    global.z = m_posZ;

    return global;
}
