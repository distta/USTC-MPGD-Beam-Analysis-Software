#include "Detector/Planar.h"
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


LocalHit Planar::GlobalToLocal(const GlobalHit& globalHit) const {
    // 假设旋转角度较小，使用简单的平移和旋转转换
    LocalHit local;
    double dx = globalHit.x - m_posX;
    double dy = globalHit.y - m_posY;
    double dz = globalHit.z - m_posZ;

    // 这里可扩展为完整的旋转矩阵
    local.u = dx * cos(m_rotZ) + dy * sin(m_rotZ);

    // 如果是二维探测器
    if (m_config.readoutPlaneAngle.size() > 1) {
        local.v = dx * -sin(m_rotZ) + dy * cos(m_rotZ);
    }

    return local;
}

GlobalHit Planar::LocalToGlobal(const LocalHit& localHit) const {
    GlobalHit global;
    double u = localHit.u;
    double v = localHit.v.value_or(0.0);

    // 简单旋转+平移
    global.x = m_posX + u * cos(m_rotZ) - v * sin(m_rotZ);
    global.y = m_posY + u * sin(m_rotZ) + v * cos(m_rotZ);
    global.z = m_posZ;

    return global;
}
