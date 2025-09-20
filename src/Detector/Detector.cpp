#include "Detector/Detector.h"
#include "Algorithms/ClusterMatcher.h"
#include "Algorithms/Clustering.h"

Detector::Detector(int id, const std::string& name, const nlohmann::json& config) : m_id(id), m_name(name) {
    // 设置探测器角色
    std::string roleStr = config.value("role", "Tracker");
    if (roleStr == "Tracker") {
        m_role = Role::Tracker;
    } else if (roleStr == "DUT") {
        m_role = Role::DUT;
    } else {
        m_role = Role::Tracker;  // 默认为Tracker
    }

    if (!config.contains("position") || !config["position"].is_array() || config["position"].size() < 3) {
        throw std::runtime_error("Detector: Missing or invalid 'position' parameter in geometry config");
    }
    if (!config.contains("rotation") || !config["rotation"].is_array() || config["rotation"].size() < 3) {
        throw std::runtime_error("Detector: Missing or invalid 'rotation' parameter in geometry config");
    }

    auto pos = config["position"];
    m_posX = pos[0].get<double>();
    m_posY = pos[1].get<double>();
    m_posZ = pos[2].get<double>();

    auto rot = config["rotation"];
    m_rotX = rot[0].get<double>();
    m_rotY = rot[1].get<double>();
    m_rotZ = rot[2].get<double>();

    if (config.contains("clusterBuilder")) {
        m_clusterBuilder = std::make_shared<Clustering>(config["clusterBuilder"]);
    } else {
        m_clusterBuilder = std::make_shared<Clustering>(nlohmann::json{});
    }

    // HitCreator
    if (config.contains("hitCreator")) {
        m_clusterMatcher = std::make_shared<ClusterMatcher>(config["clusterMatcher"]);
    } else {
        m_clusterMatcher = std::make_shared<ClusterMatcher>(nlohmann::json{});
    }
}

void Detector::Alignment(double dx, double dy, double dz, double dRotX, double dRotY, double dRotZ) {
    m_posX += dx;
    m_posY += dy;
    m_posZ += dz;
    m_rotX += dRotX;
    m_rotY += dRotY;
    m_rotZ += dRotZ;
}
