#include "Detector/Detector.h"
#include "TMath.h"
#include "TMatrixD.h"

Detector::Detector(int id, const std::string& name, const nlohmann::json& config) : m_id(id), m_name(name) {
    // 设置探测器角色
    std::string roleStr = config.value("role", "Tracker");
    if (roleStr == "Tracker") {
        m_role = Role::Tracker;
    } else if (roleStr == "DUT") {
        m_role = Role::DUT;
    } else {
        m_role = Role::Ignored;
    }

    if (!config.contains("position") || !config["position"].is_array() || config["position"].size() < 3) {
        throw std::runtime_error("Detector: Missing or invalid 'position' parameter in geometry config");
    }
    if (!config.contains("rotation") || !config["rotation"].is_array() || config["rotation"].size() < 3) {
        throw std::runtime_error("Detector: Missing or invalid 'rotation' parameter in geometry config");
    }

    auto& pos = config["position"];
    m_pos.SetXYZ(pos[0], pos[1], pos[2]);

    auto& rot = config["rotation"];
    m_rot.SetXYZ(rot[0], rot[1], rot[2]);

    // 从配置中选择算法
    if (config.contains("Algorithm")) {
        auto algoConfig = config["Algorithm"];
        m_algorithm = std::make_shared<Algorithm>(algoConfig);
    } else {
        m_algorithm = std::make_shared<Algorithm>();
    }
}

TMatrixD RotationMatrixXYZ(const TVector3& rot) {
    double cx = cos(rot.X()), sx = sin(rot.X());
    double cy = cos(rot.Y()), sy = sin(rot.Y());
    double cz = cos(rot.Z()), sz = sin(rot.Z());

    TMatrixD Rx(3, 3), Ry(3, 3), Rz(3, 3);
    Rx.UnitMatrix();
    Ry.UnitMatrix();
    Rz.UnitMatrix();

    Rx(1, 1) = cx;
    Rx(1, 2) = -sx;
    Rx(2, 1) = sx;
    Rx(2, 2) = cx;

    Ry(0, 0) = cy;
    Ry(0, 2) = sy;
    Ry(2, 0) = -sy;
    Ry(2, 2) = cy;

    Rz(0, 0) = cz;
    Rz(0, 1) = -sz;
    Rz(1, 0) = sz;
    Rz(1, 1) = cz;

    TMatrixD R = Rz * Ry * Rx;  // ZYX rotation
    return R;
}

GlobalHit Detector::LocalToGlobal(const LocalHit& aLocalHit) const {
    // Combined rotation & position
    TVector3 totalRot = GetRot();
    TVector3 totalPos = GetPos();

    TMatrixD R = RotationMatrixXYZ(totalRot);
    TVector3 local(aLocalHit.x(), aLocalHit.y(), aLocalHit.z());
    TVector3 global = totalPos + TVector3(R(0, 0) * local.X() + R(0, 1) * local.Y() + R(0, 2) * local.Z(),
                                          R(1, 0) * local.X() + R(1, 1) * local.Y() + R(1, 2) * local.Z(),
                                          R(2, 0) * local.X() + R(2, 1) * local.Y() + R(2, 2) * local.Z());

    return global;
}

LocalHit Detector::GlobalToLocal(const GlobalHit& aGlobalHit) const {
    TVector3 totalRot = GetRot();
    TVector3 totalPos = GetPos();

    TMatrixD R = RotationMatrixXYZ(totalRot);
    TMatrixD Rinv = TMatrixD(TMatrixD::kTransposed, R);

    TVector3 global(aGlobalHit.x(), aGlobalHit.y(), aGlobalHit.z());
    TVector3 localVec = global - totalPos;
    TVector3 local(Rinv(0, 0) * localVec.X() + Rinv(0, 1) * localVec.Y() + Rinv(0, 2) * localVec.Z(),
                   Rinv(1, 0) * localVec.X() + Rinv(1, 1) * localVec.Y() + Rinv(1, 2) * localVec.Z(),
                   Rinv(2, 0) * localVec.X() + Rinv(2, 1) * localVec.Y() + Rinv(2, 2) * localVec.Z());

    return local;
}

void Detector::Reconstruct() {

    if (!m_algorithm) {
        throw std::runtime_error("No clustering algorithm set for detector " + m_name);
    }

    if (m_rawData.empty()) {
        return;
    }

    m_recClusters = m_algorithm->Reconstruct(m_rawData);

    for (auto& RecCluster : m_recClusters) {
        m_localHits.push_back(GetLocalHitFromCluster(RecCluster));
    }
}
