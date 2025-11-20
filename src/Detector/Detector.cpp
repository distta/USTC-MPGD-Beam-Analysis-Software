#include "Detector/Detector.h"
#include "Algorithm/AlgorithmFactory.h"
#include "Algorithm/algorithms/ClusterBuilder.h"
#include "Algorithm/algorithms/ClusterReconstructor.h"
#include "Algorithm/algorithms/WaveformProcessor.h"
#include "TMath.h"
#include "TMatrixD.h"
#include <iostream>
#include <stdexcept>

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
        throw std::runtime_error("[Detector] Missing or invalid 'position' parameter in geometry config");
    }
    if (!config.contains("rotation") || !config["rotation"].is_array() || config["rotation"].size() < 3) {
        throw std::runtime_error("[Detector] Missing or invalid 'rotation' parameter in geometry config");
    }

    auto& pos = config["position"];
    m_pos.SetXYZ(pos[0], pos[1], pos[2]);

    auto& rot = config["rotation"];
    m_rot.SetXYZ(rot[0], rot[1], rot[2]);

    // 算法加载 - 支持Algorithms（多算法）配置
    auto& factory = AlgorithmFactory::Instance();

    if (config.contains("Algorithms")) {

        for (auto& [algoName, cfg] : config["Algorithms"].items()) {

            try {
                auto algorithm = factory.CreateAlgorithm(algoName, cfg);
                algorithm->SetDetector(this);  // 设置算法所属的探测器
                m_algorithms[algoName] = algorithm;
                std::cout << "[Detector " << m_name << "] Loaded algorithm '" << algoName << std::endl;
                algorithm->Print();
            } catch (const std::exception& e) {
                std::cerr << "[Detector " << m_name << "] Failed to create algorithm '" << algoName
                          << "': " << e.what() << std::endl;
                throw;
            }
        }

    } else {
        try {
            // 创建默认的waveformProcessor
            auto waveformProc = factory.CreateAlgorithm("WaveformProcessor", json::object());
            waveformProc->SetDetector(this);
            m_algorithms["WaveformProcessor"] = waveformProc;

            // 创建默认的clustering
            auto clustering = factory.CreateAlgorithm("ClusterBuilder", json::object());
            clustering->SetDetector(this);
            m_algorithms["ClusterBuilder"] = clustering;

            // 创建默认的reconstruction
            auto reconstruction = factory.CreateAlgorithm("ClusterReconstructor", json::object());
            reconstruction->SetDetector(this);
            m_algorithms["ClusterReconstructor"] = reconstruction;
        } catch (const std::exception& e) {
            std::cerr << "[Detector " << m_name << "] Failed to create default algorithms: " << e.what() << std::endl;
            throw;
        }
    }
}

TMatrixD RotationMatrixZYX(const TVector3& rot) {
    double cx = cos(rot.X()), sx = sin(rot.X());
    double cy = cos(rot.Y()), sy = sin(rot.Y());
    double cz = cos(rot.Z()), sz = sin(rot.Z());

    TMatrixD R(3, 3);

    // R = Rz * Ry * Rx
    R(0, 0) = cz * cy;
    R(0, 1) = cz * sy * sx - sz * cx;
    R(0, 2) = cz * sy * cx + sz * sx;

    R(1, 0) = sz * cy;
    R(1, 1) = sz * sy * sx + cz * cx;
    R(1, 2) = sz * sy * cx - cz * sx;

    R(2, 0) = -sy;
    R(2, 1) = cy * sx;
    R(2, 2) = cy * cx;

    return R;
}

GlobalHit Detector::LocalToGlobal(const LocalHit& lh) const {
    TMatrixD R = RotationMatrixZYX(GetRot());
    TVector3 p = TVector3(lh.x(), lh.y(), lh.z());

    TVector3 global(
        R(0, 0) * p.X() + R(0, 1) * p.Y() + R(0, 2) * p.Z(),
        R(1, 0) * p.X() + R(1, 1) * p.Y() + R(1, 2) * p.Z(),
        R(2, 0) * p.X() + R(2, 1) * p.Y() + R(2, 2) * p.Z());

    return global + GetPos();
}

LocalHit Detector::GlobalToLocal(const GlobalHit& gh) const {
    TMatrixD R = RotationMatrixZYX(GetRot());
    TMatrixD Rinv(TMatrixD::kTransposed, R);

    TVector3 p = TVector3(gh.x(), gh.y(), gh.z()) - GetPos();

    TVector3 local(
        Rinv(0, 0) * p.X() + Rinv(0, 1) * p.Y() + Rinv(0, 2) * p.Z(),
        Rinv(1, 0) * p.X() + Rinv(1, 1) * p.Y() + Rinv(1, 2) * p.Z(),
        Rinv(2, 0) * p.X() + Rinv(2, 1) * p.Y() + Rinv(2, 2) * p.Z());

    return local;
}

void Detector::Reconstruct() {
    if (m_rawData.empty()) {
        return;
    }

    // 阶段1：波形处理 RawData -> StripHit
    auto waveformProc = GetAlgorithm<WaveformProcessor>("WaveformProcessor");

    for (const auto& raw : m_rawData) {
        m_StripHits[raw.type].push_back(waveformProc->ProcessWaveform(raw));
    }

    // 阶段2：聚类 StripHit -> Cluster
    auto clusterBuilder = GetAlgorithm<ClusterBuilder>("ClusterBuilder");
    m_recClusters = clusterBuilder->BuildClusters(m_StripHits);

    // 阶段3：位置重建 修改Cluster的pos
    auto clusterRecon = GetAlgorithm<ClusterReconstructor>("ClusterReconstructor");
    for (auto& RecCluster : m_recClusters) {
        for (auto& cluster : RecCluster) {
            clusterRecon->ReconstructPosition(cluster);
        }
    }

    // 转换RecCluster为LocalHit
    for (auto& RecCluster : m_recClusters) {
        m_localHits.push_back(GetLocalHitFromCluster(RecCluster));
    }
}

template <typename T>
std::shared_ptr<T> Detector::GetAlgorithm(const std::string& name) const {

    auto it = m_algorithms.find(name);
    if (it == m_algorithms.end()) {
        throw std::runtime_error("Algorithm '" + name + "' not found in detector " + m_name);
    }
    auto algo = it->second;

    auto casted = std::dynamic_pointer_cast<T>(algo);
    if (!casted) {
        throw std::runtime_error("Algorithm type mismatch for: " + name);
    }
    return casted;
}
