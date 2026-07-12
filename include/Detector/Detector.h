#pragma once

#include "Algorithm/IAlgorithm.h"
#include "Config.h"
#include "Event/DataModel.h"
#include <map>

class Detector {
   public:
    enum class Role {
        Tracker,
        DUT,
        Ignored
    };
    Detector(int id, const std::string& name, const json& config);
    virtual ~Detector() = default;

    bool isDUT() const { return m_role == Role::DUT; }
    bool isTracker() const { return m_role == Role::Tracker; }

    int GetID() const { return m_id; }
    std::string GetName() const { return m_name; }

    TVector3 GetPos() const { return m_pos + m_alignPos; }
    TVector3 GetRot() const { return m_rot + m_alignRot; }
    TVector3 GetAlignPos() const { return m_alignPos; }
    TVector3 GetAlignRot() const { return m_alignRot; }

    void SetPos(const TVector3& pos) { m_pos = pos; }
    void SetRot(const TVector3& rot) { m_rot = rot; }
    void SetAlignment(double dx, double dy, double dz,
                      double dRotX, double dRotY, double dRotZ) {
        m_alignPos.SetXYZ(dx, dy, dz);
        m_alignRot.SetXYZ(dRotX, dRotY, dRotZ);
    }
    virtual const planarConfig* GetPlanarConfig() const { return nullptr; }
    virtual const planarPadConfig* GetPlanarPadConfig() const { return nullptr; }
    virtual const cylinderConfig* GetCylinderConfig() const { return nullptr; }

    // Coordinate Transform
    TVector3 LocalToGlobal(const TVector3& localPos) const;
    TVector3 GlobalToLocal(const TVector3& globalPos) const;

    // Specific Detector Geometry
    virtual TVector3 CalcHitFromTrack(const Track& track) const = 0;
    virtual std::vector<LocalHit> CalcLocalHitsFromClusters(const std::vector<Cluster>& clusters) const = 0;

    // 算法相关接口
    template <typename T>
    std::shared_ptr<T> GetAlgorithm(const std::string& name) const {
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

   protected:
    TVector3 m_pos;
    TVector3 m_rot;
    TVector3 m_alignPos = TVector3(0, 0, 0);
    TVector3 m_alignRot = TVector3(0, 0, 0);

    // 算法实例
    std::map<std::string, std::shared_ptr<IAlgorithm>> m_algorithms;

   protected:
    int m_id;
    std::string m_name;
    Role m_role;
};
