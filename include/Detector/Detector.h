#pragma once

#include "Clustering.h"
#include "DataModel.h"

class Detector {
   public:
    enum class Role {
        Tracker,
        DUT
    };
    Detector(int id, const std::string& name, const json& config);
    virtual ~Detector() = default;

    int GetID() const { return m_id; }
    std::string GetName() const { return m_name; }
    double GetPosX() const { return m_posX; }
    double GetPosY() const { return m_posY; }
    double GetPosZ() const { return m_posZ; }
    double GetRotX() const { return m_rotX; }
    double GetRotY() const { return m_rotY; }
    double GetRotZ() const { return m_rotZ; }

    bool isDUT() const { return m_role == Role::DUT; }
    bool isTracker() const { return m_role == Role::Tracker; }

    // ---------- 对齐接口 ----------
    void Alignment(double dx, double dy, double dz, double dRotX, double dRotY, double dRotZ);
    void Alignment(double dx, double dy, double dRotZ);

    virtual GlobalHit LocalToGlobal(const RecHit& recHit) const {};
    virtual double GetLocalHit(const Track& track, int type) const {};

    // ---------- 核心接口 ----------
    std::vector<RecHit> Reconstruction(const std::vector<RawData>& raws) {
        return m_clusterBuilder->Reconstruction(raws);
    };

   protected:
    double m_posX = 0, m_posY = 0, m_posZ = 0;
    double m_rotX = 0, m_rotY = 0, m_rotZ = 0;

   private:
    int m_id;
    std::string m_name;
    Role m_role;

    std::shared_ptr<Clustering> m_clusterBuilder;
};
