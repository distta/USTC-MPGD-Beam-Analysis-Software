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
    bool isDUT() const { return m_role == Role::DUT; }
    bool isTracker() const { return m_role == Role::Tracker; }

    // ---------- 对齐接口 ----------
    void Alignment(double dx, double dy, double dz, double dRotX = 0, double dRotY = 0, double dRotZ = 0);

    virtual GlobalHit LocalToGlobal(const RecHit& recHit) const = 0;

    // ---------- 核心接口 ----------
    std::vector<RecHit> Reconstruction(const std::vector<RawData>& raws) {
        std::vector<RecHit> hits = m_clusterBuilder->Reconstruction(raws);
        for (auto& hit : hits) hit.detID = m_id;
        return hits;
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
