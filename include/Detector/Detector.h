#pragma once

#include "Algorithm.h"
#include "DataModel.h"

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

    // Coordinate Transform
    GlobalHit LocalToGlobal(const LocalHit& aLocalHit) const;
    LocalHit GlobalToLocal(const GlobalHit& aGlobalHit) const;

    // Specific Detector Geometry
    virtual GlobalHit GetHitFromTrack(const Track& track) const = 0;
    virtual LocalHit GetLocalHitFromCluster(const RecCluster& cluster) const = 0;

    // Algorithms: rawData->recHits->localHits
    void Reconstruct();

    // Data Interfaces
    void SetRawData(const std::vector<RawData>& rawData) { m_rawData = rawData; }
    void AddRawData(const RawData& raw) { m_rawData.push_back(raw); }
    void AddRecCluster(const RecCluster& cluster) { m_recClusters.push_back(cluster); }
    const std::vector<RawData>& GetRawData() const { return m_rawData; }
    const std::vector<RecCluster>& GetRecClusters() const { return m_recClusters; }
    const std::vector<LocalHit>& GetLocalHits() const { return m_localHits; }
    const int GetNumOfHits() const { return m_localHits.size(); }
    void ClearData() {
        m_rawData.clear();
        m_recClusters.clear();
        m_localHits.clear();
    }

   protected:
    TVector3 m_pos;
    TVector3 m_rot;
    TVector3 m_alignPos = TVector3(0, 0, 0);
    TVector3 m_alignRot = TVector3(0, 0, 0);

    std::vector<RawData> m_rawData;
    std::vector<RecCluster> m_recClusters;
    std::vector<LocalHit> m_localHits;

    std::shared_ptr<Algorithm> m_algorithm;

   protected:
    int m_id;
    std::string m_name;
    Role m_role;
};
