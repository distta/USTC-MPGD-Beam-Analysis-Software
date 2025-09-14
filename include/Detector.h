#pragma once

#include "Algorithms.h"
#include "DataModel.h"
#include <nlohmann/json.hpp>

enum class Role {
   Tracker,
   DUT
};

class Detector {
  public:
   Detector(int id, const std::string& name)
       : m_id(id), m_name(name) {}

   int GetID() const { return m_id; }
   std::string GetName() const { return m_name; }
   bool isDUT() const { return m_role == Role::DUT; }
   bool isTracker() const { return m_role == Role::Tracker; }

   // ---------- 初始化配置 ----------
   void initialize(const nlohmann::json config);

   // ---------- 几何接口 ----------
   LocalHit GlobalToLocal(const GlobalHit& globalHit) const;
   GlobalHit LocalToGlobal(const LocalHit& localHit) const;
   void UpdateGeometry(double dx, double dy, double dz,
                       double dRotX = 0, double dRotY = 0, double dRotZ = 0);

   // ---------- 算法接口 ----------
   std::shared_ptr<IStripHitBuilder> stripHitBuilder;
   std::shared_ptr<IClusterBuilder> clusterBuilder;

   // ---------- 核心接口 ----------
   std::vector<StripHit> BuildStripHit(const std::vector<RawData>& raw) {
      return stripHitBuilder->BuildStripHit(raw);
   };
   std::vector<RecCluster> BuildClusters(const std::vector<StripHit>& strips) {
      return clusterBuilder->BuildCluster(strips);
   };
   std::vector<LocalHit> RecLocalHit(std::vector<RecCluster> clusters);

  private:
   int m_id;
   std::string m_name;
   Role m_role;

   // ---------- 几何属性 ----------
   double m_posX = 0, m_posY = 0, m_posZ = 0;
   double m_rotX = 0, m_rotY = 0, m_rotZ = 0;

   std::map<int, double> m_readoutPlaneAngle;
   std::map<int, double> m_readoutPlanePitch;
};
