#include "Detector.h"
#include <iostream>

void Detector::initialize(const nlohmann::json config) {

   if (config.contains("geometry")) {
      auto pos = config["geometry"]["position"];
      m_posX = pos[0];
      m_posY = pos[1];
      m_posZ = pos[2];

      auto rot = config["geometry"]["rotation"];
      m_rotX = rot[0];
      m_rotY = rot[1];
      m_rotZ = rot[2];

      // 读出平面
      for (auto& [planeIDStr, planeData] : config["geometry"]["readoutPlanes"].items()) {
         int planeID = std::stoi(planeIDStr);
         m_readoutPlaneAngle[planeID] = planeData["angle"];
         m_readoutPlanePitch[planeID] = planeData["pitch"];
      }
   }

   if (config.contains("algorithms")) {
      auto algoConfig = config["algorithms"];

      // StripHitBuilder
      std::string stripAlgo = algoConfig["stripHitBuilder"]["name"];
      if (stripAlgo == "DefaultStripHitBuilder") {
         stripHitBuilder = std::make_shared<DefaultStripHitBuilder>(
             algoConfig["stripHitBuilder"]["params"]);
      } else {
         throw std::runtime_error("Unknown StripHitBuilder: " + stripAlgo);
      }

      // ClusterBuilder
      std::string clusterAlgo = algoConfig["clusterBuilder"]["name"];
      if (clusterAlgo == "ChargeWeightedCluster") {
         clusterBuilder = std::make_shared<ChargeWeightedCluster>(
             algoConfig["clusterBuilder"]["params"]);
      } else {
         throw std::runtime_error("Unknown ClusterBuilder: " + clusterAlgo);
      }
   }

   std::cout << "[Detector Init] " << m_name
             << " at (" << m_posX << ", " << m_posY << ", " << m_posZ
             << "), rotation (" << m_rotX << ", " << m_rotY << ", " << m_rotZ << ")\n";
}

void Detector::UpdateGeometry(double dx, double dy, double dz,
                              double dRotX, double dRotY, double dRotZ) {
   m_posX += dx;
   m_posY += dy;
   m_posZ += dz;
   m_rotX += dRotX;
   m_rotY += dRotY;
   m_rotZ += dRotZ;
}

LocalHit Detector::GlobalToLocal(const GlobalHit& globalHit) const {
   // 假设旋转角度较小，使用简单的平移和旋转转换
   LocalHit local;
   double dx = globalHit.x - m_posX;
   double dy = globalHit.y - m_posY;
   double dz = globalHit.z - m_posZ;

   // 这里可扩展为完整的旋转矩阵
   local.u = dx * cos(m_rotZ) + dy * sin(m_rotZ);

   // 如果是二维探测器
   if (m_readoutPlaneAngle.size() > 1) {
      local.v = dx * -sin(m_rotZ) + dy * cos(m_rotZ);
   }

   return local;
}

GlobalHit Detector::LocalToGlobal(const LocalHit& localHit) const {
   GlobalHit global;
   double u = localHit.u;
   double v = localHit.v.value_or(0.0);

   // 简单旋转+平移
   global.x = m_posX + u * cos(m_rotZ) - v * sin(m_rotZ);
   global.y = m_posY + u * sin(m_rotZ) + v * cos(m_rotZ);
   global.z = m_posZ;

   return global;
}