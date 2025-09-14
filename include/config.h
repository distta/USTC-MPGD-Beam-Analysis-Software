#ifndef CONFIG_H
#define CONFIG_H

#include <map>
#include <nlohmann/json.hpp>
#include <string>
#include <variant>

struct GeometryConfig {
   std::string shape = "planar";                           // 几何形状
   std::map<int, double> pitchMap = {{0, 0.4}, {1, 0.4}};  // 间距映射
   std::map<int, double> angleMap = {{0, 0.0}, {1, 90}};   // 角度映射

   // 平面几何参数
   struct PlanarParams {
      double length = 0.0;
      double width = 0.0;
      PlanarParams() = default;
      PlanarParams(double l, double w) : length(l), width(w) {}
   };

   // 圆柱几何参数
   struct CylindricalParams {
      double radius = 0.0;
      double height = 0.0;
      CylindricalParams() = default;
      CylindricalParams(double r, double h) : radius(r), height(h) {}
   };

   // 使用 variant 存储参数
   using GeometryVariant = std::variant<PlanarParams, CylindricalParams>;
   GeometryVariant params{PlanarParams()};

   // 默认构造函数
   GeometryConfig() = default;

   // 从 JSON 构造
   GeometryConfig(const nlohmann::json& j) {
      shape = j.value("shape", "planar");

      // 根据形状设置参数
      if (shape == "cylindrical") {
         double r = j.value("radius", 0.0);
         double h = j.value("height", 0.0);
         params = CylindricalParams(r, h);
      } else {
         double l = j.value("length", 0.0);
         double w = j.value("width", 0.0);
         params = PlanarParams(l, w);
      }
   }

   // 获取平面参数
   PlanarParams getPlanarParams() const {
      return std::get<PlanarParams>(params);
   }

   // 获取圆柱参数
   CylindricalParams getCylindricalParams() const {
      return std::get<CylindricalParams>(params);
   }
};

// 合并的探测器处理器配置
struct ProcessorConfig {
   // 波形分析配置
   struct WaveformConfig {
      double cfdFraction = 0.5;         // CFD时间提取分数
      double noiseThreshold = 5.0;      // 噪声阈值
      double saturationLevel = 4000.0;  // 饱和电平
      double riseTimeStart = 0.1;       // 升时间计算起始分数
      double riseTimeEnd = 0.9;         // 升时间计算结束分数
   } waveform;

   // 聚类分析配置
   struct ClusterConfig {
      int maxGap = 2;           // 最大间隙
      int minClusterSize = 1;   // 最小聚类大小
      int maxClusterSize = 20;  // 最大聚类大小
   } cluster;
};

#endif  // CONFIG_H