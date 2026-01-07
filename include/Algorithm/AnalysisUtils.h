#pragma once

#include "Event/DataModel.h"
#include <utility>
#include <vector>

// 前向声明
class Event;
class Cluster;

/**
 * @brief 分析工具函数命名空间
 *
 * 提供通用的分析算法和辅助函数，供脚本和分析模块使用
 */
namespace AnalysisUtils {

Track FitTrack(const std::vector<TVector3>& hits);

std::pair<double, double> GetRange(const std::vector<double>& v);

double CalculateMean(const std::vector<double>& values);

double CalculateRMS(const std::vector<double>& values);

void FFTAnalyzer(Cluster& cluster, Event& evt, int detID);

}  // namespace AnalysisUtils
