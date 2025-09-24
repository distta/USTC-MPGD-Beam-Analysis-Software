#include "Clustering.h"
#include "DataModel.h"
#include <algorithm>
#include <cmath>
#include <limits>

Clustering::Clustering(const json& config) {
    m_waveformConfig.loadFrom(config);
    m_clusterConfig.loadFrom(config);
    m_reconstructionConfig.loadFrom(config);
}

std::vector<RecCluster> Clustering::BuildClusters(const std::vector<RawData>& raws) {
    std::vector<RecCluster> clusters;
    if (raws.empty()) return clusters;

    std::vector<RawData> currentGroup;
    currentGroup.push_back(raws.front());

    for (size_t i = 1; i < raws.size(); ++i) {
        if (raws[i].stripID <= raws[i - 1].stripID + 1 + m_clusterConfig.maxGap) {
            currentGroup.push_back(raws[i]);
        } else {
            // 断开 -> 完成一个 cluster
            RecCluster cluster;
            cluster.type = currentGroup.front().type;

            for (const auto& raw : currentGroup) {
                StripHit stripData;
                stripData.stripID = raw.stripID;
                stripData.type = raw.type;
                processWaveform(raw, stripData);
                if (stripData.isValid)
                    cluster.strips.push_back(std::move(stripData));
            }

            if (processCluster(cluster)) {
                clusters.push_back(std::move(cluster));
            }

            // 开启新的 cluster
            currentGroup.clear();
            currentGroup.push_back(raws[i]);
        }
    }

    // 处理最后一个 cluster
    if (!currentGroup.empty()) {
        RecCluster cluster;
        cluster.type = currentGroup.front().type;

        for (const auto& raw : currentGroup) {
            StripHit stripData;
            stripData.stripID = raw.stripID;
            stripData.type = raw.type;
            processWaveform(raw, stripData);
            if (stripData.isValid)
                cluster.strips.push_back(std::move(stripData));
        }

        if (processCluster(cluster)) {
            clusters.push_back(std::move(cluster));
        }
    }

    return clusters;
}

void Clustering::processWaveform(const RawData& rawData, StripHit& stripData) {
    const auto& waveform = rawData.adc;
    const size_t nSamples = waveform.size();

    // 空波形直接返回
    if (nSamples == 0) {
        stripData.isValid = false;
        return;
    }

    // 1. 找峰值 + 计算电荷 + 找阈值上下限
    int peakTime = 0;
    double inducedCharge = 0.0;
    int firstOverTh = -1, lastOverTh = -1;
    const int noiseTh = m_waveformConfig.noiseThreshold;

    static const double weights[5] = {0.25, 0.5, 1.0, 0.5, 0.25};
    std::vector<double> smooth(nSamples, 0.0);

    for (size_t i = 0; i < nSamples; ++i) {
        double sumw = 0;
        double s = 0;
        for (int m = -2; m <= 2; ++m) {
            int idx = static_cast<int>(i) + m;
            if (idx >= 0 && idx < static_cast<int>(nSamples)) {
                double w = weights[m + 2];
                s += w * static_cast<double>(waveform[idx]);
                sumw += w;
            }
        }
        if (sumw > 0)
            smooth[i] = s / sumw;
        else
            smooth[i] = 0.0;
    }

    // --- 2) 找峰值（基于平滑波形） ---
    double peakAmpD = -1.0;
    int peakIdx = -1;
    for (size_t i = 0; i < nSamples; ++i) {
        if (smooth[i] > peakAmpD) {
            peakAmpD = smooth[i];
            peakIdx = static_cast<int>(i);
        }
    }
    int peakAmp = static_cast<int>(std::round(peakAmpD));

    // --- 3) 计算目标 CFD 电平（基线假设为 0 或可计算） ---
    double baseline = 0.0;  // 如果你需要更精确基线，可以在前若干个采样点平均
    size_t baselineSamples = std::min<size_t>(5, nSamples);
    double bsum = 0;
    for (size_t i = 0; i < baselineSamples; ++i) bsum += waveform[i];
    baseline = bsum / baselineSamples;

    // 计算电荷（简单积分）
    for (size_t i = 0; i < nSamples; ++i) {
        if (waveform[i] > baseline) {
            inducedCharge += (waveform[i] - baseline);
        }
    }

    double targetY = baseline + m_waveformConfig.cfdFraction * (peakAmpD - baseline);

    // 若 targetY 不在 [min(smooth区间), max(smooth区间)] 说明 CFD 不可用
    double minSmooth = *std::min_element(smooth.begin(), smooth.end());
    double maxSmooth = *std::max_element(smooth.begin(), smooth.end());
    if (targetY < minSmooth || targetY > maxSmooth) {
        stripData.isValid = false;
        return;
    }

    // --- 4) 从前向后找首次穿越点（在峰前） ---
    // 我们希望找到上升沿的 crossing 点：从第一个超过 threshold 的索引向前找线性插值位置
    int crossingIdx = -1;
    for (int i = 1; i <= peakIdx; ++i) {
        if (smooth[i] >= targetY && smooth[i - 1] < targetY) {
            crossingIdx = i;
            break;
        }
    }
    if (crossingIdx == -1) {
        // 没找到，尝试从峰后向前找（安全降级）
        for (int i = peakIdx + 1; i < static_cast<int>(nSamples); ++i) {
            if (smooth[i] >= targetY && smooth[i - 1] < targetY) {
                crossingIdx = i;
                break;
            }
        }
    }
    if (crossingIdx == -1) {
        stripData.isValid = false;
        return;
    }

    // --- 5) 线性插值得到更精确时间 ---
    double y1 = smooth[crossingIdx - 1];
    double y2 = smooth[crossingIdx];
    double x1 = static_cast<double>(crossingIdx - 1);
    double x2 = static_cast<double>(crossingIdx);
    double frac = 0.0;
    if (y2 != y1)
        frac = (targetY - y1) / (y2 - y1);
    else
        frac = 0.0;
    double fitTimeSamples = x1 + frac;  // time in sample units

    // --- 6) 估计 rise time (10%-90%) 简单近似 ---
    // 近似：rise time ≈ (log(9) * 2 * tau)；也可以用 (t90 - t10) 在样点上估算
    // 计算 t10 和 t90 by interpolation
    double y10 = baseline + 0.1 * (peakAmpD - baseline);
    double y90 = baseline + 0.9 * (peakAmpD - baseline);

    auto findCross = [&](double target) -> double {
        // search from rising edge start to peakIdx
        for (int i = 1; i <= peakIdx; ++i) {
            if (smooth[i] >= target && smooth[i - 1] < target) {
                double yy1 = smooth[i - 1], yy2 = smooth[i];
                double xx1 = i - 1, xx2 = i;
                double f = (yy2 - yy1) != 0 ? (target - yy1) / (yy2 - yy1) : 0.0;
                return xx1 + f;
            }
        }
        return -1.0;
    };

    double t10 = findCross(y10);
    double t90 = findCross(y90);
    double riseTimeSamples = (t10 >= 0 && t90 >= 0) ? (t90 - t10) : (m_waveformConfig.timePitch * 0.0);

    // --- 7) time error 简单估计（由采样间距和斜率估算） ---
    double slope = (y2 - y1) / (x2 - x1);  // in ADC per sample
    double timeErrorSamples = 0.0;
    if (slope != 0.0) {
        // propagate error: sigma_y ~ 1 ADC? (经验值) -> sigma_t = sigma_y / slope
        double sigmaY = 1.0;  // 你可以把噪声 RMS 放到配置中
        timeErrorSamples = std::abs(sigmaY / slope);
    }

    stripData.amplitude = peakAmpD - baseline;
    stripData.charge = inducedCharge;  // 仍然可用先前计算的 inducedCharge
    stripData.peakTime = peakIdx;
    stripData.time = fitTimeSamples * m_waveformConfig.timePitch;  // 转为真实时间单位
    stripData.riseTime = riseTimeSamples * m_waveformConfig.timePitch;
    stripData.timeError = timeErrorSamples * m_waveformConfig.timePitch;
    stripData.isSaturated = (peakAmpD > m_waveformConfig.saturationLevel);
    stripData.isValid = true;
}

bool Clustering::processCluster(RecCluster& cluster) {
    cluster.size = static_cast<int>(cluster.strips.size());

    if (cluster.size < m_clusterConfig.minClusterSize ||
        cluster.size > m_clusterConfig.maxClusterSize) {
        return false;
    }

    // 一次遍历计算所有属性
    cluster.charge = 0.0;
    cluster.maxAmplitude = 0.0;
    cluster.time = std::numeric_limits<double>::max();
    cluster.range = cluster.strips.back().stripID - cluster.strips.front().stripID + 1;

    for (const auto& strip : cluster.strips) {
        cluster.charge += strip.charge;
        cluster.maxAmplitude = std::max(cluster.maxAmplitude, strip.amplitude);
        cluster.time = std::min(cluster.time, strip.time);
    }

    // 根据方法进行位置重建
    switch (m_reconstructionConfig.method) {
        case ReconstructionMethod::ChargeWeighted:
            reconstructChargeWeighted(cluster);
            break;
        case ReconstructionMethod::UTPC:
            reconstructUTPC(cluster);
            break;
        default:
            reconstructChargeWeighted(cluster);
            break;
    }

    return true;
}

void Clustering::MatchClusters(std::vector<RecCluster>& clustersU, std::vector<RecCluster>& clustersV, std::vector<RecHit>& recHits) {

    for (auto& u : clustersU) {
        for (auto& v : clustersV) {

            RecHit hit;
            hit.cluster.push_back(u);
            hit.cluster.push_back(v);
            recHits.push_back(std::move(hit));
        }
    }
}

void Clustering::reconstructChargeWeighted(RecCluster& cluster) {
    if (cluster.charge <= 0.0) {
        cluster.pos = -999.0;
        return;
    }

    double weightedSum = 0.0;
    for (const auto& strip : cluster.strips) {
        weightedSum += strip.stripID * strip.charge;
    }

    cluster.pos = weightedSum / cluster.charge;
}

void Clustering::reconstructUTPC(RecCluster& cluster) {
    if (cluster.strips.empty()) {
        cluster.pos = -999.0;
        return;
    }

    double sum = 0.0;
    for (const auto& strip : cluster.strips) {
        sum += strip.stripID;
    }

    cluster.pos = sum / cluster.strips.size();
}

std::vector<RecHit> Clustering::Reconstruction(const std::vector<RawData>& raws) {
    std::vector<RecHit> recHits;
    if (raws.empty()) return recHits;

    // -----------------------------
    // Step 1: 按 type 对 RawData 分组
    // -----------------------------
    std::map<int, std::vector<RawData>> typeMap;
    for (const auto& rd : raws) {
        typeMap[rd.type].push_back(rd);
    }

    // -----------------------------
    // Step 2: 每种 type 进行聚类
    // -----------------------------
    std::map<int, std::vector<RecCluster>> clustersByType;

    for (auto& [type, rawGroup] : typeMap) {
        // 按条号排序，保证聚类算法稳定
        std::sort(rawGroup.begin(), rawGroup.end(),
                  [](const RawData& a, const RawData& b) { return a.stripID < b.stripID; });

        // 调用 Clustering 完成该 type 内部的聚类
        auto clusters = BuildClusters(rawGroup);
        clustersByType[type] = std::move(clusters);
    }

    // -----------------------------
    // Step 3: 根据 type 数量生成 RecHit
    // -----------------------------
    if (clustersByType.empty()) return recHits;

    if (clustersByType.size() == 1) {
        // --- 1D 探测器，直接输出每个 cluster 为一个 RecHit ---
        const auto& clusters = clustersByType.begin()->second;
        for (const auto& cluster : clusters) {
            RecHit hit;
            hit.cluster.push_back(cluster);
            recHits.push_back(std::move(hit));
        }
    } else {
        // --- 2D 探测器，需要进行匹配，例如 U-V 匹配 ---
        // 假设只有两种类型：U / V
        auto it = clustersByType.begin();
        auto& clustersU = it->second;
        ++it;
        auto& clustersV = it->second;

        MatchClusters(clustersU, clustersV, recHits);
    }

    return recHits;
}
