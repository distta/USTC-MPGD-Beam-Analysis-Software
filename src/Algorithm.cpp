#include "Algorithm.h"
#include "Config.h"
#include "DataModel.h"
#include <algorithm>
#include <cmath>
#include <limits>

Algorithm::Algorithm(const json& config) {
    m_waveformConfig.loadFrom(config);
    m_clusterConfig.loadFrom(config);
    m_reconstructionConfig.loadFrom(config);
}

std::vector<Cluster> Algorithm::BuildClusters(const std::vector<StripHit>& stripHits) {
    std::vector<Cluster> clusters;
    if (stripHits.empty()) return clusters;

    std::vector<StripHit> currentGroup;
    currentGroup.push_back(stripHits.front());

    for (size_t i = 1; i < stripHits.size(); ++i) {
        if (stripHits[i].stripID <= stripHits[i - 1].stripID + 1 + m_clusterConfig.maxGap) {
            currentGroup.push_back(stripHits[i]);
        } else {
            // 断开 -> 完成一个 cluster
            Cluster cluster;
            cluster.type = currentGroup.front().type;

            for (const auto& aHit : currentGroup) {
                cluster.strips.push_back(std::move(aHit));
            }

            if (processCluster(cluster)) {
                clusters.push_back(std::move(cluster));
            }

            currentGroup.clear();
            currentGroup.push_back(stripHits[i]);
        }
    }

    // 处理最后一个 cluster
    if (!currentGroup.empty()) {
        Cluster cluster;
        cluster.type = currentGroup.front().type;

        for (const auto& aHit : currentGroup) {
            cluster.strips.push_back(std::move(aHit));
        }

        if (processCluster(cluster)) {
            clusters.push_back(std::move(cluster));
        }
    }

    return clusters;
}

StripHit Algorithm::processFastWaveform(const RawData& rawData) {
    // 快速波形处理算法的简单实现
    const auto& waveform = rawData.adc;
    const size_t nSamples = waveform.size();

    if (nSamples == 0) {
        StripHit stripData;
        stripData.isValid = false;
        return stripData;
    }

    StripHit stripData;

    // 简单的峰值查找和电荷积分
    int peakAmp = 0;
    int peakIdx = 0;
    double inducedCharge = 0.0;
    const int noiseTh = m_waveformConfig.noiseThreshold;

    for (size_t i = 0; i < nSamples; ++i) {
        if (waveform[i] > peakAmp) {
            peakAmp = waveform[i];
            peakIdx = static_cast<int>(i);
        }

        // 积分电荷计算
        if (waveform[i] > noiseTh) {
            inducedCharge += (waveform[i] - noiseTh);
        }
    }

    // 基本的时间估算（峰值位置）
    double fitTimeSamples = static_cast<double>(peakIdx);

    stripData.amp = static_cast<double>(peakAmp);
    stripData.charge = inducedCharge;
    stripData.peakTime = peakIdx;
    stripData.time = fitTimeSamples * m_waveformConfig.timePitch;
    stripData.riseTime = 0.0;   // 快速算法不计算上升时间
    stripData.timeError = 0.0;  // 快速算法不计算时间误差
    stripData.isSaturated = (peakAmp > m_waveformConfig.saturationLevel);
    stripData.isValid = true;

    return stripData;
}

StripHit Algorithm::processWaveform(const RawData& rawData) {

    if (m_waveformConfig.mode == "Fast") {
        return processFastWaveform(rawData);
    }
    const auto& waveform = rawData.adc;
    const size_t nSamples = waveform.size();

    StripHit stripData;

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
        return stripData;
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
        stripData.isValid = false;
        return stripData;
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

    stripData.type = rawData.type;
    stripData.stripID = rawData.stripID;
    stripData.amp = peakAmpD - baseline;
    stripData.charge = inducedCharge;
    stripData.peakTime = peakIdx;
    stripData.time = fitTimeSamples * m_waveformConfig.timePitch;
    stripData.riseTime = riseTimeSamples * m_waveformConfig.timePitch;
    stripData.timeError = timeErrorSamples * m_waveformConfig.timePitch;
    stripData.isSaturated = (peakAmpD > m_waveformConfig.saturationLevel);

    stripData.isValid = true;
    return stripData;
}

bool Algorithm::processCluster(Cluster& cluster) {
    cluster.size = static_cast<int>(cluster.strips.size());

    if (cluster.size < m_clusterConfig.minClusterSize ||
        cluster.size > m_clusterConfig.maxClusterSize) {
        return false;
    }

    // 一次遍历计算所有属性
    cluster.charge = 0.0;
    cluster.maxAmp = 0.0;
    cluster.time = std::numeric_limits<double>::max();
    cluster.range = cluster.strips.back().stripID - cluster.strips.front().stripID + 1;

    for (const auto& strip : cluster.strips) {
        cluster.charge += strip.charge;
        cluster.maxAmp = std::max(cluster.maxAmp, strip.amp);
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

std::vector<RecCluster> Algorithm::MatchClusters(std::map<int, std::vector<Cluster>>& clustersByType) {
    std::vector<RecCluster> recClusters;

    if (clustersByType.empty()) return recClusters;

    const auto& refType = clustersByType.begin()->first;
    const auto& refClusters = clustersByType.begin()->second;

    int matchCounter = 0;

    for (const auto& refCl : refClusters) {
        RecCluster rec;
        rec.push_back(refCl);

        double refCharge = refCl.charge;

        for (const auto& [type, cls] : clustersByType) {
            if (type == refType) continue;
            if (cls.empty()) continue;

            const Cluster* best = nullptr;
            double bestDiff = std::numeric_limits<double>::max();

            for (const auto& c : cls) {
                double diff = std::fabs(c.charge - refCharge);
                if (diff < bestDiff) {
                    bestDiff = diff;
                    best = &c;
                }
            }
            if (best && bestDiff < m_clusterConfig.MaxChargeDiff) rec.push_back(*best);
        }

        if (rec.size() > 1) {
            for (auto& c : rec) c.matchID = matchCounter;
            matchCounter++;
            recClusters.push_back(std::move(rec));
        }
    }

    return recClusters;
}

void Algorithm::reconstructChargeWeighted(Cluster& cluster) {
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

void Algorithm::reconstructUTPC(Cluster& cluster) {
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

std::vector<RecCluster> Algorithm::Reconstruct(const std::vector<RawData>& raws) {
    std::vector<RecCluster> recClusters;
    if (raws.empty()) return recClusters;

    // -----------------------------
    // Step 1: 按 type 对 RawData 分组
    // -----------------------------
    std::map<int, std::vector<StripHit>> stripHitsByType;
    for (const auto& rd : raws) {
        StripHit aStripHit = processWaveform(rd);
        if (aStripHit.isValid) stripHitsByType[rd.type].push_back(aStripHit);
    }

    // -----------------------------
    // Step 2: 每种 type 进行聚类
    // -----------------------------
    std::map<int, std::vector<Cluster>> clustersByType;

    for (auto& [type, stripHits] : stripHitsByType) {

        std::sort(stripHits.begin(), stripHits.end(),
                  [](const StripHit& a, const StripHit& b) { return a.stripID < b.stripID; });

        clustersByType[type] = BuildClusters(stripHits);
    }

    return MatchClusters(clustersByType);
}
