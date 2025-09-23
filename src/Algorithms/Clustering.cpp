#include "Algorithms/Clustering.h"
#include "DataModel.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include <algorithm>
#include <cmath>
#include <limits>

Clustering::Clustering(const json& config) {
    m_waveformConfig.loadFrom(config);
    m_clusterConfig.loadFrom(config);
    m_reconstructionConfig.loadFrom(config);
}

std::vector<std::vector<RawData>> Clustering::preClustering(const std::vector<RawData>& rawData) const {
    std::vector<std::vector<RawData>> preClusters;

    if (rawData.empty()) {
        return preClusters;
    }

    // 按类型和stripID排序
    std::vector<RawData> sortedData(rawData);
    std::sort(sortedData.begin(), sortedData.end(), [](const RawData& a, const RawData& b) {
        return (a.type == b.type) ? (a.stripID < b.stripID) : (a.type < b.type);
    });

    // 仅根据连续性进行分组
    std::vector<RawData> currentCluster;

    for (const auto& data : sortedData) {
        if (currentCluster.empty()) {
            // 第一个数据，直接添加
            currentCluster.push_back(data);
        } else {
            const RawData& lastData = currentCluster.back();

            bool isContinuous = (data.type == lastData.type) &&
                                (data.stripID <= lastData.stripID + m_clusterConfig.maxGap + 1);

            if (*std::max_element(data.adc.begin(), data.adc.end()) < m_waveformConfig.noiseThreshold) continue;

            if (isContinuous) {
                currentCluster.push_back(data);
            } else {
                if (currentCluster.size() >= m_clusterConfig.minClusterSize && currentCluster.size() <= m_clusterConfig.maxClusterSize) {
                    preClusters.push_back(currentCluster);
                }
                currentCluster.clear();
                currentCluster.push_back(data);
            }
        }
    }

    // 处理最后一个聚类
    if (!currentCluster.empty()) {
        preClusters.push_back(currentCluster);
    }

    return preClusters;
}

RecCluster Clustering::BuildCluster(const std::vector<RawData>& raws) {
    RecCluster aCluster;
    if (raws.empty()) return aCluster;

    aCluster.type = raws.front().type;
    for (const auto& raw : raws) {
        StripHit stripData;
        stripData.stripID = raw.stripID;
        stripData.type = raw.type;

        processWaveform(raw, stripData);
        if (!stripData.isValid) continue;

        aCluster.strips.push_back(std::move(stripData));
    }

    if (processCluster(aCluster))
        return aCluster;
    else
        return RecCluster();
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
    int peakAmp = 0, peakTime = 0, inducedCharge = 0;
    int firstOverTh = -1, lastOverTh = -1;
    const int noiseTh = m_waveformConfig.noiseThreshold;

    for (size_t i = 0; i < nSamples; ++i) {
        const int amp = waveform[i];

        if (amp > peakAmp) {
            peakAmp = amp;
            peakTime = static_cast<int>(i);
        }

        if (amp > noiseTh) {
            if (firstOverTh == -1) firstOverTh = static_cast<int>(i);
            lastOverTh = static_cast<int>(i);
            inducedCharge += amp;
        }
    }

    // 阈值以下，无有效信号
    if (firstOverTh == -1) {
        stripData.isValid = false;
        return;
    }

    // 2. 平滑滤波 + 构造TGraph
    static const double weights[5] = {0.25, 0.5, 1.0, 0.5, 0.25};

    TGraph signalGraph;
    signalGraph.Set(nSamples);

    for (size_t i = 0; i < nSamples; ++i) {
        double smoothedWave = 0.0;
        // 窗口 [-2,2]
        for (int m = -2; m <= 2; ++m) {
            int idx = static_cast<int>(i) + m;
            if (idx >= 0 && idx < static_cast<int>(nSamples)) {
                smoothedWave += waveform[idx] * weights[m + 2];
            }
        }
        signalGraph.SetPoint(i, static_cast<double>(i), smoothedWave);
    }

    // 3. 拟合
    static TF1 riseFunc("riseFunc", "[0]/(1+exp(-(x-[2])/[3]))+[1]", -5, 30);
    riseFunc.SetParameters(peakAmp, 0, (peakTime + firstOverTh) * 0.5, 15);

    riseFunc.SetParLimits(1, -30, 30);
    riseFunc.SetParLimits(0, peakAmp - 50, peakAmp + 50);

    const double fitStart = std::max(0.0, static_cast<double>(firstOverTh) - 5);
    const double fitEnd   = std::min(static_cast<double>(nSamples - 1), peakTime + 1.0);

    // 使用返回值判断拟合是否成功
    int fitStatus = signalGraph.Fit(&riseFunc, "QNR", "", fitStart, fitEnd);
    if (fitStatus != 0) {
        stripData.isValid = false;
        return;
    }

    // 4. 计算目标Y值
    const double fitMin = riseFunc.Eval(fitStart);
    const double fitMax = riseFunc.Eval(fitEnd);
    const double targetY = m_waveformConfig.cfdFraction * riseFunc.GetParameter(0) + riseFunc.GetParameter(1);

    if (targetY < fitMin || targetY > fitMax) {
        stripData.isValid = false;
        return;
    }

    // 5. 计算拟合时间
    const double fitTime = riseFunc.GetX(targetY, fitStart, fitEnd);

    // 计算rise time与time error
    const double tau = riseFunc.GetParameter(3);
    const double tauErr = riseFunc.GetParError(3);

    const double riseTime = tau * 2.0 * std::log(9.0);
    const double timeError = std::abs(std::log(1.0 / m_waveformConfig.cfdFraction - 1.0) * tauErr);

    if (timeError > 3) {
        stripData.isValid = false;
        return;
    }

    // 6. 填充结果
    stripData.amplitude    = peakAmp;
    stripData.charge       = inducedCharge;
    stripData.peakTime     = peakTime * m_waveformConfig.timePitch;
    stripData.time         = fitTime * m_waveformConfig.timePitch;
    stripData.riseTime     = riseTime * m_waveformConfig.timePitch;
    stripData.timeError    = timeError * m_waveformConfig.timePitch;
    stripData.isSaturated  = (peakAmp > m_waveformConfig.saturationLevel);
    stripData.isValid      = true;
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