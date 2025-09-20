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

std::vector<RecCluster> Clustering::BuildCluster(const std::vector<RawData>& raw) {
    std::vector<RecCluster> clusters;
    if (raw.empty()) return clusters;

    std::vector<RawData> sortedRaw(raw);
    std::sort(sortedRaw.begin(), sortedRaw.end(), [](const RawData& a, const RawData& b) {
        return (a.type == b.type) ? (a.stripID < b.stripID) : (a.type < b.type);
    });

    RecCluster currentCluster;
    int lastStripID = -999;

    auto flushCluster = [&](RecCluster& cluster) {
        if (!cluster.strips.empty()) {
            cluster.type = cluster.strips.front().type;
            processCluster(cluster, clusters);
            cluster.strips.clear();
        }
    };

    for (const auto& raw : sortedRaw) {
        StripHit stripData;
        stripData.stripID = raw.stripID;
        stripData.type = raw.type;

        processWaveform(raw, stripData);
        if (!stripData.isValid) continue;

        if (!currentCluster.strips.empty() &&
            (stripData.stripID > lastStripID + m_clusterConfig.maxGap ||
             stripData.type != currentCluster.strips.front().type)) {
            flushCluster(currentCluster);
        }

        currentCluster.strips.push_back(std::move(stripData));
        lastStripID = raw.stripID;
    }

    flushCluster(currentCluster);

    return clusters;
}

void Clustering::processWaveform(const RawData& rawData, StripHit& stripData) {

    const auto& waveform = rawData.adc;
    if (waveform.empty()) {
        stripData.isValid = false;
        return;
    }

    int peakAmp = 0, peakTime = 0, inducedCharge = 0;
    int firstOverTh = -1, lastOverTh = -1;

    for (size_t i = 0; i < waveform.size(); ++i) {
        const int amp = waveform[i];

        if (amp > peakAmp) {
            peakAmp = amp;
            peakTime = static_cast<int>(i);
        }

        if (amp > m_waveformConfig.noiseThreshold) {
            if (firstOverTh == -1) firstOverTh = static_cast<int>(i);
            lastOverTh = static_cast<int>(i);
            inducedCharge += amp;
        }
    }

    if (firstOverTh == -1) {
        stripData.isValid = false;
        return;
    }

    TGraph signalGraph;
    static const double weights[5] = {0.25, 0.5, 1.0, 0.5, 0.25};
    for (size_t i = 0; i < waveform.size(); ++i) {
        double smoothedWave = 0.0;
        for (int m = -2; m <= 2; ++m) {
            const int idx = static_cast<int>(i) + m;
            if (idx >= 0 && idx < static_cast<int>(waveform.size())) {
                smoothedWave += waveform[idx] / weights[m + 2];
            }
        }
        signalGraph.AddPoint(static_cast<double>(i), smoothedWave);
    }

    TF1 riseFunc("riseFunc", "[0]/(1+exp(-(x-[2])/[3]))+[1]", -5, 30);
    riseFunc.SetParameter(0, peakAmp);
    riseFunc.SetParameter(1, 0);
    riseFunc.SetParameter(2, (peakTime + firstOverTh) * 0.5);
    riseFunc.SetParameter(3, 15);

    riseFunc.SetParLimits(1, -30, 30);
    riseFunc.SetParLimits(0, peakAmp - 50, peakAmp + 50);

    const double fitStart = -5;
    const double fitEnd = peakTime + 1;

    signalGraph.Fit(&riseFunc, "q", "", fitStart, fitEnd);

    double fitMin = riseFunc.Eval(fitStart);
    double fitMax = riseFunc.Eval(fitEnd);
    double TargetY = m_waveformConfig.cfdFraction * riseFunc.GetParameter(0) + riseFunc.GetParameter(1);
    if (TargetY < fitMin || TargetY > fitMax) {
        stripData.isValid = false;
        return;
    }

    const double fitTime = riseFunc.GetX(TargetY, fitStart, fitEnd);
    const double riseTime = riseFunc.GetParameter(3) * 2 * std::log(9);
    const double timeError = std::abs(
        std::log(1.0 / m_waveformConfig.cfdFraction - 1) * riseFunc.GetParError(3));

    if (timeError > 10) {
        stripData.isValid = false;
        return;
    }

    stripData.amplitude = peakAmp;
    stripData.charge = inducedCharge;
    stripData.peakTime = peakTime * m_waveformConfig.timePitch;
    stripData.time = fitTime * m_waveformConfig.timePitch;
    stripData.riseTime = riseTime * m_waveformConfig.timePitch;
    stripData.timeError = timeError * m_waveformConfig.timePitch;
    stripData.isSaturated = peakAmp > m_waveformConfig.saturationLevel;
    stripData.isValid = true;
}

void Clustering::processCluster(RecCluster& cluster, std::vector<RecCluster>& clusters) {
    cluster.size = static_cast<int>(cluster.strips.size());

    if (cluster.size < m_clusterConfig.minClusterSize ||
        cluster.size > m_clusterConfig.maxClusterSize) {
        return;
    }

    // 一次遍历计算所有属性
    cluster.charge = 0.0;
    cluster.maxAmplitude = 0.0;
    cluster.time = std::numeric_limits<double>::max();
    cluster.range = cluster.strips.back().stripID - cluster.strips.front().stripID;

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

    clusters.push_back(std::move(cluster));
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