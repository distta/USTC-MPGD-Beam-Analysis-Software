#include "algorithms/WaveformProcessor.h"
#include "AlgorithmFactory.h"
#include "DetectorFrame.h"
#include <TF1.h>
#include <TGraph.h>
#include <algorithm>
#include <cmath>

REGISTER_ALGORITHM("WaveformProcessor", WaveformProcessor)

bool WaveformProcessor::Process(DetectorFrame& frame) {
    const auto& rawData = frame.Raw();
    if (rawData.empty()) return false;

    auto& stripHits = frame.GetMutableStripHits();
    stripHits.clear();
    stripHits.reserve(rawData.size());

    // 处理每个RawData，生成StripHit
    for (size_t i = 0; i < rawData.size(); ++i) {
        StripHit sh = ProcessWaveform(rawData[i]);
        sh.rawIndices = static_cast<int>(i);
        stripHits.push_back(sh);
    }

    // 排序：首先按type升序，相同type内按stripID升序
    std::sort(stripHits.begin(), stripHits.end(),
              [](const StripHit& a, const StripHit& b) {
                  if (a.type != b.type) return a.type < b.type;
                  return a.ID < b.ID;
              });

    return true;
}

StripHit WaveformProcessor::ProcessWaveform(const RawData& rawData) {
    if (m_config.mode == "Fit") {
        return processWaveformLeadingEdgeFit(rawData);
    }
    return processWaveformDefault(rawData);
}

StripHit WaveformProcessor::processWaveformLeadingEdgeFit(const RawData& rawData) {

    const auto& waveform = rawData.adc;
    const size_t nSamples = waveform.size();

    StripHit stripData;
    stripData.isValid = true;
    stripData.ID = rawData.stripID;
    stripData.type = rawData.type;

    if (nSamples == 0) {
        stripData.isValid = false;
        return stripData;
    }

    // 1. 找峰值 + 计算电荷 + 找阈值上下限
    int peakAmp = 0, peakTime = 0, inducedCharge = 0;
    int firstOverTh = -1, lastOverTh = -1;
    const int noiseTh = m_config.noiseThreshold;

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
    if (peakAmp < noiseTh || peakAmp > m_config.saturationLevel) {
        stripData.isValid = false;
    }

    TGraph signalGraph;
    signalGraph.Set(nSamples);

    for (size_t i = 0; i < nSamples; ++i) {
        signalGraph.SetPoint(i, static_cast<double>(i), waveform[i]);
    }

    // 3. 拟合
    static TF1 riseFunc("riseFunc", "[0]/(1+exp(-(x-[2])/[3]))+[1]", -5, 30);
    riseFunc.SetParameters(peakAmp, 0, (peakTime + firstOverTh) * 0.5, 15);
    riseFunc.SetParLimits(1, -30, 30);
    riseFunc.SetParLimits(0, peakAmp - 50, peakAmp + 50);

    const double fitStart = std::max(0.0, static_cast<double>(firstOverTh) - 5);
    const double fitEnd = std::min(static_cast<double>(nSamples - 1), peakTime + 1.0);

    int fitStatus = signalGraph.Fit(&riseFunc, "QNR", "", fitStart, fitEnd);
    if (fitStatus != 0) {
        stripData.isValid = false;
    }

    // 4. 计算目标Y值
    const double fitMin = riseFunc.GetMinimum(fitStart, fitEnd);
    const double fitMax = riseFunc.GetMaximum(fitStart, fitEnd);
    double targetY = m_config.cfdFraction * riseFunc.GetParameter(0) + riseFunc.GetParameter(1);
    double fitTime;

    if (targetY < fitMin || targetY > fitMax) {
        fitTime = firstOverTh;
    } else {
        fitTime = riseFunc.GetX(targetY, fitStart, fitEnd);
    }

    if (fitTime < 2 || fitTime > 16)
        stripData.isValid = false;

    // 5. 计算拟合时间
    const double tau = riseFunc.GetParameter(3);
    const double tauErr = riseFunc.GetParError(3);
    const double riseTime = tau * 2.0 * std::log(9.0);
    const double timeError = std::abs(std::log(1.0 / m_config.cfdFraction - 1.0) * tauErr);

    // 6. 填充结果
    stripData.amp = peakAmp;
    stripData.charge = inducedCharge;
    stripData.peakTime = peakTime * m_config.timePitch;
    stripData.time = fitTime * m_config.timePitch;
    stripData.riseTime = riseTime * m_config.timePitch;
    stripData.timeError = timeError * m_config.timePitch;
    stripData.isSaturated = (peakAmp > m_config.saturationLevel);

    return stripData;
}

StripHit WaveformProcessor::processWaveformDefault(const RawData& rawData) {
    const auto& waveform = rawData.adc;
    const size_t nSamples = waveform.size();

    StripHit stripData;
    stripData.isValid = true;
    stripData.ID = rawData.stripID;
    stripData.type = rawData.type;

    if (nSamples == 0) {
        stripData.isValid = false;
        return stripData;
    }

    const int noiseTh = m_config.noiseThreshold;
    static const double weights[5] = {0.25, 0.5, 1.0, 0.5, 0.25};
    std::vector<double> smooth(nSamples, 0.0);

    // 平滑滤波
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
        smooth[i] = (sumw > 0) ? (s / sumw) : 0.0;
    }

    // 找峰值
    double peakAmpD = -1.0;
    int peakIdx = -1;
    for (size_t i = 0; i < nSamples; ++i) {
        if (smooth[i] > peakAmpD) {
            peakAmpD = smooth[i];
            peakIdx = static_cast<int>(i);
        }
    }
    int peakAmp = static_cast<int>(std::round(peakAmpD));

    if (peakAmp < m_config.noiseThreshold) {
        stripData.isValid = false;
    }

    // 计算基线!!!!
    double baseline = 0.0;
    size_t baselineSamples = std::min<size_t>(5, nSamples);
    double bsum = 0;
    for (size_t i = 0; i < baselineSamples; ++i) bsum += waveform[i];
    baseline = bsum / baselineSamples;

    // 计算电荷
    double inducedCharge = 0.0;
    for (size_t i = 0; i < nSamples; ++i) {
        if (waveform[i] > baseline) {
            inducedCharge += (waveform[i] - baseline);
        }
    }

    // CFD时间提取
    double targetY = baseline + m_config.cfdFraction * (peakAmpD - baseline);
    double minSmooth = *std::min_element(smooth.begin(), smooth.end());
    double maxSmooth = *std::max_element(smooth.begin(), smooth.end());

    if (targetY < minSmooth || targetY > maxSmooth) {
        stripData.isValid = false;
    }

    // 找穿越点
    int crossingIdx = -1;
    for (int i = 1; i <= peakIdx; ++i) {
        if (smooth[i] >= targetY && smooth[i - 1] < targetY) {
            crossingIdx = i;
            break;
        }
    }

    if (crossingIdx == -1) {
        stripData.isValid = false;
    }

    // 线性插值
    double fitTime = 0.0;
    if (stripData.isValid) {
        double y1 = smooth[crossingIdx - 1];
        double y2 = smooth[crossingIdx];
        double frac = (targetY - y1) / (y2 - y1);
        fitTime = (crossingIdx - 1) + frac;
    }

    // 填充结果
    stripData.amp = peakAmp - baseline;
    stripData.charge = inducedCharge;
    stripData.peakTime = peakIdx * m_config.timePitch;
    stripData.time = fitTime * m_config.timePitch;
    stripData.riseTime = 0;
    stripData.timeError = 0;
    stripData.isSaturated = (peakAmp > m_config.saturationLevel);

    return stripData;
}
