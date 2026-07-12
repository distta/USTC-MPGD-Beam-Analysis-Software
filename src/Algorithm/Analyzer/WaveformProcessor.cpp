#include "Algorithm/Analyzer/WaveformProcessor.h"
#include "AlgorithmFactory.h"
#include "DetectorFrame.h"
#include "TVirtualFFT.h"
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

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
    } else if (m_config.mode == "Mode1") {
        return processWaveformMode1(rawData);
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
    // riseFunc.SetParLimits(1, -30, 30);
    // riseFunc.SetParLimits(0, peakAmp - 50, peakAmp + 50);

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

    if (fitTime * m_config.timePitch < m_config.timeWindowStart || fitTime * m_config.timePitch > m_config.timeWindowEnd)
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

    // 找峰值
    int peakAmp = 0, peakTime = 0, inducedCharge = 0;
    int firstOverTh = -1, lastOverTh = -1;

    double baseline = 0.0;
    // size_t baselineSamples = std::min<size_t>(5, nSamples);
    // double bsum = 0;
    // for (size_t i = 0; i < baselineSamples; ++i) bsum += waveform[i];
    // baseline = bsum / baselineSamples;

    for (size_t i = 0; i < nSamples; ++i) {
        const int amp = waveform[i] - baseline;

        if (amp > peakAmp) {
            peakAmp = amp;
            peakTime = static_cast<int>(i);
        }

        if (amp > noiseTh) {
            if (firstOverTh == -1) firstOverTh = static_cast<int>(i);
            lastOverTh = static_cast<int>(i);
            inducedCharge += amp - baseline;
        }
    }

    // 阈值以下，无有效信号
    if (peakAmp < noiseTh || peakAmp > m_config.saturationLevel) {
        stripData.isValid = false;
    }

    if (peakAmp < m_config.noiseThreshold) {
        stripData.isValid = false;
    }

    // 计算电荷
    for (size_t i = 0; i < nSamples; ++i) {
        if (waveform[i] > noiseTh) {
            inducedCharge += (waveform[i] - noiseTh);
        }
    }

    // Constant-fraction leading-edge time with sample-to-sample interpolation.
    const double targetY = baseline + m_config.cfdFraction * peakAmp;
    int crossingIdx = -1;
    for (int i = 1; i <= peakTime; ++i) {
        if (waveform[i - 1] < targetY && waveform[i] >= targetY) {
            crossingIdx = i;
            break;
        }
    }

    double fitTime = std::numeric_limits<double>::quiet_NaN();
    if (crossingIdx < 0) {
        stripData.isValid = false;
    } else {
        const double y1 = waveform[crossingIdx - 1];
        const double y2 = waveform[crossingIdx];
        if (y2 == y1) {
            stripData.isValid = false;
        } else {
            fitTime = (crossingIdx - 1) + (targetY - y1) / (y2 - y1);
        }
    }

    if (stripData.isValid &&
        (fitTime * m_config.timePitch < m_config.timeWindowStart ||
         fitTime * m_config.timePitch > m_config.timeWindowEnd))
        stripData.isValid = false;

    stripData.amp = peakAmp;
    stripData.charge = inducedCharge;
    stripData.peakTime = peakTime * m_config.timePitch;
    stripData.time = fitTime * m_config.timePitch;
    stripData.riseTime = 0;
    stripData.timeError = 0;
    stripData.isSaturated = (peakAmp > m_config.saturationLevel);

    return stripData;
}

StripHit WaveformProcessor::processWaveformMode1(const RawData& rawData) {

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

    const double fitStart = std::max(0.0, static_cast<double>(firstOverTh) - 5);

    static TF1 sigFunc("sigFunc", "[0]*(x-[1]>0?TMath::Power(x-[1],[4])*(exp(-(x-[1])/[2])):0) + [3]", 0, 30);
    sigFunc.SetParameters(peakAmp, firstOverTh - 1, double(peakTime) / firstOverTh, 0, 1.0);

    int fitSigStatus = signalGraph.Fit(&sigFunc, "QNR", "", fitStart, 30);
    // if (fitSigStatus != 0) {
    //     stripData.isValid = false;  // 拟合失败
    //     return stripData;
    // }

    double n = sigFunc.GetParameter(4);
    double tau = sigFunc.GetParameter(2);
    double baseline = sigFunc.GetParameter(3);
    double deltaT = sigFunc.GetParameter(1);

    if (n < 0 || tau < 0 || deltaT < 0) {
        stripData.isValid = false;  // 拟合参数异常
        return stripData;
    }

    double peakAmpD = sigFunc.GetParameter(0) * std::pow(n * tau, n) * std::exp(-n);
    double endTime = sigFunc.GetX(m_config.cfdFraction * peakAmpD + baseline, deltaT + n * tau, 100);
    double fitTime = sigFunc.GetX(m_config.cfdFraction * peakAmpD + baseline, 0, deltaT + n * tau);
    double riseTime = sigFunc.GetX((1 - m_config.cfdFraction) * peakAmpD + baseline, 0, deltaT + n * tau) - fitTime;
    // std::cout << "Fit Time: " << fitTime << ", End Time: " << endTime << std::endl;
    // std::cout << "Peak Amp: " << peakAmpD + baseline << ", Baseline: " << peakAmp << std::endl;
    // std::cout << "n: " << n << " Delta T: " << deltaT << ", Tau: " << tau << " baseline: " << baseline << std::endl;

    // 5. 计算拟合时间
    // const double tau = riseFunc.GetParameter(3);
    // const double tauErr = riseFunc.GetParError(3);
    // const double riseTime = tau * 2.0 * std::log(9.0);
    // const double timeError = std::abs(std::log(1.0 / m_config.cfdFraction - 1.0) * tauErr);

    // 6. 填充结果
    stripData.amp = peakAmp;
    stripData.charge = inducedCharge;
    stripData.peakTime = peakTime * m_config.timePitch;
    stripData.time = fitTime * m_config.timePitch;
    stripData.width = (endTime - fitTime) * m_config.timePitch;
    stripData.riseTime = riseTime * m_config.timePitch;
    // stripData.timeError = timeError * m_config.timePitch;
    stripData.isSaturated = (peakAmp > m_config.saturationLevel);

    // if (std::isnan(fitTime)) {
    // std::cout << "Fit Time: " << fitTime << ", End Time: " << endTime << std::endl;
    // std::cout << "Peak Amp: " << peakAmpD + baseline << ", Baseline: " << peakAmp << std::endl;
    // std::cout << "Delta T: " << deltaT << ", Tau: " << tau << std::endl;
    // TFile outFile("waveform_fit.root", "UPDATE");
    // signalGraph.Write(("signalGraph_" + std::to_string(rawData.type) + "_" + std::to_string(rawData.stripID)).c_str());
    // sigFunc.Write(("sigFunc_" + std::to_string(rawData.type) + "_" + std::to_string(rawData.stripID)).c_str());
    // outFile.Close();
    // }
    return stripData;
}
