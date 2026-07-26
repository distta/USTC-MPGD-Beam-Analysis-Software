#include "Algorithm/Analyzer/HitProcessor.h"
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

REGISTER_ALGORITHM("HitProcessor", HitProcessor)

bool HitProcessor::Process(DetectorFrame& frame) {
   const auto& rawData = frame.Raw();
   if (rawData.empty()) return false;

   auto& channelHits = frame.GetMutableChannelHits();
   channelHits.clear();
   channelHits.reserve(rawData.size());

   // 处理每个RawData，生成ChannelHit
   for (size_t i = 0; i < rawData.size(); ++i) {
      ChannelHit sh = ProcessHit(rawData[i]);
      sh.rawIndices = static_cast<int>(i);
      channelHits.push_back(sh);
   }

   // 排序：type -> row(id1) -> column/strip(id0)
   std::sort(channelHits.begin(), channelHits.end(),
             [](const ChannelHit& a, const ChannelHit& b) {
                if (a.type != b.type) return a.type < b.type;
                if (a.id1 != b.id1) return a.id1 < b.id1;
                return a.id0 < b.id0;
             });

   return true;
}

ChannelHit HitProcessor::ProcessHit(const RawData& rawData) {
   if (rawData.HasDirectHit()) {
      return processDirectHit(rawData);
   }
   if (m_config.mode == "Fit") {
      return processWaveformLeadingEdgeFit(rawData);
   } else if (m_config.mode == "Mode1") {
      return processWaveformMode1(rawData);
   }
   return processWaveformDefault(rawData);
}

ChannelHit HitProcessor::processDirectHit(const RawData& rawData) {
   ChannelHit channelData{};
   channelData.id0 = rawData.id0;
   channelData.id1 = rawData.id1;
   channelData.type = rawData.type;
   channelData.amp = rawData.adcValue;
   channelData.charge = rawData.adcValue;
   channelData.peakTime = 0;
   channelData.time = rawData.hitTimeNs;
   channelData.riseTime = 0.0;
   channelData.width = 0.0;
   channelData.timeError = 0.0;
   channelData.isSaturated =
       rawData.adcValue >= m_config.saturationLevel;
   channelData.isValid = rawData.adcValue >= m_config.threshold;
   return channelData;
}

ChannelHit HitProcessor::processWaveformLeadingEdgeFit(const RawData& rawData) {

   const auto& waveform = rawData.adc;
   const size_t nSamples = waveform.size();

   ChannelHit channelData{};
   channelData.isValid = true;
   channelData.id0 = rawData.id0;
   channelData.id1 = rawData.id1;
   channelData.type = rawData.type;

   if (nSamples == 0) {
      channelData.isValid = false;
      return channelData;
   }

   // 1. 找峰值 + 计算电荷 + 找阈值上下限
   int peakAmp = 0, peakTime = 0, inducedCharge = 0;
   int firstOverTh = -1, lastOverTh = -1;
   const int noiseTh = m_config.threshold;

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
      channelData.isValid = false;
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
      channelData.isValid = false;
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
      channelData.isValid = false;

   // 5. 计算拟合时间
   const double tau = riseFunc.GetParameter(3);
   const double tauErr = riseFunc.GetParError(3);
   const double riseTime = tau * 2.0 * std::log(9.0);
   const double crossingTimeError = riseFunc.GetParError(2);
   const double fractionDerivative =
       std::log(1.0 / m_config.cfdFraction - 1.0);
   const double timeError =
       std::hypot(crossingTimeError, fractionDerivative * tauErr);

   // 6. 填充结果
   channelData.amp = peakAmp;
   channelData.charge = inducedCharge;
   channelData.peakTime = peakTime * m_config.timePitch;
   channelData.time = fitTime * m_config.timePitch;
   channelData.riseTime = riseTime * m_config.timePitch;
   channelData.timeError = timeError * m_config.timePitch;
   channelData.isSaturated = (peakAmp > m_config.saturationLevel);

   return channelData;
}

ChannelHit HitProcessor::processWaveformDefault(const RawData& rawData) {
   const auto& waveform = rawData.adc;
   const size_t nSamples = waveform.size();

   ChannelHit channelData{};
   channelData.isValid = true;
   channelData.id0 = rawData.id0;
   channelData.id1 = rawData.id1;
   channelData.type = rawData.type;

   if (nSamples == 0) {
      channelData.isValid = false;
      return channelData;
   }

   const int noiseTh = m_config.threshold;

   // 找峰值
   int peakAmp = 0, peakTime = 0, inducedCharge = 0;
   int firstOverTh = -1, lastOverTh = -1;

   double baseline = 0.0;
//    size_t baselineSamples = std::min<size_t>(3, nSamples);
//    double bsum = 0;
//    for (size_t i = 0; i < baselineSamples; ++i) bsum += waveform[i];
//    baseline = bsum / baselineSamples;

   for (size_t i = 0; i < nSamples; ++i) {
      const int amp = waveform[i] - baseline;

      if (amp > peakAmp) {
         peakAmp = amp;
         peakTime = static_cast<int>(i);
      }

      if (amp > noiseTh) {
         if (firstOverTh == -1) firstOverTh = static_cast<int>(i);
         lastOverTh = static_cast<int>(i);
      }
   }

   // 阈值以下，无有效信号
   if (peakAmp < noiseTh || peakAmp > m_config.saturationLevel) {
      channelData.isValid = false;
   }

   if (peakAmp < m_config.threshold) {
      channelData.isValid = false;
   }

   // 计算电荷
   for (size_t i = 0; i < nSamples; ++i) {
      if (waveform[i] > noiseTh) {
         inducedCharge += (waveform[i] - noiseTh - baseline);
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
      channelData.isValid = false;
   } else {
      const double y1 = waveform[crossingIdx - 1];
      const double y2 = waveform[crossingIdx];
      if (y2 == y1) {
         channelData.isValid = false;
      } else {
         fitTime = (crossingIdx - 1) + (targetY - y1) / (y2 - y1);
      }
   }

   if (channelData.isValid &&
       (fitTime * m_config.timePitch < m_config.timeWindowStart ||
        fitTime * m_config.timePitch > m_config.timeWindowEnd))
      channelData.isValid = false;

   channelData.amp = peakAmp;
   channelData.charge = inducedCharge;
   channelData.peakTime = peakTime * m_config.timePitch;
   channelData.time = fitTime * m_config.timePitch;
   channelData.riseTime = 0;
   channelData.timeError = 0;
   channelData.isSaturated = (peakAmp > m_config.saturationLevel);

   return channelData;
}

ChannelHit HitProcessor::processWaveformMode1(const RawData& rawData) {

   const auto& waveform = rawData.adc;
   const size_t nSamples = waveform.size();

   ChannelHit channelData{};
   channelData.isValid = true;
   channelData.id0 = rawData.id0;
   channelData.id1 = rawData.id1;
   channelData.type = rawData.type;

   if (nSamples == 0) {
      channelData.isValid = false;
      return channelData;
   }

   // 1. 找峰值 + 计算电荷 + 找阈值上下限
   int peakAmp = 0, peakTime = 0, inducedCharge = 0;
   int firstOverTh = -1, lastOverTh = -1;
   const int noiseTh = m_config.threshold;

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
      channelData.isValid = false;
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
   //     channelData.isValid = false;  // 拟合失败
   //     return channelData;
   // }

   double n = sigFunc.GetParameter(4);
   double tau = sigFunc.GetParameter(2);
   double baseline = sigFunc.GetParameter(3);
   double deltaT = sigFunc.GetParameter(1);

   if (n < 0 || tau < 0 || deltaT < 0) {
      channelData.isValid = false;  // 拟合参数异常
      return channelData;
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
   channelData.amp = peakAmp;
   channelData.charge = inducedCharge;
   channelData.peakTime = peakTime * m_config.timePitch;
   channelData.time = fitTime * m_config.timePitch;
   channelData.width = (endTime - fitTime) * m_config.timePitch;
   channelData.riseTime = riseTime * m_config.timePitch;
   // channelData.timeError = timeError * m_config.timePitch;
   channelData.isSaturated = (peakAmp > m_config.saturationLevel);

   // if (std::isnan(fitTime)) {
   // std::cout << "Fit Time: " << fitTime << ", End Time: " << endTime << std::endl;
   // std::cout << "Peak Amp: " << peakAmpD + baseline << ", Baseline: " << peakAmp << std::endl;
   // std::cout << "Delta T: " << deltaT << ", Tau: " << tau << std::endl;
   // TFile outFile("waveform_fit.root", "UPDATE");
   // signalGraph.Write(("signalGraph_" + std::to_string(rawData.type) + "_" + std::to_string(rawData.id0)).c_str());
   // sigFunc.Write(("sigFunc_" + std::to_string(rawData.type) + "_" + std::to_string(rawData.id0)).c_str());
   // outFile.Close();
   // }
   return channelData;
}
