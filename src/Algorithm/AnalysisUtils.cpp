#include "Algorithm/AnalysisUtils.h"
#include "Event/DetectorFrame.h"

#include <TFile.h>
#include <TH1D.h>
#include <TVirtualFFT.h>

#include <cmath>
#include <iostream>

using namespace std;

namespace AnalysisUtils {

Track FitTrack(const vector<TVector3>& hits) {
    Track t{};
    size_t n = hits.size();
    if (n < 2) return t;

    double sumX = 0, sumY = 0, sumZ = 0;
    for (auto& p : hits) {
        sumX += p.X();
        sumY += p.Y();
        sumZ += p.Z();
    }
    double mX = sumX / n, mY = sumY / n, mZ = sumZ / n;

    double Szz = 0, Szx = 0, Szy = 0;
    for (auto& p : hits) {
        double dz = p.Z() - mZ;
        Szz += dz * dz;
        Szx += dz * (p.X() - mX);
        Szy += dz * (p.Y() - mY);
    }

    t.kx = Szx / Szz;
    t.ky = Szy / Szz;
    t.bx = mX - t.kx * mZ;
    t.by = mY - t.ky * mZ;

    double chi2 = 0;
    for (auto& p : hits) {
        double dx = p.X() - (t.kx * p.Z() + t.bx);
        double dy = p.Y() - (t.ky * p.Z() + t.by);
        chi2 += dx * dx + dy * dy;
    }
    t.chi2 = chi2 / (2 * n - 4);
    return t;
}

std::pair<double, double> GetRange(const std::vector<double>& v) {
    if (v.size() < 3) return std::make_pair(0.0, 1.0);
    const double k = 4;

    // -----------------------------
    // 1st pass: raw mean / sigma
    // -----------------------------
    double sum1 = 0, sq1 = 0;
    for (double x : v) sum1 += x;
    double mean1 = sum1 / v.size();

    for (double x : v) sq1 += (x - mean1) * (x - mean1);
    double sigma1 = std::sqrt(sq1 / v.size());

    double low1 = mean1 - k * sigma1;
    double high1 = mean1 + k * sigma1;

    double sum2 = 0, sq2 = 0;
    int n2 = 0;

    for (double x : v) {
        if (x >= low1 && x <= high1) {
            sum2 += x;
            n2++;
        }
    }

    double mean2 = sum2 / n2;

    for (double x : v) {
        if (x >= low1 && x <= high1)
            sq2 += (x - mean2) * (x - mean2);
    }

    double sigma2 = std::sqrt(sq2 / n2);

    return std::make_pair(mean2 - k * sigma2, mean2 + k * sigma2);
}

void FFTAnalyzer(Cluster& cluster, Event& evt, int detID) {

    // ========= 静态累积对象 =========
    static TH1D* hSpectrum = nullptr;
    static int nAccumulated = 0;
    int nSamples = 26;

    // ---------- 1. 初始化直方图（只做一次） ----------
    if (!hSpectrum) {
        const double dt = 25;                     // ns
        const double df = 1.0 / (nSamples * dt);  // GHz
        const int nFreq = nSamples / 2 + 1;

        hSpectrum = new TH1D(
            "hSpectrum",
            "Accumulated Power Spectrum;Frequency;Power",
            nFreq, 0, nFreq * df);
    }

    for (int stripIndex : cluster.channelHitIndices) {
        ChannelHit hit = evt.detectorFramesMap[detID]->GetChannelHit(stripIndex);
        auto rawData = evt.detectorFramesMap[detID]->GetRawFromChannel(hit);
        const auto& waveform = rawData->adc;

        std::vector<double> timeSignal(nSamples);
        for (int i = 0; i < nSamples; ++i)
            timeSignal[i] = static_cast<double>(waveform[i]);

        TVirtualFFT* fft = TVirtualFFT::FFT(1, &nSamples, "R2C M K");
        fft->SetPoints(timeSignal.data());
        fft->Transform();

        const int nFreq = nSamples / 2 + 1;

        // ---------- 3. 累积功率谱 ----------
        double re = 0.0, im = 0.0;
        for (int k = 0; k < nFreq; ++k) {
            fft->GetPointComplex(k, re, im);

            double power = re * re + im * im;

            // 单边谱修正
            if (k != 0 && k != nSamples / 2)
                power *= 2.0;

            power /= (nSamples * nSamples);

            hSpectrum->AddBinContent(k + 1, power);
        }

        delete fft;
        ++nAccumulated;
    }

    // ---------- 2. FFT ----------

    if (nAccumulated % 10000 == 0) {
        hSpectrum->Scale(1.0 / nAccumulated);

        TFile outFile("spectrum.root", "RECREATE");
        hSpectrum->Write();
        outFile.Close();

        std::cout << "[Mode1] Power spectrum written, N = "
                  << nAccumulated << std::endl;
    }
}

double CalculateMean(const std::vector<double>& values) {
    if (values.empty()) return 0.0;

    double sum = 0.0;
    for (double val : values) {
        sum += val;
    }
    return sum / values.size();
}

double CalculateRMS(const std::vector<double>& values) {
    if (values.empty()) return 0.0;

    double mean = CalculateMean(values);
    double sumSq = 0.0;
    for (double val : values) {
        sumSq += (val - mean) * (val - mean);
    }
    return std::sqrt(sumSq / values.size());
}

}  // namespace AnalysisUtils
