#include "Algorithm/Analyzer/ClusterReconstructor.h"
#include "AlgorithmFactory.h"
#include "DataModel.h"
#include "Detector/PlanarPad.h"
#include "DetectorFrame.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include <chrono>
#include <cmath>

REGISTER_ALGORITHM("ClusterReconstructor", ClusterReconstructor)

bool ClusterReconstructor::Process(DetectorFrame& frame) {
    auto& clusters = frame.GetMutableClusters();
    if (clusters.empty()) return false;

    const auto& channelHits = frame.ChannelHits();

    // 遍历每个Cluster，调用内部重建逻辑更新pos字段
    for (auto& cluster : clusters) {
        if (m_detector && m_detector->GetPlanarPadConfig()) {
            reconstructPadChargeWeighted(cluster, channelHits);
        } else if (m_config.method == ReconstructionMethod::UTPC) {
            reconstructUTPC(cluster, channelHits, 0);
        } else if (m_config.method == ReconstructionMethod::RawUTPC) {
            reconstructRawUTPC(cluster, channelHits, 0);
        } else if (m_config.method == ReconstructionMethod::ChargeWeighted) {
            reconstructChargeWeighted(cluster, channelHits);
        } else {
            throw std::runtime_error("Unknown reconstruction method");
        }
    }

    return true;
}

void ClusterReconstructor::reconstructPadChargeWeighted(
    Cluster& cluster, const std::vector<ChannelHit>& channelHits) {
    const auto* detector = dynamic_cast<const PlanarPad*>(m_detector);
    if (!detector || cluster.channelHitIndices.empty()) return;

    TVector3 weightedPosition;
    double totalWeight = 0.0;
    for (int index : cluster.channelHitIndices) {
        if (index < 0 || index >= static_cast<int>(channelHits.size())) continue;
        const auto& hit = channelHits[index];
        if (!hit.HasID1()) continue;
        const auto& config = detector->GetPadConfig();
        const int row = hit.id1;
        const int column = hit.id0;
        if (row < 0 || row >= config.rows ||
            column < 0 || column >= config.columns) continue;
        const double value = m_config.weightSource == "amp"
                                 ? hit.amp
                                 : hit.charge;
        double weight = m_config.weightPower == 1.0
                            ? value
                            : std::pow(std::max(0.0, value),
                                       m_config.weightPower);
        if (!(weight > 0.0)) weight = 1.0;
        weightedPosition += weight * detector->PadPosition(row, column);
        totalWeight += weight;
    }
    if (totalWeight <= 0.0) return;
    cluster.localPosition = (1.0 / totalWeight) * weightedPosition;
    cluster.hasLocalPosition = true;
}

void ClusterReconstructor::reconstructChargeWeighted(Cluster& cluster, const std::vector<ChannelHit>& channelHits) {
    if (cluster.channelHitIndices.empty()) {
        cluster.pos = 0.0;
        return;
    }

    // 电荷加权位置重建（以stripID为单位）
    double weightedSum = 0.0;
    double totalCharge = 0.0;

    for (int idx : cluster.channelHitIndices) {
        const auto& strip = channelHits[idx];
        const double value = m_config.weightSource == "amp"
                                 ? strip.amp
                                 : strip.charge;
        const double weight = m_config.weightPower == 1.0
                                  ? value
                                  : std::pow(std::max(0.0, value),
                                             m_config.weightPower);
        weightedSum += strip.id0 * weight;
        totalCharge += weight;
    }

    if (totalCharge > 0) {
        cluster.pos = weightedSum / totalCharge;
    } else {
        int centerIdx = cluster.size / 2;
        cluster.pos = channelHits[cluster.channelHitIndices[centerIdx]].id0;
    }
}

double LinearFitIter(const std::vector<ChannelHit>& channelHits, double x0) {

    double lorentzAngle = 0. / 180 * TMath::Pi();
    double velocity = 0.021;
    double gasGap = 5;

    double weightedSum = 0.0;
    double totalCharge = 0.0;
    for (auto& strip : channelHits) {
        double weight = pow(strip.amp, 1);
        weightedSum += strip.id0 * weight;
        totalCharge += weight;
    }
    double ccRecPos = weightedSum / totalCharge;

    static TF1 disCorFunc = TF1("fitFunc", "pol3", -2, 2);
    disCorFunc.SetParameters(-81.77373196181657, 9.060528533349935, -6.620452354084467, 21.23345582870448);

    static TF1 TACorFun("TACorFun", "[0]+[1]/TMath::Power(x-[2],[3])", 0, 30);
    TACorFun.SetParameters(0, 10.31123357358735, -1.971738638688841e-8, 1);

    static TF1 tCorFunc = TF1("fitFunc", "pol3", 0, 400);
    tCorFunc.SetParameters(131.16578571441693, -2.057779313411427, 0.01023940912182441, -0.00001630064845924149);

    static TF1 tzFunc = TF1("fitFunc", "pol3", 0, 400);
    tzFunc.SetParameters(-0.8446133570883949, -0.0034434242618225676, 0.00015110013297083415, -2.607600601582346e-7);

    TGraphErrors* track = new TGraphErrors();

    for (int i = 0; i < channelHits.size(); i++) {
        auto& strip = channelHits[i];
        double x = strip.id0;

        if (i + 1 < channelHits.size() && channelHits[i + 1].id0 == strip.id0 + 1 && abs(channelHits[i + 1].time - strip.time) < 15) {
            x = (strip.id0 * strip.charge + channelHits[i + 1].id0 * channelHits[i + 1].charge) / (strip.charge + channelHits[i + 1].charge);
        }

        if (i - 1 >= 0 && channelHits[i - 1].id0 == strip.id0 - 1 && abs(channelHits[i - 1].time - strip.time) < 15) {
            x = (strip.id0 * strip.charge + channelHits[i - 1].id0 * channelHits[i - 1].charge) / (strip.charge + channelHits[i - 1].charge);
        }
        double dis = (x - x0) * 0.4;
        double disCor = disCorFunc.Eval(dis);
        double tCor = tCorFunc.Eval(strip.time);
        double y = (strip.time + disCor + tCor) * velocity;
        // y = tzFunc.Eval(strip.time);

        if (y < -0.4 || y > gasGap + 0.4) continue;

        track->SetPoint(i, x, y);
        track->SetPointError(i, 0, strip.timeError * velocity);
    }

    int index = track->GetN();
    track->SetPoint(index, ccRecPos, gasGap / 2);
    track->SetPointError(index, 0.26, 0);

    track->Fit("pol1", "q");
    double k = track->GetFunction("pol1")->GetParameter(1);
    double b = track->GetFunction("pol1")->GetParameter(0);

    delete track;

    return (gasGap / 2 - b) / k;
}

void ClusterReconstructor::reconstructUTPC(Cluster& cluster, const std::vector<ChannelHit>& channelHits, double t0) {
    if (cluster.channelHitIndices.empty()) {
        cluster.pos = 0.0;
        return;
    }

    double weightedSum = 0.0;
    double totalCharge = 0.0;

    for (int idx : cluster.channelHitIndices) {
        const auto& strip = channelHits[idx];
        const double value = m_config.weightSource == "amp"
                                 ? strip.amp
                                 : strip.charge;
        const double weight = m_config.weightPower == 1.0
                                  ? value
                                  : std::pow(std::max(0.0, value),
                                             m_config.weightPower);
        weightedSum += strip.id0 * weight;
        totalCharge += weight;
    }
    double ccRecPos = weightedSum / totalCharge;

    if (cluster.type == 0) {
        cluster.pos = ccRecPos;
        return;
    }

    if (cluster.channelHitIndices.size() < 3) {
        cluster.pos = ccRecPos;
        return;
    }

    std::vector<ChannelHit> modifiedStrips;
    modifiedStrips.reserve(cluster.channelHitIndices.size());
    for (int idx : cluster.channelHitIndices) {
        modifiedStrips.push_back(channelHits[idx]);
    }

    // ---------------------- 剔除时间骤降的坏点 ----------------------
    std::sort(modifiedStrips.begin(), modifiedStrips.end(),
              [](const auto& a, const auto& b) { return a.id0 < b.id0; });

    double pos = LinearFitIter(modifiedStrips, ccRecPos);
    double lastPos = -99;
    int iter = 1;
    while (abs(pos - lastPos) > 0.01 && iter < 5) {
        lastPos = pos;
        pos = (LinearFitIter(modifiedStrips, pos) + pos * iter) / (iter + 1);
        iter++;
    }

    cluster.pos = pos;
}

void ClusterReconstructor::reconstructRawUTPC(Cluster& cluster, const std::vector<ChannelHit>& channelHits, double t0) {
    if (cluster.channelHitIndices.empty()) {
        cluster.pos = 0.0;
        return;
    }

    TGraph* track = new TGraph();

    double weightedSum = 0.0;
    double totalCharge = 0.0;

    for (int idx : cluster.channelHitIndices) {
        const auto& strip = channelHits[idx];
        double weight = pow(strip.charge, 1);
        weightedSum += strip.id0 * weight;
        totalCharge += weight;
    }
    double ccRecPos = weightedSum / totalCharge;

    if (cluster.type == 0) {
        cluster.pos = ccRecPos;
        return;
    }

    if (cluster.channelHitIndices.size() < 3) {
        cluster.pos = ccRecPos;
        return;
    }

    double lorentzAngle = 0. / 180 * TMath::Pi();
    double velocity = 0.022754;
    double gasGap = 5;

    std::vector<ChannelHit> modifiedStrips;
    modifiedStrips.reserve(cluster.channelHitIndices.size());
    for (int idx : cluster.channelHitIndices) {
        modifiedStrips.push_back(channelHits[idx]);
    }

    std::sort(modifiedStrips.begin(), modifiedStrips.end(),
              [](const auto& a, const auto& b) { return a.id0 < b.id0; });

    for (auto& strip : modifiedStrips) {
        double deltaT = strip.time - 100;
        double x0 = strip.id0;
        double y0 = deltaT * velocity;

        int index = track->GetN();
        track->SetPoint(index, x0, y0);
    }

    if (track->GetN() < 3) {
        cluster.pos = ccRecPos;
        return;
    }

    track->Fit("pol1", "q");
    double k = track->GetFunction("pol1")->GetParameter(1);
    double b = track->GetFunction("pol1")->GetParameter(0);

    cluster.pos = (gasGap / 2 - b) / k;
}
