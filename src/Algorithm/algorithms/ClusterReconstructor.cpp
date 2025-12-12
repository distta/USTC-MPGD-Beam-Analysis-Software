#include "algorithms/ClusterReconstructor.h"
#include "AlgorithmFactory.h"
#include "DataModel.h"
#include "DetectorFrame.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include <cmath>

REGISTER_ALGORITHM("ClusterReconstructor", ClusterReconstructor)

bool ClusterReconstructor::Process(DetectorFrame& frame) {
    auto& clusters = frame.GetMutableClusters();
    if (clusters.empty()) return false;

    const auto& stripHits = frame.StripHits();

    // 遍历每个Cluster，调用内部重建逻辑更新pos字段
    for (auto& cluster : clusters) {
        ReconstructPosition(cluster, stripHits);
    }

    return true;
}

void ClusterReconstructor::ReconstructPosition(Cluster& cluster, const std::vector<StripHit>& stripHits) {
    if (m_config.method == ReconstructionMethod::UTPC) {
        reconstructUTPC(cluster, stripHits);
    } else {
        reconstructChargeWeighted(cluster, stripHits);
    }
}

void ClusterReconstructor::reconstructChargeWeighted(Cluster& cluster, const std::vector<StripHit>& stripHits) {
    if (cluster.stripHitIndices.empty()) {
        cluster.pos = 0.0;
        return;
    }

    // 电荷加权位置重建（以stripID为单位）
    double weightedSum = 0.0;
    double totalCharge = 0.0;

    for (int idx : cluster.stripHitIndices) {
        const auto& strip = stripHits[idx];
        double weight = pow(strip.amp, 1);
        weightedSum += strip.ID * weight;
        totalCharge += weight;
    }

    if (totalCharge > 0) {
        cluster.pos = weightedSum / totalCharge;
    } else {
        int centerIdx = cluster.size / 2;
        cluster.pos = stripHits[cluster.stripHitIndices[centerIdx]].ID;
    }
}

void ClusterReconstructor::reconstructUTPC(Cluster& cluster, const std::vector<StripHit>& stripHits) {
    if (cluster.stripHitIndices.empty()) {
        cluster.pos = 0.0;
        return;
    }

    TGraphErrors* track = new TGraphErrors();

    double weightedSum = 0.0;
    double totalCharge = 0.0;

    for (int idx : cluster.stripHitIndices) {
        const auto& strip = stripHits[idx];
        double weight = pow(strip.charge, 1);
        weightedSum += strip.ID * weight;
        totalCharge += weight;
    }
    double ccRecPos = weightedSum / totalCharge;

    if (cluster.type == 0) {
        cluster.pos = ccRecPos;
        return;
    }

    if (cluster.stripHitIndices.size() < 3) {
        cluster.pos = ccRecPos;
        return;
    }

    double lorentzAngle = 0. / 180 * TMath::Pi();
    double velocity = 0.021;
    double gasGap = 5;

    static TF1 disCorFunc = TF1("fitFunc", "pol3", -2, 2);
    disCorFunc.SetParameters(82.6668, -8.664, 7.29361, -21.4624);

    static TF1 chCorFunc = TF1("chFitFunc", "pol5", 0, 1000);
    chCorFunc.SetParameters(4.27534, -0.115551, 0.000690723, -1.56102e-06, 1.56769e-09, -5.79673e-13);
    // ccRecPos -= 5 * tan(lorentzAngle) * 0.5;

    std::vector<StripHit> modifiedStrips;
    modifiedStrips.reserve(cluster.stripHitIndices.size());
    for (int idx : cluster.stripHitIndices) {
        modifiedStrips.push_back(stripHits[idx]);
    }

    // for (auto& strip : modifiedStrips) {
    //     double disCor = disCorFunc.Eval((strip.ID - ccRecPos) * 0.4);
    //     double chCor = strip.amp > 1000 ? chCorFunc.Eval(1000) : chCorFunc.Eval(strip.amp);
    //     strip.time = strip.time - disCor - chCor;
    // }

    // ---------------------- 剔除时间骤降的坏点 ----------------------
    std::sort(modifiedStrips.begin(), modifiedStrips.end(),
              [](const auto& a, const auto& b) { return a.ID < b.ID; });

    // double dropThreshold = 20.0;  // 时间骤降阈值，可以自己调

    // std::vector<StripHit> cleaned;
    // cleaned.reserve(modifiedStrips.size());

    // for (size_t i = 0; i < modifiedStrips.size(); i++) {
    //     if (i == 0) {
    //         cleaned.push_back(modifiedStrips[i]);
    //         continue;
    //     }

    //     double tPrev = cleaned.back().time;  // 已通过筛选的前一条
    //     double tNow = modifiedStrips[i].time;

    //     if (tNow < tPrev - dropThreshold) {
    //         continue;
    //     }

    //     cleaned.push_back(modifiedStrips[i]);
    // }

    for (auto& strip : modifiedStrips) {
        double deltaT = strip.time;
        double x0 = strip.ID - deltaT * velocity * TMath::Tan(lorentzAngle);
        double y0 = deltaT * velocity;

        if (y0 < -0.4 || y0 > gasGap + 0.4) continue;

        int index = track->GetN();
        track->SetPoint(index, x0, y0);
        track->SetPointError(index, 0, strip.timeError * velocity);
    }

    track->AddPoint(ccRecPos, gasGap / 2);  // add CC point
    int index = track->GetN();
    track->SetPoint(index, ccRecPos, gasGap / 2);
    track->SetPointError(index, 0.3, 0);

    if (track->GetN() < 3) {
        cluster.pos = ccRecPos;
        return;
    }

    track->Fit("pol1", "q");
    double k = track->GetFunction("pol1")->GetParameter(1);
    double b = track->GetFunction("pol1")->GetParameter(0);

    // track->SetMarkerStyle(20);
    // track->SetLineWidth(0);  // 0 = no line
    // TFile* outFile = new TFile("1.root", "UPDATE");
    // if (outFile) {
    //     outFile->cd();
    //     track->Write();
    // }
    // outFile->Close();

    cluster.pos = (gasGap / 2 - b) / k;
}
