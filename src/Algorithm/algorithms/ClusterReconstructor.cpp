#include "algorithms/ClusterReconstructor.h"
#include "AlgorithmFactory.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include <cmath>

REGISTER_ALGORITHM("ClusterReconstructor", ClusterReconstructor)

void ClusterReconstructor::ReconstructPositions(std::vector<RecCluster>& recClusters) {
    for (auto& recCluster : recClusters) {
        for (auto& cluster : recCluster) {
            ReconstructPosition(cluster);
        }
    }
}

void ClusterReconstructor::ReconstructPosition(Cluster& cluster) {
    if (m_config.method == ReconstructionMethod::UTPC) {
        reconstructUTPC(cluster);
    } else {
        reconstructChargeWeighted(cluster);
    }
}

void ClusterReconstructor::reconstructChargeWeighted(Cluster& cluster) {
    if (cluster.strips.empty()) {
        cluster.pos = 0.0;
        return;
    }

    // 电荷加权位置重建（以stripID为单位）
    double weightedSum = 0.0;
    double totalCharge = 0.0;

    for (const auto& strip : cluster.strips) {
        double weight = pow(strip.charge, 1);
        weightedSum += strip.stripID * weight;
        totalCharge += weight;
    }

    if (totalCharge > 0) {
        // 直接输出stripID（不乘以pitch）
        cluster.pos = weightedSum / totalCharge;
    } else {
        // 如果总电荷为0，使用中心条的stripID
        int centerIdx = cluster.size / 2;
        cluster.pos = cluster.strips[centerIdx].stripID;
    }
}

void ClusterReconstructor::reconstructUTPC(Cluster& cluster) {
    if (cluster.strips.empty()) {
        cluster.pos = 0.0;
        return;
    }

    TGraphErrors* track = new TGraphErrors();

    double weightedSum = 0.0;
    double totalCharge = 0.0;

    for (const auto& strip : cluster.strips) {
        double weight = pow(strip.charge, 1);
        weightedSum += strip.stripID * weight;
        totalCharge += weight;
    }
    double ccRecPos = weightedSum / totalCharge;

    if (cluster.type == 0) {
        cluster.pos = ccRecPos;
        return;
    }

    if(cluster.strips.size() < 3) {
        cluster.pos = ccRecPos;
        return;
    }

    double lorentzAngle = 0. / 180 * TMath::Pi();
    double velocity = 0.02;
    double gasGap = 5;
    double tmin = 80;
    double bias = 0;

    // ccRecPos -= 5 * tan(lorentzAngle) * 0.5;
    for (const auto& strip : cluster.strips) {
        double deltaT = strip.time - tmin + bias;
        double x0 = strip.stripID - deltaT * velocity * TMath::Tan(lorentzAngle);
        double y0 = deltaT * velocity;

        int index = track->GetN();
        track->SetPoint(index, x0, y0);
        track->SetPointError(index, 0, strip.timeError * velocity);
    }

    track->AddPoint(ccRecPos, gasGap / 2);  // add CC point

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
