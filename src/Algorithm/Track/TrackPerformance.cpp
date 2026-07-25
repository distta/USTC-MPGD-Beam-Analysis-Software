#include "Algorithm/Track/TrackPerformance.h"
#include "Detector/Detector.h"
#include "Event/DetectorFrame.h"

#include <TCanvas.h>
#include <TDirectory.h>
#include <TF1.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TPaveText.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

namespace Tracking {
namespace {

struct HitMapBinning {
    int binsX = 128;
    int binsY = 128;
    double xMin = 0.0;
    double xMax = 1.0;
    double yMin = 0.0;
    double yMax = 1.0;
};

HitMapBinning MakeHitMapBinning(const Detector& detector) {
    HitMapBinning result;
    if (const auto* config = detector.GetPlanarConfig()) {
        const auto xPitch = config->readoutPlanePitch.find(0);
        const auto xStrips = config->readoutPlaneStripNumber.find(0);
        if (xPitch != config->readoutPlanePitch.end() &&
            xStrips != config->readoutPlaneStripNumber.end()) {
            result.xMax = xPitch->second * xStrips->second;
        }
        const auto yPitch = config->readoutPlanePitch.find(1);
        const auto yStrips = config->readoutPlaneStripNumber.find(1);
        if (yPitch != config->readoutPlanePitch.end() &&
            yStrips != config->readoutPlaneStripNumber.end()) {
            result.yMax = yPitch->second * yStrips->second;
        }
    } else if (const auto* config = detector.GetPlanarPadConfig()) {
        result.xMin = -0.5 * config->pitchX;
        result.xMax = (config->columns - 0.5) * config->pitchX;
        result.yMin = -0.5 * config->pitchY;
        result.yMax = (config->rows - 0.5) * config->pitchY;
    }
    return result;
}

struct GaussianFit {
    double mean{};
    double sigma{};
    double chi2Ndf{};
};

GaussianFit FitSingleGaussian(TH1D* histogram) {
    GaussianFit result;
    if (!histogram || histogram->GetEntries() < 20) return result;
    histogram->Fit("gaus", "Q0");
    auto* gaussian = histogram->GetFunction("gaus");
    if (!gaussian) return result;
    gaussian->ResetBit(TF1::kNotDraw);
    result.mean = gaussian->GetParameter(1);
    result.sigma = std::abs(gaussian->GetParameter(2));
    result.chi2Ndf = gaussian->GetNDF() > 0 ? gaussian->GetChisquare() / gaussian->GetNDF() : 0.0;
    return result;
}

void DrawGaussianFit(TH1D* histogram, const GaussianFit& fit) {
    histogram->SetStats(false);
    histogram->Draw();
    if (auto* gaussian = histogram->GetFunction("gaus")) {
        gaussian->SetLineColor(kRed + 1);
        gaussian->SetLineWidth(2);
        gaussian->Draw("SAME");
    }
    auto* label = new TPaveText(0.62, 0.72, 0.88, 0.88, "NDC");
    label->SetFillColor(kWhite);
    label->SetBorderSize(1);
    label->SetTextAlign(12);
    std::ostringstream mean, sigma, quality;
    mean << "Gaussian #mu = " << std::setprecision(4) << fit.mean << " mm";
    sigma << "Gaussian #sigma = " << std::setprecision(4) << fit.sigma << " mm";
    quality << "#chi^{2}/ndf = " << std::setprecision(4) << fit.chi2Ndf;
    label->AddText(mean.str().c_str());
    label->AddText(sigma.str().c_str());
    label->AddText(quality.str().c_str());
    label->Draw();
}

}  // namespace

PerformanceAnalyzer::PerformanceAnalyzer(TDirectory* output,
                                         std::vector<std::shared_ptr<Detector>> detectors,
                                         std::vector<std::shared_ptr<Detector>> referenceDetectors,
                                         const Config& config,
                                         double residualRange)
    : m_output(output),
      m_detectors(std::move(detectors)),
      m_referenceDetectors(std::move(referenceDetectors)),
      m_config(config) {
    auto* resolutionDirectory = m_output->mkdir("Resolution");
    {
        TDirectory::TContext resolutionContext(resolutionDirectory);
        m_commonEquivalentHitX = new TH1D(
            "hCommonEquivalentHitX",
            "Common equivalent hit residual X;Equivalent hit residual X [mm];Entries",
            240, -residualRange, residualRange);
        m_commonEquivalentHitY = new TH1D(
            "hCommonEquivalentHitY",
            "Common equivalent hit residual Y;Equivalent hit residual Y [mm];Entries",
            240, -residualRange, residualRange);
    }

    for (const auto& detector : m_detectors) {
        const int id = detector->GetID();
        auto* detectorDir = m_output->mkdir(("Detector_" + std::to_string(id)).c_str());
        auto* hitsDir = detectorDir->mkdir("Hits");
        auto* residualDir = detectorDir->mkdir("Residuals");
        auto& h = m_detectorHistograms[id];
        {
            TDirectory::TContext context(hitsDir);
            h.clusterMultiplicity = new TH1D(
                "hClusterMultiplicity",
                "Cluster multiplicity per selected event;Clusters/event;Events",
                101, -0.5, 100.5);
            const auto binning = MakeHitMapBinning(*detector);
            h.hitMap = new TH2D(
                "hHitMap", "Selected local-hit map;Local X [mm];Local Y [mm]",
                binning.binsX, binning.xMin, binning.xMax,
                binning.binsY, binning.yMin, binning.yMax);
        }
        {
            TDirectory::TContext context(residualDir);
            h.residualX = new TH1D("hUnbiasedResidualX", "Excluded Tracker hit - remaining-Tracker fit, X;X_{excluded}-X_{remaining fit} [mm];Entries", 240, -residualRange, residualRange);
            h.residualY = new TH1D("hUnbiasedResidualY", "Excluded Tracker hit - remaining-Tracker fit, Y;Y_{excluded}-Y_{remaining fit} [mm];Entries", 240, -residualRange, residualRange);
        }
    }
}

void PerformanceAnalyzer::RecordEvent(const Event& event) {
    std::vector<TVector3> globalHits;
    globalHits.reserve(m_detectors.size());
    for (const auto& detector : m_detectors) {
        const auto& frame = event.detectorFramesMap.at(detector->GetID());
        if (frame->LocalHits().size() != 1) return;
        globalHits.push_back(detector->LocalToGlobal(frame->LocalHits().front().localPos));
    }

    for (size_t detectorIndex = 0; detectorIndex < m_detectors.size(); ++detectorIndex) {
        const auto& detector = m_detectors[detectorIndex];
        const auto& frame = event.detectorFramesMap.at(detector->GetID());
        const auto& hit = frame->LocalHits().front();
        auto& h = m_detectorHistograms[detector->GetID()];
        h.clusterMultiplicity->Fill(frame->Clusters().size());
        h.hitMap->Fill(hit.localPos.X(), hit.localPos.Y());
        std::vector<TVector3> otherHits;
        otherHits.reserve(m_detectors.size() - 1);
        for (size_t otherIndex = 0; otherIndex < globalHits.size(); ++otherIndex) {
            if (otherIndex != detectorIndex) otherHits.push_back(globalHits[otherIndex]);
        }
        if (otherHits.size() < 2) continue;

        const auto fit = FitWeighted(otherHits, m_config.resolutionX, m_config.resolutionY);
        const auto predicted = detector->GlobalToLocal(detector->CalcHitFromTrack(fit));
        const double rx = hit.localPos.X() - predicted.X();
        const double ry = hit.localPos.Y() - predicted.Y();

        double sumZ = 0.0, sumZZ = 0.0;
        for (const auto& otherHit : otherHits) {
            sumZ += otherHit.Z();
            sumZZ += otherHit.Z() * otherHit.Z();
        }
        const double n = static_cast<double>(otherHits.size());
        const double meanZ = sumZ / n;
        const double szz = sumZZ - sumZ * sumZ / n;
        double predictionVarianceScale = 0.0;
        if (szz > 1e-12) {
            const double dz = globalHits[detectorIndex].Z() - meanZ;
            predictionVarianceScale = 1.0 / n + dz * dz / szz;
        }
        const double scale = std::sqrt(1.0 + predictionVarianceScale);
        h.residualX->Fill(rx);
        h.residualY->Fill(ry);
        m_commonEquivalentHitX->Fill(rx / scale);
        m_commonEquivalentHitY->Fill(ry / scale);
    }
}

std::pair<double, double> PerformanceAnalyzer::EstimateHitResolution() {
    const auto fitX = FitSingleGaussian(m_commonEquivalentHitX);
    const auto fitY = FitSingleGaussian(m_commonEquivalentHitY);
    return {fitX.sigma, fitY.sigma};
}

void PerformanceAnalyzer::Write() {
    for (auto& [id, h] : m_detectorHistograms) {
        const auto residualFitX = FitSingleGaussian(h.residualX);
        const auto residualFitY = FitSingleGaussian(h.residualY);
        TDirectory::TContext residualContext(h.residualX->GetDirectory());
        TCanvas residualCanvas(
            ("cLeaveOneOutResidual_Detector" + std::to_string(id)).c_str(),
            ("Remaining trackers vs excluded Tracker " + std::to_string(id)).c_str(),
            1400, 600);
        residualCanvas.Divide(2, 1);
        residualCanvas.cd(1);
        DrawGaussianFit(h.residualX, residualFitX);
        residualCanvas.cd(2);
        DrawGaussianFit(h.residualY, residualFitY);
        residualCanvas.Write();
    }

    const auto commonFitX = FitSingleGaussian(m_commonEquivalentHitX);
    const auto commonFitY = FitSingleGaussian(m_commonEquivalentHitY);
    const double sigmaHitX = commonFitX.sigma;
    const double sigmaHitY = commonFitY.sigma;
    if (sigmaHitX <= 0.0 || sigmaHitY <= 0.0 || m_detectors.size() < 3) {
        m_output->Write();
        return;
    }

    double meanZ = 0.0;
    for (const auto& detector : m_detectors) meanZ += detector->GetPos().Z();
    meanZ /= static_cast<double>(m_detectors.size());
    double szz = 0.0;
    for (const auto& detector : m_detectors) {
        const double dz = detector->GetPos().Z() - meanZ;
        szz += dz * dz;
    }
    if (szz <= 1e-12) {
        m_output->Write();
        return;
    }

    const auto pointingScale = [&](double z) {
        const double dz = z - meanZ;
        return std::sqrt(1.0 / static_cast<double>(m_detectors.size()) + dz * dz / szz);
    };

    double zMin = m_detectors.front()->GetPos().Z();
    double zMax = zMin;
    for (const auto& detector : m_detectors) {
        zMin = std::min(zMin, detector->GetPos().Z());
        zMax = std::max(zMax, detector->GetPos().Z());
    }
    for (const auto& detector : m_referenceDetectors) {
        zMin = std::min(zMin, detector->GetPos().Z());
        zMax = std::max(zMax, detector->GetPos().Z());
    }
    const double margin = std::max(50.0, 0.05 * (zMax - zMin));
    zMin -= margin;
    zMax += margin;
    constexpr int graphPoints = 301;
    TGraph resolutionX(graphPoints), resolutionY(graphPoints);
    resolutionX.SetName("gTrackingResolutionXVsZ");
    resolutionY.SetName("gTrackingResolutionYVsZ");
    resolutionX.SetTitle("Three-tracker pointing resolution;z [mm];Pointing resolution [#mum]");
    resolutionY.SetTitle("Three-tracker pointing resolution;z [mm];Pointing resolution [#mum]");
    double maximumResolution = 0.0;
    for (int point = 0; point < graphPoints; ++point) {
        const double z = zMin + (zMax - zMin) * point / (graphPoints - 1);
        const double scale = pointingScale(z);
        const double x = 1000.0 * sigmaHitX * scale;
        const double y = 1000.0 * sigmaHitY * scale;
        resolutionX.SetPoint(point, z, x);
        resolutionY.SetPoint(point, z, y);
        maximumResolution = std::max({maximumResolution, x, y});
    }

    auto* directory = m_output->GetDirectory("Resolution");
    TDirectory::TContext context(directory);
    TCanvas commonFitCanvas("cCommonHitGaussianFit", "Common single-hit Gaussian fit", 1400, 600);
    commonFitCanvas.Divide(2, 1);
    commonFitCanvas.cd(1);
    DrawGaussianFit(m_commonEquivalentHitX, commonFitX);
    commonFitCanvas.cd(2);
    DrawGaussianFit(m_commonEquivalentHitY, commonFitY);
    commonFitCanvas.Write();
    resolutionX.Write();
    resolutionY.Write();

    TCanvas curveCanvas("cTrackingResolutionVsZ", "Tracking pointing resolution vs z", 1200, 800);
    curveCanvas.SetLeftMargin(0.12);
    curveCanvas.SetRightMargin(0.05);
    curveCanvas.SetBottomMargin(0.12);
    resolutionX.SetLineColor(kBlue + 1);
    resolutionX.SetLineWidth(3);
    resolutionY.SetLineColor(kRed + 1);
    resolutionY.SetLineWidth(3);
    resolutionX.Draw("AL");
    resolutionX.GetYaxis()->SetRangeUser(0.0, 1.15 * maximumResolution);
    resolutionY.Draw("L SAME");
    TLegend legend(0.72, 0.77, 0.92, 0.90);
    legend.AddEntry(&resolutionX, "X", "l");
    legend.AddEntry(&resolutionY, "Y", "l");
    legend.Draw();
    curveCanvas.Write();

    std::ostringstream report;
    report << std::fixed << std::setprecision(2)
           << "[Resolution] hit(X,Y)=(" << 1000.0 * sigmaHitX << ", "
           << 1000.0 * sigmaHitY << ") um, angle(X,Y)=("
           << 1.0e6 * sigmaHitX / std::sqrt(szz) << ", "
           << 1.0e6 * sigmaHitY / std::sqrt(szz) << ") urad";
    for (const auto& detector : m_referenceDetectors) {
        const double scale = pointingScale(detector->GetPos().Z());
        report << ", DUT" << detector->GetID() << "@z=" << detector->GetPos().Z()
               << " mm pointing(X,Y)=(" << 1000.0 * sigmaHitX * scale << ", "
               << 1000.0 * sigmaHitY * scale << ") um";
    }
    std::cout << report.str() << '\n';
    m_output->Write();
}

}  // namespace Tracking
