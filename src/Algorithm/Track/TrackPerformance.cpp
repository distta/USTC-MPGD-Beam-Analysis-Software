#include "Algorithm/TrackPerformance.h"
#include "Detector/Detector.h"
#include "Event/DetectorFrame.h"

#include <TDirectory.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>

#include <algorithm>
#include <cmath>
#include <string>

namespace Tracking {
namespace {

double Sigma68(std::vector<double> values) {
    if (values.empty()) return 0.0;
    const auto middle = values.begin() + values.size() / 2;
    std::nth_element(values.begin(), middle, values.end());
    const double median = *middle;
    for (double& value : values) value = std::abs(value - median);
    const size_t index = std::min(values.size() - 1, static_cast<size_t>(0.6827 * values.size()));
    std::nth_element(values.begin(), values.begin() + index, values.end());
    return values[index];
}

void FillSummary(TH1D* summary, TH1D* x, TH1D* y,
                 const std::vector<double>& valuesX, const std::vector<double>& valuesY) {
    const double sx68 = Sigma68(valuesX), sy68 = Sigma68(valuesY);
    summary->SetBinContent(1, x->GetMean());
    summary->SetBinContent(2, x->GetRMS());
    summary->SetBinContent(3, sx68);
    summary->SetBinContent(5, y->GetMean());
    summary->SetBinContent(6, y->GetRMS());
    summary->SetBinContent(7, sy68);
    if (x->GetEntries() > 20 && sx68 > 0) {
        x->Fit("gaus", "Q0", "", x->GetMean() - 2 * sx68, x->GetMean() + 2 * sx68);
        if (x->GetFunction("gaus")) summary->SetBinContent(4, std::abs(x->GetFunction("gaus")->GetParameter(2)));
    }
    if (y->GetEntries() > 20 && sy68 > 0) {
        y->Fit("gaus", "Q0", "", y->GetMean() - 2 * sy68, y->GetMean() + 2 * sy68);
        if (y->GetFunction("gaus")) summary->SetBinContent(8, std::abs(y->GetFunction("gaus")->GetParameter(2)));
    }
}

}  // namespace

PerformanceAnalyzer::PerformanceAnalyzer(TDirectory* output,
                                         std::vector<std::shared_ptr<Detector>> detectors, const Config& config,
                                         int totalEvents, double residualRange)
    : m_output(output), m_detectors(std::move(detectors)), m_config(config), m_totalEvents(totalEvents), m_residualRange(residualRange) {
    auto* global = m_output->mkdir("Global");
    TDirectory::TContext globalContext(global);
    m_tracksPerEvent = new TH1D("hTracksPerEvent", "Reconstructed tracks per event;Tracks;Events", config.maxTracks + 1, -0.5, config.maxTracks + 0.5);
    m_seedCandidates = new TH1D("hSeedCandidates", "Seed candidates per event;Candidates;Events", 101, -0.5, config.maxCandidates + 0.5);
    m_finalCandidates = new TH1D("hFinalCandidates", "Quality-selected candidates per event;Candidates;Events", 101, -0.5, 100.5);
    m_conflictNodes = new TH1D("hConflictSearchNodes", "Conflict-search nodes per event;Visited nodes;Events", 100, 0, config.conflictSearchNodes);
    m_chi2 = new TH1D("hTrackChi2Ndf", "Track fit quality;#chi^{2}/ndf;Tracks", 200, 0, config.maxChi2Ndf);
    m_layersPerTrack = new TH1D("hLayersPerTrack", "Layers used per track;Layers;Tracks", detectors.size() + 1, -0.5, detectors.size() + 0.5);
    m_kx = new TH1D("hTrackSlopeX", "Track slope X;k_{x};Tracks", 200, -0.1, 0.1);
    m_ky = new TH1D("hTrackSlopeY", "Track slope Y;k_{y};Tracks", 200, -0.1, 0.1);
    m_bx = new TH1D("hTrackInterceptX", "Track intercept X;b_{x} [mm];Tracks", 256, -10, 120);
    m_by = new TH1D("hTrackInterceptY", "Track intercept Y;b_{y} [mm];Tracks", 256, -10, 120);
    m_angle = new TH1D("hTrackAngle", "Track polar angle;atan(#sqrt{k_{x}^{2}+k_{y}^{2}}) [rad];Tracks", 200, 0, 0.15);

    for (const auto& detector : m_detectors) {
        const int id = detector->GetID();
        auto* detectorDir = m_output->mkdir(("Detector_" + std::to_string(id)).c_str());
        auto* hitsDir = detectorDir->mkdir("Hits");
        auto* residualDir = detectorDir->mkdir("Residuals");
        double xMax = 1.0, yMax = 1.0;
        const auto& detConfig = detector->getConfig();
        if (detConfig.readoutPlaneStripNumber.count(0)) xMax = detConfig.readoutPlaneStripNumber.at(0) * detConfig.readoutPlanePitch.at(0);
        if (detConfig.readoutPlaneStripNumber.count(1)) yMax = detConfig.readoutPlaneStripNumber.at(1) * detConfig.readoutPlanePitch.at(1);
        auto& h = m_detectorHistograms[id];
        {
            TDirectory::TContext context(hitsDir);
            h.hitMultiplicity = new TH1D("hHitMultiplicity", "Local-hit multiplicity;Local hits/event;Events", 65, -0.5, 64.5);
            h.clusterMultiplicity = new TH1D("hClusterMultiplicity", "Cluster multiplicity;Clusters/event;Events", 101, -0.5, 100.5);
            h.clusterSize = new TH1D("hClusterSize", "Cluster size;Strips/cluster;Clusters", 21, -0.5, 20.5);
            h.allHitX = new TH1D("hAllHitX", "All local hits X;Local X [mm];Hits", 256, 0, xMax);
            h.allHitY = new TH1D("hAllHitY", "All local hits Y;Local Y [mm];Hits", 256, 0, yMax);
            h.selectedHitX = new TH1D("hSelectedHitX", "Track-associated hits X;Local X [mm];Hits", 256, 0, xMax);
            h.selectedHitY = new TH1D("hSelectedHitY", "Track-associated hits Y;Local Y [mm];Hits", 256, 0, yMax);
            h.allHitXY = new TH2D("hAllHitXY", "All local hits;Local X [mm];Local Y [mm]", 128, 0, xMax, 128, 0, yMax);
            h.selectedHitXY = new TH2D("hSelectedHitXY", "Track-associated hits;Local X [mm];Local Y [mm]", 128, 0, xMax, 128, 0, yMax);
        }
        {
            TDirectory::TContext context(residualDir);
            h.residualX = new TH1D("hUnbiasedResidualX", "Unbiased residual X;X_{meas}-X_{pred} [mm];Entries", 240, -residualRange, residualRange);
            h.residualY = new TH1D("hUnbiasedResidualY", "Unbiased residual Y;Y_{meas}-Y_{pred} [mm];Entries", 240, -residualRange, residualRange);
            h.residualXY = new TH2D("hUnbiasedResidualXY", "Unbiased residuals;Residual X [mm];Residual Y [mm]", 120, -residualRange, residualRange, 120, -residualRange, residualRange);
            h.pullX = new TH1D("hPullX", "Unbiased pull X;Residual X/#sigma_{X};Entries", 200, -10, 10);
            h.pullY = new TH1D("hPullY", "Unbiased pull Y;Residual Y/#sigma_{Y};Entries", 200, -10, 10);
            h.residualXVsHitX = new TH2D("hResidualXVsHitX", "Residual X vs hit X;Local X [mm];Residual X [mm]", 128, 0, xMax, 120, -residualRange, residualRange);
            h.residualYVsHitY = new TH2D("hResidualYVsHitY", "Residual Y vs hit Y;Local Y [mm];Residual Y [mm]", 128, 0, yMax, 120, -residualRange, residualRange);
            h.residualXVsSlope = new TH2D("hResidualXVsSlopeX", "Residual X vs slope X;k_{x};Residual X [mm]", 100, -0.1, 0.1, 120, -residualRange, residualRange);
            h.residualYVsSlope = new TH2D("hResidualYVsSlopeY", "Residual Y vs slope Y;k_{y};Residual Y [mm]", 100, -0.1, 0.1, 120, -residualRange, residualRange);
            h.residualXVsEvent = new TH2D("hResidualXVsEvent", "Residual X vs event;Event ID;Residual X [mm]", 200, 0, totalEvents, 120, -residualRange, residualRange);
            h.residualYVsEvent = new TH2D("hResidualYVsEvent", "Residual Y vs event;Event ID;Residual Y [mm]", 200, 0, totalEvents, 120, -residualRange, residualRange);
            h.summary = new TH1D("hResidualSummary", "Residual summary;Metric;Value [mm]", 8, 0.5, 8.5);
            const char* labels[] = {"meanX", "RMSX", "sigma68X", "gausSigmaX", "meanY", "RMSY", "sigma68Y", "gausSigmaY"};
            for (int bin = 1; bin <= 8; ++bin) h.summary->GetXaxis()->SetBinLabel(bin, labels[bin - 1]);
        }
    }
}

void PerformanceAnalyzer::RecordEvent(const Event& event, const std::vector<Result>& tracks,
                                      const ReconstructionStats& stats) {
    m_tracksPerEvent->Fill(tracks.size());
    m_seedCandidates->Fill(stats.seedCandidates);
    m_finalCandidates->Fill(stats.finalCandidates);
    m_conflictNodes->Fill(stats.conflictSearchNodes);
    for (const auto& detector : m_detectors) {
        auto& h = m_detectorHistograms[detector->GetID()];
        const auto& frame = event.detectorFramesMap.at(detector->GetID());
        h.hitMultiplicity->Fill(frame->LocalHits().size());
        h.clusterMultiplicity->Fill(frame->Clusters().size());
        for (const auto& cluster : frame->Clusters()) h.clusterSize->Fill(cluster.size);
        for (const auto& hit : frame->LocalHits()) {
            h.allHitX->Fill(hit.localPos.X());
            h.allHitY->Fill(hit.localPos.Y());
            h.allHitXY->Fill(hit.localPos.X(), hit.localPos.Y());
        }
    }
    for (const auto& result : tracks) {
        m_chi2->Fill(result.track.chi2);
        m_layersPerTrack->Fill(result.hitIndices.size());
        m_kx->Fill(result.track.kx);
        m_ky->Fill(result.track.ky);
        m_bx->Fill(result.track.bx);
        m_by->Fill(result.track.by);
        m_angle->Fill(std::atan(std::hypot(result.track.kx, result.track.ky)));
        for (const auto& [id, hitIndex] : result.hitIndices) {
            auto& h = m_detectorHistograms[id];
            const auto& hit = event.detectorFramesMap.at(id)->LocalHits().at(hitIndex);
            h.selectedHitX->Fill(hit.localPos.X());
            h.selectedHitY->Fill(hit.localPos.Y());
            h.selectedHitXY->Fill(hit.localPos.X(), hit.localPos.Y());
            std::vector<TVector3> otherHits;
            for (const auto& [otherID, otherIndex] : result.hitIndices) {
                if (otherID == id) continue;
                const auto detector = std::find_if(m_detectors.begin(), m_detectors.end(), [otherID](const auto& d) { return d->GetID() == otherID; });
                otherHits.push_back((*detector)->LocalToGlobal(event.detectorFramesMap.at(otherID)->LocalHits().at(otherIndex).localPos));
            }
            if (otherHits.size() < 2) continue;
            const auto fit = FitWeighted(otherHits, m_config.resolutionX, m_config.resolutionY);
            const auto detector = std::find_if(m_detectors.begin(), m_detectors.end(), [id](const auto& d) { return d->GetID() == id; });
            const auto predicted = (*detector)->GlobalToLocal((*detector)->CalcHitFromTrack(fit));
            const double rx = hit.localPos.X() - predicted.X(), ry = hit.localPos.Y() - predicted.Y();
            h.residualX->Fill(rx);
            h.residualY->Fill(ry);
            h.residualXY->Fill(rx, ry);
            h.pullX->Fill(rx / m_config.resolutionX);
            h.pullY->Fill(ry / m_config.resolutionY);
            h.residualXVsHitX->Fill(hit.localPos.X(), rx);
            h.residualYVsHitY->Fill(hit.localPos.Y(), ry);
            h.residualXVsSlope->Fill(result.track.kx, rx);
            h.residualYVsSlope->Fill(result.track.ky, ry);
            h.residualXVsEvent->Fill(event.eventID, rx);
            h.residualYVsEvent->Fill(event.eventID, ry);
            h.residualValuesX.push_back(rx);
            h.residualValuesY.push_back(ry);
        }
    }
}

void PerformanceAnalyzer::RecordAlignment(const std::vector<AlignmentIteration>& history) {
    if (history.empty()) return;
    auto* directory = m_output->mkdir("Alignment");
    TDirectory::TContext context(directory);
    const int n = history.size();
    auto* before = new TH1D("hLossBefore", "Alignment loss before update;Iteration;Robust loss", n, 0.5, n + 0.5);
    auto* after = new TH1D("hLossAfter", "Alignment loss after update;Iteration;Robust loss", n, 0.5, n + 0.5);
    auto* tracks = new TH1D("hTracksUsed", "Tracks used by alignment;Iteration;Tracks", n, 0.5, n + 0.5);
    auto* shift = new TH1D("hMaxShiftStep", "Maximum shift update;Iteration;Shift [mm]", n, 0.5, n + 0.5);
    auto* rotation = new TH1D("hMaxRotationStep", "Maximum rotation update;Iteration;Rotation [rad]", n, 0.5, n + 0.5);
    std::map<int, TH1D*> dx, dy, rz;
    for (const auto& detector : m_detectors) {
        const auto suffix = std::to_string(detector->GetID());
        dx[detector->GetID()] = new TH1D(("hAlignDx_" + suffix).c_str(), ("Detector " + suffix + " alignment dx;Iteration;dx [mm]").c_str(), n, 0.5, n + 0.5);
        dy[detector->GetID()] = new TH1D(("hAlignDy_" + suffix).c_str(), ("Detector " + suffix + " alignment dy;Iteration;dy [mm]").c_str(), n, 0.5, n + 0.5);
        rz[detector->GetID()] = new TH1D(("hAlignRz_" + suffix).c_str(), ("Detector " + suffix + " alignment rz;Iteration;rz [rad]").c_str(), n, 0.5, n + 0.5);
    }
    for (int i = 0; i < n; ++i) {
        before->SetBinContent(i + 1, history[i].lossBefore);
        after->SetBinContent(i + 1, history[i].lossAfter);
        tracks->SetBinContent(i + 1, history[i].tracks);
        shift->SetBinContent(i + 1, history[i].maxShift);
        rotation->SetBinContent(i + 1, history[i].maxRotation);
        for (const auto& detector : m_detectors) {
            const int id = detector->GetID();
            dx[id]->SetBinContent(i + 1, history[i].alignPosition.at(id).X());
            dy[id]->SetBinContent(i + 1, history[i].alignPosition.at(id).Y());
            rz[id]->SetBinContent(i + 1, history[i].alignRotation.at(id).Z());
        }
    }
}

void PerformanceAnalyzer::Write() {
    for (auto& [id, h] : m_detectorHistograms)
        FillSummary(h.summary, h.residualX, h.residualY, h.residualValuesX, h.residualValuesY);
    m_output->Write();
}

}  // namespace Tracking
