#include "Script/TrackAnalysisScript.h"
#include "Detector/DetectorFactory.h"
#include "Event/DetectorFrame.h"
#include "Script/Base/RawDataParser.h"
#include "Script/Base/ScriptFactory.h"
#include "Terminal.h"

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TTree.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <sstream>
#include <utility>
#include <vector>

using namespace std;

namespace {

struct AxisTrackSamples {
    vector<double> clusterCharge;
    vector<int> clusterSize;
    vector<double> hitADC;
};

struct DetectorTrackSamples {
    int typeX = -1;
    int typeY = -1;
    AxisTrackSamples x;
    AxisTrackSamples y;
};

pair<int, int> FindPlanarAxisTypes(const Detector& detector) {
    const auto* config = detector.GetPlanarConfig();
    if (!config) return {-1, -1};

    int typeX = -1;
    int typeY = -1;
    double bestX = numeric_limits<double>::infinity();
    double bestY = numeric_limits<double>::infinity();
    constexpr double degreesToRadians = 3.14159265358979323846 / 180.0;
    for (const int type : config->readoutPlaneType) {
        const auto angleIt = config->readoutPlaneAngle.find(type);
        if (angleIt == config->readoutPlaneAngle.end()) continue;
        const double angle = angleIt->second * degreesToRadians;
        const double xScore = abs(sin(angle));
        const double yScore = abs(cos(angle));
        if (xScore < bestX) {
            bestX = xScore;
            typeX = type;
        }
        if (yScore < bestY) {
            bestY = yScore;
            typeY = type;
        }
    }
    if (typeX == typeY) return {-1, -1};
    return {typeX, typeY};
}

AxisTrackSamples* AxisSamplesForType(DetectorTrackSamples& samples, int type) {
    if (type == samples.typeX) return &samples.x;
    if (type == samples.typeY) return &samples.y;
    return nullptr;
}

pair<double, double> CommonRange(const vector<double>& x,
                                 const vector<double>& y) {
    double minimum = numeric_limits<double>::infinity();
    double maximum = -numeric_limits<double>::infinity();
    const auto update = [&](const vector<double>& values) {
        for (const double value : values) {
            if (!isfinite(value)) continue;
            minimum = min(minimum, value);
            maximum = max(maximum, value);
        }
    };
    update(x);
    update(y);
    if (!isfinite(minimum) || !isfinite(maximum)) return {0.0, 1.0};

    double low = min(0.0, minimum);
    double high = max(0.0, maximum);
    if (high <= low) return {low - 0.5, high + 0.5};
    const double padding = 0.05 * (high - low);
    if (low < 0.0) low -= padding;
    high += padding;
    return {low, high};
}

void StyleAxisHistogram(TH1D& histogram, Color_t color) {
    histogram.SetDirectory(nullptr);
    histogram.SetStats(false);
    histogram.SetLineColor(color);
    histogram.SetMarkerColor(color);
    histogram.SetLineWidth(2);
}

void WriteOverlayCanvas(TDirectory* directory, TH1D& histogramX,
                        TH1D& histogramY, const string& canvasName,
                        const string& canvasTitle) {
    if (!directory) return;
    TDirectory::TContext context(directory);

    TCanvas canvas(canvasName.c_str(), canvasTitle.c_str(), 900, 700);
    histogramX.SetMaximum(max(1.0, 1.15 * max(histogramX.GetMaximum(),
                                             histogramY.GetMaximum())));
    histogramX.Draw("HIST");
    histogramY.Draw("HIST SAME");
    TLegend legend(0.68, 0.75, 0.88, 0.88);
    legend.SetBorderSize(1);
    legend.SetFillColor(kWhite);
    legend.AddEntry(&histogramX,
                    ("X (N=" + to_string(
                         static_cast<long long>(histogramX.GetEntries())) + ")").c_str(),
                    "l");
    legend.AddEntry(&histogramY,
                    ("Y (N=" + to_string(
                         static_cast<long long>(histogramY.GetEntries())) + ")").c_str(),
                    "l");
    legend.Draw();
    canvas.Write();
}

void WriteTrackHitDistributions(
    TFile& output, const vector<shared_ptr<Detector>>& trackers,
    const map<int, DetectorTrackSamples>& allSamples) {
    auto* propertiesDirectory = output.mkdir("TrackHitProperties");
    if (!propertiesDirectory) return;

    for (const auto& detector : trackers) {
        const auto sampleIt = allSamples.find(detector->GetID());
        if (sampleIt == allSamples.end()) continue;
        const auto& samples = sampleIt->second;
        auto* detectorDirectory = propertiesDirectory->mkdir(
            ("Detector_" + to_string(detector->GetID())).c_str());
        if (!detectorDirectory) continue;

        const auto chargeRange =
            CommonRange(samples.x.clusterCharge, samples.y.clusterCharge);
        TH1D chargeX(
            "hClusterChargeX",
            ("Detector " + to_string(detector->GetID()) +
             " selected-track cluster charge;Cluster charge [ADC];Entries").c_str(),
            100, chargeRange.first, chargeRange.second);
        TH1D chargeY(
            "hClusterChargeY",
            ("Detector " + to_string(detector->GetID()) +
             " selected-track cluster charge;Cluster charge [ADC];Entries").c_str(),
            100, chargeRange.first, chargeRange.second);

        int maximumClusterSize = 1;
        for (const int size : samples.x.clusterSize)
            maximumClusterSize = max(maximumClusterSize, size);
        for (const int size : samples.y.clusterSize)
            maximumClusterSize = max(maximumClusterSize, size);
        TH1D sizeX(
            "hClusterSizeX",
            ("Detector " + to_string(detector->GetID()) +
             " selected-track cluster size;Cluster size [channels];Entries").c_str(),
            maximumClusterSize, 0.5, maximumClusterSize + 0.5);
        TH1D sizeY(
            "hClusterSizeY",
            ("Detector " + to_string(detector->GetID()) +
             " selected-track cluster size;Cluster size [channels];Entries").c_str(),
            maximumClusterSize, 0.5, maximumClusterSize + 0.5);

        const auto adcRange = CommonRange(samples.x.hitADC, samples.y.hitADC);
        TH1D adcX(
            "hHitADCX",
            ("Detector " + to_string(detector->GetID()) +
             " selected-track hit amplitude;Hit amplitude [ADC];Entries").c_str(),
            100, adcRange.first, adcRange.second);
        TH1D adcY(
            "hHitADCY",
            ("Detector " + to_string(detector->GetID()) +
             " selected-track hit amplitude;Hit amplitude [ADC];Entries").c_str(),
            100, adcRange.first, adcRange.second);

        StyleAxisHistogram(chargeX, kBlue + 1);
        StyleAxisHistogram(chargeY, kRed + 1);
        StyleAxisHistogram(sizeX, kBlue + 1);
        StyleAxisHistogram(sizeY, kRed + 1);
        StyleAxisHistogram(adcX, kBlue + 1);
        StyleAxisHistogram(adcY, kRed + 1);
        for (const double value : samples.x.clusterCharge) chargeX.Fill(value);
        for (const double value : samples.y.clusterCharge) chargeY.Fill(value);
        for (const int value : samples.x.clusterSize) sizeX.Fill(value);
        for (const int value : samples.y.clusterSize) sizeY.Fill(value);
        for (const double value : samples.x.hitADC) adcX.Fill(value);
        for (const double value : samples.y.hitADC) adcY.Fill(value);

        WriteOverlayCanvas(
            detectorDirectory, chargeX, chargeY, "cClusterChargeXY",
            "Selected-track cluster charge: X vs Y");
        WriteOverlayCanvas(
            detectorDirectory, sizeX, sizeY, "cClusterSizeXY",
            "Selected-track cluster size: X vs Y");
        WriteOverlayCanvas(
            detectorDirectory, adcX, adcY, "cHitADCXY",
            "Selected-track hit amplitude: X vs Y");
    }
}

}  // namespace

void TrackAnalysisScript::LoadConfig(const json& config) {
    m_runAlignment = config.value(
        "runAlignment",
        config.value("performAlignment", false));
    m_saveValidationData = config.value("saveValidationData", true);
    m_debug = config.value("debug", false);
    m_performanceHistograms = config.value("performanceHistograms", true);
    m_useEstimatedResolution = config.value("useEstimatedResolution", true);
    m_residualHistogramRange = config.value("residualHistogramRange", 2.0);
    m_tracking.resolutionX = config.value("resolutionX", 0.12);
    m_tracking.resolutionY = config.value("resolutionY", 0.12);
    m_tracking.gateSigma = config.value("gateSigma", 3.0);
    m_tracking.maxChi2Ndf = config.value("maxChi2Ndf", 25.0);
    m_tracking.maxBranchesPerLayer = config.value("maxBranchesPerLayer", 3);
    m_tracking.maxCandidates = config.value("maxCandidates", 4000);
    m_tracking.maxTracks = config.value("maxTracksPerEvent", 32);
    m_tracking.conflictSearchNodes = config.value("conflictSearchNodes", 20000);
    m_alignment.maxIterations = config.value("alignmentIterations", 30);
    m_alignment.minTracks = config.value("alignmentMinTracks", 20);
    m_alignment.maxFunctionCalls = config.value("alignmentMaxFunctionCalls", 600);
    m_alignment.maxShiftStep = config.value("alignmentMaxShiftStep", 0.20);
    m_alignment.maxRotationStep = config.value("alignmentMaxRotationStep", 0.001);
    m_alignment.relativeLossTolerance = std::max(
        0.0, config.value("alignmentRelativeLossTolerance", 1e-4));
    m_alignment.convergencePatience = std::max(
        1, config.value("alignmentConvergencePatience", 3));
    m_alignment.debug = m_debug;
}

void TrackAnalysisScript::Print() const {
    if (Terminal::Verbose()) {
        Terminal::Detail(
            "single-hit calibration · multi-hit reconstruction · alignment " +
            string(m_runAlignment ? "on" : "off") +
            " · estimated resolution " +
            string(m_useEstimatedResolution ? "on" : "off"));
    }
}

bool TrackAnalysisScript::Execute() {
    auto parser = GetParser();
    if (!parser) {
        cerr << "[TrackAnalysis] parser not set\n";
        return false;
    }
    auto trackers = DetectorFactory::GetInstance().GetDetectorsByRole(Detector::Role::Tracker);
    sort(trackers.begin(), trackers.end(), [](const auto& a, const auto& b) { return a->GetPos().Z() < b->GetPos().Z(); });
    if (trackers.size() < 3) {
        cerr << "[TrackAnalysis] at least three trackers are required\n";
        return false;
    }
    Print();

    const Long64_t total = parser->GetTotalEvents();
    auto loadEvent = [&](Long64_t i) -> optional<Event> {
        const auto rawHits = parser->LoadEvent(i);
        if (rawHits.empty()) return nullopt;
        Event event{.eventID = static_cast<int>(i)};
        for (const auto& detector : trackers) {
            auto frame = make_shared<DetectorFrame>(*detector);
            const auto raw = rawHits.find(detector->GetID());
            if (raw != rawHits.end()) frame->SetRawData(raw->second);
            frame->Process();
            event.detectorFramesMap[detector->GetID()] = std::move(frame);
        }
        return event;
    };

    vector<Event> calibrationEvents;
    calibrationEvents.reserve(total);
    for (Long64_t i = 0; i < total; ++i) {
        auto event = loadEvent(i);
        if (!event) continue;
        bool singleHitOnEveryTracker = true;
        for (const auto& detector : trackers) {
            if (event->detectorFramesMap.at(detector->GetID())->LocalHits().size() != 1) {
                singleHitOnEveryTracker = false;
                break;
            }
        }
        if (!singleHitOnEveryTracker) continue;
        calibrationEvents.push_back(std::move(*event));
    }
    if (Terminal::Verbose()) {
        Terminal::Detail("calibration " +
                         Terminal::Count(calibrationEvents.size()) + "/" +
                         Terminal::Count(total) + " events");
    }
    if (calibrationEvents.empty()) return false;

    Tracking::Reconstructor calibrationReconstructor(trackers, m_tracking);
    if (m_runAlignment) {
        Tracking::Aligner aligner(trackers, calibrationReconstructor, m_alignment);
        if (!aligner.Run(calibrationEvents)) cerr << "[TrackAnalysis] alignment did not converge; using last valid geometry\n";
        ostringstream geometry;
        geometry << defaultfloat << setprecision(10)
                 << "[TrackAlign] final geometry (config.json format)\n";
        for (const auto& detector : trackers) {
            const auto pos = detector->GetPos();
            const auto rot = detector->GetRot();
            geometry << "  detector " << detector->GetID() << " ("
                     << detector->GetName() << "):\n"
                     << "      \"position\": [ " << pos.X() << ", " << pos.Y()
                     << ", " << pos.Z() << " ],\n"
                     << "      \"rotation\": [ " << rot.X() << ", " << rot.Y()
                     << ", " << rot.Z() << " ],\n";
        }
        cout << geometry.str();
    }

    const string outputPath = GetOutputDir() + "TrackInfo.root";
    auto output = unique_ptr<TFile>(TFile::Open(outputPath.c_str(), "RECREATE"));
    if (!output || output->IsZombie()) return false;
    TTree tracksTree("Tracks", "Track info");
    Int_t eventID = 0;
    ULong64_t rawEventID = 0;
    Int_t trackIndex = 0;
    Track track{};
    tracksTree.Branch("eventID", &eventID);
    tracksTree.Branch("rawEventID", &rawEventID);
    tracksTree.Branch("trackIndex", &trackIndex);
    tracksTree.Branch("track", &track);

    unique_ptr<TTree> validation;
    Int_t detID = 0;
    Double_t resX = 0, resY = 0, hitX = 0, hitY = 0;
    vector<Int_t> clusterIndices;
    vector<ChannelHit> channelHits;
    vector<Cluster> clusters;
    if (m_saveValidationData) {
        validation = make_unique<TTree>("TrackerValidation", "Tracker QA");
        validation->Branch("eventID", &eventID);
        validation->Branch("rawEventID", &rawEventID);
        validation->Branch("trackIndex", &trackIndex);
        validation->Branch("detID", &detID);
        validation->Branch("resX", &resX);
        validation->Branch("resY", &resY);
        validation->Branch("hitX", &hitX);
        validation->Branch("hitY", &hitY);
        validation->Branch("clusterIndices", &clusterIndices);
        validation->Branch("channelHits", &channelHits);
        validation->Branch("clusters", &clusters);
    }

    unique_ptr<Tracking::PerformanceAnalyzer> performance;
    if (m_performanceHistograms || m_useEstimatedResolution) {
        auto* performanceDir = output->mkdir("Performance");
        const auto referenceDetectors =
            DetectorFactory::GetInstance().GetDetectorsByRole(Detector::Role::DUT);
        performance = make_unique<Tracking::PerformanceAnalyzer>(
            performanceDir, trackers, referenceDetectors, m_tracking,
            m_residualHistogramRange);
    }

    if (performance && m_useEstimatedResolution) {
        for (const auto& event : calibrationEvents) performance->RecordEvent(event);
        const auto [estimatedX, estimatedY] = performance->EstimateHitResolution();
        if (estimatedX > 0.0 && estimatedY > 0.0) {
            m_tracking.resolutionX = estimatedX;
            m_tracking.resolutionY = estimatedY;
            if (Terminal::Verbose()) {
                ostringstream estimate;
                estimate << fixed << setprecision(1)
                         << "estimated hit resolution "
                         << 1000.0 * estimatedX << " × "
                         << 1000.0 * estimatedY << " µm";
                Terminal::Detail(estimate.str());
            }
        } else {
            cerr << "[TrackAnalysis] hit-resolution fit failed; using configured resolution\n";
        }
        if (m_performanceHistograms) performance->Reset();
    }

    // Calibration retains full detector frames, including raw waveforms. Free
    // them before the all-event reconstruction pass to keep memory bounded.
    calibrationEvents.clear();
    calibrationEvents.shrink_to_fit();
    Tracking::Reconstructor reconstructor(trackers, m_tracking);

    map<int, DetectorTrackSamples> trackHitSamples;
    for (const auto& detector : trackers) {
        const auto [typeX, typeY] = FindPlanarAxisTypes(*detector);
        if (typeX < 0 || typeY < 0) {
            cerr << "[TrackAnalysis] detector " << detector->GetID()
                 << " has no distinct planar X/Y readout; "
                    "track-hit distributions will be skipped\n";
            continue;
        }
        trackHitSamples.emplace(
            detector->GetID(),
            DetectorTrackSamples{.typeX = typeX, .typeY = typeY});
    }

    size_t savedTracks = 0, trackedEvents = 0, analyzedEvents = 0, multiHitEvents = 0;
    for (Long64_t i = 0; i < total; ++i) {
        auto loaded = loadEvent(i);
        if (!loaded) continue;
        auto& event = *loaded;
        ++analyzedEvents;
        const bool hasMultipleHits = any_of(
            trackers.begin(), trackers.end(), [&](const auto& detector) {
                return event.detectorFramesMap.at(detector->GetID())->LocalHits().size() > 1;
            });
        if (hasMultipleHits) ++multiHitEvents;
        const ULong64_t currentRawEventID = parser->GetCurrentEventID();
        auto results = reconstructor.Reconstruct(event);
        if (!results.empty()) ++trackedEvents;
        for (size_t resultIndex = 0; resultIndex < results.size(); ++resultIndex) {
            const auto& result = results[resultIndex];
            eventID = event.eventID;
            rawEventID = currentRawEventID;
            trackIndex = static_cast<Int_t>(resultIndex);
            track = result.track;
            tracksTree.Fill();
            ++savedTracks;
            if (performance && m_performanceHistograms)
                performance->RecordTrack(event, result);

            for (const auto& [id, hitIndex] : result.hitIndices) {
                const auto sampleIt = trackHitSamples.find(id);
                if (sampleIt == trackHitSamples.end()) continue;
                const auto& frame = event.detectorFramesMap.at(id);
                const auto& localHit = frame->LocalHits().at(hitIndex);
                const auto& frameClusters = frame->Clusters();
                const auto& frameChannelHits = frame->ChannelHits();
                for (const int clusterIndex : localHit.clusterIndices) {
                    if (clusterIndex < 0 ||
                        clusterIndex >= static_cast<int>(frameClusters.size()))
                        continue;
                    const auto& selectedCluster = frameClusters[clusterIndex];
                    auto* axis =
                        AxisSamplesForType(sampleIt->second, selectedCluster.type);
                    if (!axis) continue;
                    if (isfinite(selectedCluster.charge))
                        axis->clusterCharge.push_back(selectedCluster.charge);
                    axis->clusterSize.push_back(selectedCluster.size);
                    for (const int channelHitIndex :
                         selectedCluster.channelHitIndices) {
                        if (channelHitIndex < 0 ||
                            channelHitIndex >=
                                static_cast<int>(frameChannelHits.size()))
                            continue;
                        const auto& selectedHit =
                            frameChannelHits[channelHitIndex];
                        if (selectedHit.isValid && isfinite(selectedHit.amp))
                            axis->hitADC.push_back(selectedHit.amp);
                    }
                }
            }

            if (!validation) continue;
            for (const auto& [id, hitIndex] : result.hitIndices) {
                const auto detector = DetectorFactory::GetInstance().GetDetector(id);
                const auto& frame = event.detectorFramesMap.at(id);
                const auto& localHit = frame->LocalHits().at(hitIndex);
                vector<TVector3> otherHits;
                for (const auto& [otherID, otherIndex] : result.hitIndices) {
                    if (otherID == id) continue;
                    const auto other = DetectorFactory::GetInstance().GetDetector(otherID);
                    otherHits.push_back(other->LocalToGlobal(event.detectorFramesMap.at(otherID)->LocalHits().at(otherIndex).localPos));
                }
                const auto unbiased = Tracking::FitWeighted(otherHits, m_tracking.resolutionX, m_tracking.resolutionY);
                const auto predicted = detector->GlobalToLocal(detector->CalcHitFromTrack(unbiased));
                detID = id;
                hitX = localHit.localPos.X();
                hitY = localHit.localPos.Y();
                resX = hitX - predicted.X();
                resY = hitY - predicted.Y();
                clusterIndices.assign(localHit.clusterIndices.begin(), localHit.clusterIndices.end());
                channelHits = frame->ChannelHits();
                clusters = frame->Clusters();
                validation->Fill();
            }
        }
    }
    ostringstream trackingSummary;
    const double trackedFraction =
        analyzedEvents > 0
            ? 100.0 * static_cast<double>(trackedEvents) / analyzedEvents
            : 0.0;
    trackingSummary << Terminal::Count(trackedEvents) << '/'
                    << Terminal::Count(analyzedEvents) << " events ("
                    << fixed << setprecision(1) << trackedFraction
                    << "%) · " << Terminal::Count(savedTracks) << " tracks";
    Terminal::Detail(trackingSummary.str());
    output->cd();
    tracksTree.Write();
    if (validation) validation->Write();
    if (performance && m_performanceHistograms) performance->Write();
    WriteTrackHitDistributions(*output, trackers, trackHitSamples);
    if (Terminal::Verbose()) {
        for (const auto& detector : trackers) {
            const auto sampleIt = trackHitSamples.find(detector->GetID());
            if (sampleIt == trackHitSamples.end()) continue;
            const auto& samples = sampleIt->second;
            ostringstream detail;
            detail << "tracker " << detector->GetID()
                   << " clusters X/Y "
                   << Terminal::Count(samples.x.clusterSize.size()) << '/'
                   << Terminal::Count(samples.y.clusterSize.size())
                   << " · hits X/Y "
                   << Terminal::Count(samples.x.hitADC.size()) << '/'
                   << Terminal::Count(samples.y.hitADC.size());
            Terminal::Detail(Terminal::Muted(detail.str()));
        }
    }
    tracksTree.SetDirectory(nullptr);
    if (validation) validation->SetDirectory(nullptr);
    output->Close();
    return true;
}

REGISTER_SCRIPT("TrackAnalysis", TrackAnalysisScript);
