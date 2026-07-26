#include "Script/TrackAnalysisScript.h"
#include "Detector/DetectorFactory.h"
#include "Event/DetectorFrame.h"
#include "Script/Base/RawDataParser.h"
#include "Script/Base/ScriptFactory.h"

#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <optional>
#include <sstream>

using namespace std;

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
    cout << "[Tracking] single-hit calibration, multi-hit reconstruction, alignment="
         << (m_runAlignment ? "on" : "off")
         << ", estimated-resolution=" << (m_useEstimatedResolution ? "on" : "off") << '\n';
}

bool TrackAnalysisScript::Execute() {
    const auto start = chrono::steady_clock::now();
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
    cout << "[Tracking] input=" << parser->GetTotalEvents()
         << " events, trackers=" << trackers.size() << '\n';
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
    cout << "[Tracking] calibration events=" << calibrationEvents.size() << '/' << total << '\n';
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
            cout << "[Tracking] estimated hit resolution X=" << 1000.0 * estimatedX
                 << " um, Y=" << 1000.0 * estimatedY << " um\n";
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
    cout << "[Tracking] reconstructed=" << trackedEvents << '/' << analyzedEvents
         << " events, multi-hit input=" << multiHitEvents
         << ", tracks=" << savedTracks << '\n';
    output->cd();
    tracksTree.Write();
    if (validation) validation->Write();
    if (performance && m_performanceHistograms) performance->Write();
    tracksTree.SetDirectory(nullptr);
    if (validation) validation->SetDirectory(nullptr);
    output->Close();
    cout << "[Tracking] output=TrackInfo.root, elapsed="
         << fixed << setprecision(2) << chrono::duration<double>(chrono::steady_clock::now() - start).count() << " s\n";
    return true;
}

REGISTER_SCRIPT("TrackAnalysis", TrackAnalysisScript);
