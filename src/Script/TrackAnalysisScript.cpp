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

using namespace std;

void TrackAnalysisScript::LoadConfig(const json& config) {
    m_saveValidationData = config.value("saveValidationData", true);
    m_progressInterval = max(1, config.value("progressInterval", 10000));
    m_debug = config.value("debug", true);
    m_performanceHistograms = config.value("performanceHistograms", true);
    m_residualHistogramRange = config.value("residualHistogramRange", 2.0);
    m_tracking.resolutionX = config.value("resolutionX", 0.12);
    m_tracking.resolutionY = config.value("resolutionY", 0.12);
    m_tracking.gateSigma = config.value("gateSigma", 5.0);
    m_tracking.maxChi2Ndf = config.value("maxChi2Ndf", 25.0);
    m_tracking.maxBranchesPerLayer = config.value("maxBranchesPerLayer", 3);
    m_tracking.maxCandidates = config.value("maxCandidates", 4000);
    m_tracking.maxTracks = config.value("maxTracksPerEvent", 32);
    m_tracking.conflictSearchNodes = config.value("conflictSearchNodes", 20000);
    m_alignment.maxIterations = config.value("alignmentIterations", 20);
    m_alignment.minTracks = config.value("alignmentMinTracks", 20);
    m_alignment.maxEvents = config.value("alignmentMaxEvents", 10000);
    m_alignment.maxTracks = config.value("alignmentMaxTracks", 5000);
    m_alignment.maxFunctionCalls = config.value("alignmentMaxFunctionCalls", 600);
    m_alignment.maxShiftStep = config.value("alignmentMaxShiftStep", 0.20);
    m_alignment.maxRotationStep = config.value("alignmentMaxRotationStep", 0.001);
    m_alignment.debug = m_debug;
}

void TrackAnalysisScript::Print() const {
    cout << "TrackAnalysisScript: gate=" << m_tracking.gateSigma << " sigma, resolution=("
         << m_tracking.resolutionX << ", " << m_tracking.resolutionY << ") mm, maxCandidates="
         << m_tracking.maxCandidates << ", maxTracks/event=" << m_tracking.maxTracks << '\n';
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
    cout << "[TrackAnalysis] trackers=" << trackers.size() << " events=" << parser->GetTotalEvents() << '\n';
    Print();

    const Long64_t total = parser->GetTotalEvents();
    vector<Event> events;
    map<int, ULong64_t> rawEventIDs;
    events.reserve(total);
    const int requiredLayers = trackers.size() <= 3 ? trackers.size() : trackers.size() - 1;
    for (Long64_t i = 0; i < total; ++i) {
        if (i % m_progressInterval == 0 || i + 1 == total)
            cout << "\r[TrackAnalysis] decode " << i + 1 << '/' << total << flush;
        const auto rawHits = parser->LoadEvent(i);
        if (rawHits.empty()) continue;
        Event event{.eventID = static_cast<int>(i)};
        rawEventIDs[event.eventID] = parser->GetCurrentEventID();
        int validLayers = 0;
        for (const auto& detector : trackers) {
            auto frame = make_shared<DetectorFrame>(*detector);
            const auto raw = rawHits.find(detector->GetID());
            if (raw != rawHits.end()) frame->SetRawData(raw->second);
            if (frame->Process()) ++validLayers;
            event.detectorFramesMap[detector->GetID()] = std::move(frame);
        }
        if (validLayers >= requiredLayers) events.push_back(std::move(event));
    }
    cout << "\n[TrackAnalysis] reconstructable events=" << events.size() << '/' << total << '\n';
    if (events.empty()) return false;

    Tracking::Reconstructor reconstructor(trackers, m_tracking);
    vector<Tracking::AlignmentIteration> alignmentHistory;
    cout << "[TrackAnalysis] perform tracker alignment? (y/n): " << flush;
    char alignChoice = 'n';
    cin >> alignChoice;
    if (alignChoice == 'y' || alignChoice == 'Y') {
        Tracking::Aligner aligner(trackers, reconstructor, m_alignment);
        if (!aligner.Run(events)) cerr << "[TrackAnalysis] alignment did not converge; using last valid geometry\n";
        alignmentHistory = aligner.History();
        cout << "[TrackAlign] final geometry\n";
        for (const auto& detector : trackers) {
            const auto pos = detector->GetPos();
            const auto rot = detector->GetRot();
            cout << "  id=" << detector->GetID() << " pos=(" << pos.X() << ',' << pos.Y() << ',' << pos.Z()
                 << ") rotZ=" << rot.Z() << '\n';
        }
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
    vector<StripHit> stripHits;
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
        validation->Branch("stripHits", &stripHits);
        validation->Branch("clusters", &clusters);
    }

    unique_ptr<Tracking::PerformanceAnalyzer> performance;
    if (m_performanceHistograms) {
        auto* performanceDir = output->mkdir("Performance");
        performance = make_unique<Tracking::PerformanceAnalyzer>(
            performanceDir, trackers, m_tracking, total, m_residualHistogramRange);
        performance->RecordAlignment(alignmentHistory);
    }

    size_t savedTracks = 0, trackedEvents = 0;
    for (size_t i = 0; i < events.size(); ++i) {
        Tracking::ReconstructionStats reconstructionStats;
        auto results = reconstructor.Reconstruct(events[i], &reconstructionStats);
        if (performance) performance->RecordEvent(events[i], results, reconstructionStats);
        if (!results.empty()) ++trackedEvents;
        for (size_t resultIndex = 0; resultIndex < results.size(); ++resultIndex) {
            const auto& result = results[resultIndex];
            eventID = events[i].eventID;
            rawEventID = rawEventIDs.at(events[i].eventID);
            trackIndex = static_cast<Int_t>(resultIndex);
            track = result.track;
            tracksTree.Fill();
            ++savedTracks;

            if (!validation) continue;
            for (const auto& [id, hitIndex] : result.hitIndices) {
                const auto detector = DetectorFactory::GetInstance().GetDetector(id);
                const auto& frame = events[i].detectorFramesMap.at(id);
                const auto& localHit = frame->LocalHits().at(hitIndex);
                vector<TVector3> otherHits;
                for (const auto& [otherID, otherIndex] : result.hitIndices) {
                    if (otherID == id) continue;
                    const auto other = DetectorFactory::GetInstance().GetDetector(otherID);
                    otherHits.push_back(other->LocalToGlobal(events[i].detectorFramesMap.at(otherID)->LocalHits().at(otherIndex).localPos));
                }
                const auto unbiased = Tracking::FitWeighted(otherHits, m_tracking.resolutionX, m_tracking.resolutionY);
                const auto predicted = detector->GlobalToLocal(detector->CalcHitFromTrack(unbiased));
                detID = id;
                hitX = localHit.localPos.X();
                hitY = localHit.localPos.Y();
                resX = hitX - predicted.X();
                resY = hitY - predicted.Y();
                clusterIndices.assign(localHit.clusterIndices.begin(), localHit.clusterIndices.end());
                stripHits = frame->StripHits();
                clusters = frame->Clusters();
                validation->Fill();
            }
        }
        if (m_debug && ((i + 1) % m_progressInterval == 0 || i + 1 == events.size()))
            cout << "\r[TrackAnalysis] fit " << i + 1 << '/' << events.size() << " tracks=" << savedTracks << flush;
    }
    cout << "\n[TrackAnalysis] tracked events=" << trackedEvents << '/' << events.size()
         << " tracks=" << savedTracks << '\n';
    output->cd();
    tracksTree.Write();
    if (validation) validation->Write();
    if (performance) performance->Write();
    cout << "[TrackAnalysis] output=" << outputPath << " time="
         << fixed << setprecision(2) << chrono::duration<double>(chrono::steady_clock::now() - start).count() << " s\n";
    return true;
}

REGISTER_SCRIPT("TrackAnalysis", TrackAnalysisScript);
