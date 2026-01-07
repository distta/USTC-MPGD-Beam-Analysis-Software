#include "Script/TrackAnalysisScript.h"
#include "Algorithm/AnalysisUtils.h"
#include "Detector/DetectorFactory.h"
#include "Event/DataModel.h"
#include "Event/DetectorFrame.h"
#include "Script/Base/RawDataParser.h"
#include "Script/Base/ScriptFactory.h"

#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>

using namespace std;
using AnalysisUtils::FitTrack;
using AnalysisUtils::GetRange;

void TrackAnalysisScript::LoadConfig(const json& config) {
    m_saveValidationData = config.value("saveValidationData", true);
    m_progressInterval = config.value("progressInterval", 10000);
}

void TrackAnalysisScript::Print() const {
    cout << "TrackAnalysisScript Configuration:" << endl;
    cout << "  Save Validation Data: " << (m_saveValidationData ? "Yes" : "No") << endl;
    cout << "  Progress Interval: " << m_progressInterval << endl;
}

bool TrackAnalysisScript::Execute() {
    using Clock = chrono::high_resolution_clock;
    auto t1 = Clock::now();

    cout << "\n========================================" << endl;
    cout << " Track Analysis Script" << endl;
    cout << "========================================" << endl;

    auto parser = GetParser();
    if (!parser) {
        cerr << "Error: Parser not set!" << endl;
        return false;
    }

    auto& factory = DetectorFactory::GetInstance();
    const auto trackers = factory.GetDetectorsByRole(Detector::Role::Tracker);

    // 初始化m_trackerIDs
    m_trackerIDs.clear();
    for (const auto& tracker : trackers) {
        m_trackerIDs.push_back(tracker->GetID());
    }
    sort(m_trackerIDs.begin(), m_trackerIDs.end());
    cout << "[Track Analysis] Initialized with " << m_trackerIDs.size() << " trackers" << endl;

    const Long64_t total = parser->GetTotalEvents();
    std::vector<Event> events;  // Script本地数据
    events.reserve(total);

    // Event Loop
    cout << "\nProcessing events..." << endl;
    for (Long64_t i = 0; i < total; ++i) {
        if ((i % m_progressInterval) == 0 || i == total - 1) {
            cout << "\r[Track Analysis] Processed " << (i + 1) << "/" << total << flush;
        }

        auto rawHits = parser->LoadEvent(i);
        if (rawHits.empty())
            continue;

        Event evt{.eventID = int(i)};
        bool validEvent = true;

        for (auto& det : trackers) {
            auto detEvt = make_shared<DetectorFrame>(*det);
            detEvt->SetRawData(rawHits[det->GetID()]);

            if (!detEvt->Process()) {
                validEvent = false;
                break;
            }
            const int id = det->GetID();
            evt.detectorFramesMap[id] = move(detEvt);
        }

        if (validEvent)
            events.push_back(move(evt));
    }

    cout << "\n[Track Analysis] Valid events: " << events.size() << endl;

    if (events.empty()) {
        cerr << "Error: No valid events!" << endl;
        return false;
    }

    string trackFile = GetOutputDir() + "TrackInfo.root";
    TFile* f = new TFile(trackFile.c_str(), "RECREATE");

    // 直接调用私有方法
    RunTrackerAlign(events, f);

    f->cd();

    // Track Tree
    TTree* tTrack = new TTree("Tracks", "Track info");
    Int_t eventID;
    Track track;
    double t0;
    tTrack->Branch("eventID", &eventID);
    tTrack->Branch("track", &track);
    tTrack->Branch("t0", &t0);

    // Validation Tree
    TTree* tVal = nullptr;
    Int_t detID;
    Double_t resX, resY;
    Double_t hitX, hitY;
    vector<Int_t> clusterIndices;
    vector<StripHit> stripHits;
    vector<Cluster> clusters;

    if (m_saveValidationData) {
        tVal = new TTree("TrackerValidation", "Tracker QA");
        tVal->Branch("eventID", &eventID);
        tVal->Branch("detID", &detID);
        tVal->Branch("resX", &resX);
        tVal->Branch("resY", &resY);
        tVal->Branch("hitX", &hitX);
        tVal->Branch("hitY", &hitY);
        tVal->Branch("clusterIndices", &clusterIndices);
        tVal->Branch("stripHits", &stripHits);
        tVal->Branch("clusters", &clusters);
    }

    int saved = 0;
    int totalEvents = events.size();

    cout << "\nSaving track data..." << endl;
    for (Event& evt : events) {
        vector<double> t0Vec;
        auto [bestTrack, hitIndices, success] = FindBestTrack(evt);

        if (!success)
            continue;

        if (m_saveValidationData) {
            for (const auto& det : trackers) {
                int tid = det->GetID();
                detID = tid;

                int hitIdx = hitIndices[tid];
                const auto& detFrame = evt.detectorFramesMap[tid];

                stripHits = detFrame->StripHits();
                clusters = detFrame->Clusters();
                const LocalHit& localHit = detFrame->LocalHits()[hitIdx];

                for (auto cluster : clusters)
                    t0Vec.push_back(cluster.time);

                auto detector = factory.GetDetector(tid);
                TVector3 predG = detector->CalcHitFromTrack(bestTrack);
                TVector3 predL = detector->GlobalToLocal(predG);

                hitX = localHit.localPos.X();
                hitY = localHit.localPos.Y();
                clusterIndices.clear();
                clusterIndices = localHit.clusterIndices;

                resX = hitX - predL.X();
                resY = hitY - predL.Y();

                tVal->Fill();
            }
        } else {
            for (const auto& det : trackers) {
                const auto& detFrame = evt.detectorFramesMap[det->GetID()];
                for (auto cluster : detFrame->Clusters()) {
                    t0Vec.push_back(cluster.time);
                }
            }
        }

        eventID = evt.eventID;
        track = bestTrack;
        // t0 = accumulate(t0Vec.begin(), t0Vec.end(), 0.0) / t0Vec.size();
        t0 = t0Vec[0];
        tTrack->Fill();
        saved++;
    }

    tTrack->Write();
    if (tVal) tVal->Write();
    f->Close();
    delete f;

    cout << "[Track Analysis] Saved " << saved << " / " << totalEvents << " events" << endl;
    cout << "Output: " << trackFile << endl;

    auto t2 = Clock::now();
    double sec = chrono::duration<double>(t2 - t1).count();

    cout << "\n========================================" << endl;
    cout << "Track Analysis Complete" << endl;
    cout << "Time: " << fixed << setprecision(2) << sec << " seconds" << endl;
    cout << "========================================" << endl;

    return true;
}

// ========== 从TrackAnalysis迁移的方法实现 ==========

map<int, pair<double, double>> TrackAnalysisScript::ComputeTrackError(const vector<Event>& events, TFile* file) {

    auto& factory = DetectorFactory::GetInstance();
    map<int, vector<double>> residX, residY;

    // 计算残差
    for (auto& e : events) {
        vector<GlobalHit> hits;
        for (int tid : m_trackerIDs) {
            auto detector = factory.GetDetector(tid);
            const LocalHit& hit = e.detectorFramesMap.at(tid)->LocalHits().at(0);
            hits.push_back(detector->LocalToGlobal(hit.localPos));
        }

        Track t = FitTrack(hits);

        for (size_t i = 0; i < m_trackerIDs.size(); ++i) {
            int tid = m_trackerIDs[i];
            auto detector = factory.GetDetector(tid);

            GlobalHit pred = detector->CalcHitFromTrack(t);
            GlobalHit meas = hits[i];

            residX[tid].push_back(meas.X() - pred.X());
            residY[tid].push_back(meas.Y() - pred.Y());
        }
    }

    // create directory
    int index = 0;
    std::string dirname = "TrackError";
    while (file->GetDirectory(dirname.c_str()) != nullptr) {
        ++index;
        dirname = "TrackError_" + std::to_string(index);
    }

    file->mkdir(dirname.c_str());
    file->cd(dirname.c_str());

    map<int, pair<double, double>> sigmas;

    for (int tid : m_trackerIDs) {
        auto& vx = residX[tid];
        auto& vy = residY[tid];

        auto [xmin, xmax] = GetRange(vx);
        TH1D* hx = new TH1D(Form("hFinalResX_%d", tid),
                            Form("Tracker %d Final Residual X;#DeltaX [mm];Events", tid),
                            200, xmin, xmax);
        for (double r : vx) hx->Fill(r);
        hx->Fit("gaus", "Q");
        hx->Write();
        sigmas[tid].first = hx->GetFunction("gaus")->GetParameter(2);

        auto [ymin, ymax] = GetRange(vy);
        TH1D* hy = new TH1D(Form("hFinalResY_%d", tid),
                            Form("Tracker %d Final Residual Y;#DeltaY [mm];Events", tid),
                            200, ymin, ymax);
        for (double r : vy) hy->Fill(r);
        hy->Fit("gaus", "Q");
        hy->Write();
        sigmas[tid].second = hy->GetFunction("gaus")->GetParameter(2);
    }

    return sigmas;
}

void TrackAnalysisScript::AlignTrackers(const vector<Event>& events) {
    auto& factory = DetectorFactory::GetInstance();

    const int nParPerDet = 3;  // dx, dy, dRotZ
    const UInt_t nPar = (m_trackerIDs.size() - 1) * nParPerDet;

    auto minim = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    minim->SetTolerance(0.01);
    minim->SetPrintLevel(0);

    // χ² 定义
    auto chi2 = [&](const double* par) {
        for (size_t i = 1; i < m_trackerIDs.size(); ++i) {
            int tid = m_trackerIDs[i];
            double dx = par[(i - 1) * nParPerDet];
            double dy = par[(i - 1) * nParPerDet + 1];
            double dRotZ = par[(i - 1) * nParPerDet + 2];
            auto detector = factory.GetDetector(tid);
            detector->SetAlignment(dx, dy, 0, 0, 0, dRotZ);
        }

        double chi2Sum = 0;
        int nevt = 0;

        for (auto& e : events) {
            vector<GlobalHit> hits;
            for (int tid : m_trackerIDs) {
                auto detector = factory.GetDetector(tid);
                hits.push_back(detector->LocalToGlobal(e.detectorFramesMap.at(tid)->LocalHits().at(0).localPos));
            }

            Track t = FitTrack(hits);
            chi2Sum += t.chi2;
            nevt++;
        }
        return chi2Sum > 0 ? chi2Sum / nevt : 1e9;
    };

    ROOT::Math::Functor f(chi2, nPar);
    minim->SetFunction(f);

    for (UInt_t i = 0; i < nPar; i++)
        minim->SetVariable(i, Form("p%d", i), 0, 0.001);

    minim->Minimize();

    const double* par = minim->X();
    for (size_t i = 1; i < m_trackerIDs.size(); ++i) {
        int tid = m_trackerIDs[i];
        double dx = par[(i - 1) * 3];
        double dy = par[(i - 1) * 3 + 1];
        double dRotZ = par[(i - 1) * 3 + 2];
        auto detector = factory.GetDetector(tid);
        detector->SetAlignment(dx, dy, 0, 0, 0, dRotZ);
    }

    delete minim;
}

pair<double, double> TrackAnalysisScript::ComputePredictionError(int targetDetID) {
    auto& factory = DetectorFactory::GetInstance();

    // Seed tracker IDs
    const int seed1 = m_seedTrackerIDs[0];
    const int seed2 = m_seedTrackerIDs[1];

    // Their intrinsic resolutions (assumed independent)
    const auto [s1x, s1y] = m_sigmaMap[seed1];
    const auto [s2x, s2y] = m_sigmaMap[seed2];

    // Target intrinsic resolution
    const auto [stx, sty] = m_sigmaMap[targetDetID];

    // Z positions
    auto det1 = factory.GetDetector(seed1);
    auto det2 = factory.GetDetector(seed2);
    auto detT = factory.GetDetector(targetDetID);

    const double z1 = det1->GetPos().Z();
    const double z2 = det2->GetPos().Z();
    const double zt = detT->GetPos().Z();

    // Linear interpolation/extrapolation factor α
    const double L = (z2 - z1);
    if (std::abs(L) < 1e-12) {
        std::cerr << "[ComputePredictionError] Seed trackers have identical z!" << std::endl;
        return {0, 0};
    }
    const double alpha = (zt - z1) / L;

    // === Exact track prediction variance (no approximation) ===
    const double varTrackX = (1 - alpha) * (1 - alpha) * s1x * s1x + alpha * alpha * s2x * s2x;
    const double varTrackY = (1 - alpha) * (1 - alpha) * s1y * s1y + alpha * alpha * s2y * s2y;

    // === Combine with target detector intrinsic resolution ===
    const double sigmaPredX = std::sqrt(varTrackX + stx * stx);
    const double sigmaPredY = std::sqrt(varTrackY + sty * sty);

    return {sigmaPredX, sigmaPredY};
}

tuple<Track, map<int, int>, bool> TrackAnalysisScript::FindBestTrack(const Event& event) {
    auto& factory = DetectorFactory::GetInstance();

    // 选择seed tracker, hit数量最小
    vector<pair<int, int>> hitCounts;  // (hitCount, trackerIndex)
    for (size_t i = 0; i < m_trackerIDs.size(); ++i) {
        int tid = m_trackerIDs[i];
        int count = event.detectorFramesMap.at(tid)->LocalHits().size();
        hitCounts.push_back({count, i});
    }

    sort(hitCounts.begin(), hitCounts.end());

    // 选择击中数最少的两个tracker
    int seedTrackerId1 = m_trackerIDs[hitCounts[0].second];
    int seedTrackerId2 = m_trackerIDs[hitCounts[1].second];

    m_seedTrackerIDs = {seedTrackerId1, seedTrackerId2};

    const auto& hits_seed1 = event.detectorFramesMap.at(seedTrackerId1)->LocalHits();
    const auto& hits_seed2 = event.detectorFramesMap.at(seedTrackerId2)->LocalHits();

    auto det1 = factory.GetDetector(seedTrackerId1);
    auto det2 = factory.GetDetector(seedTrackerId2);

    // 贪婪搜索：遍历seed组合
    for (size_t i1 = 0; i1 < hits_seed1.size(); ++i1) {
        for (size_t i2 = 0; i2 < hits_seed2.size(); ++i2) {

            GlobalHit globalHit1 = det1->LocalToGlobal(hits_seed1[i1].localPos);
            GlobalHit globalHit2 = det2->LocalToGlobal(hits_seed2[i2].localPos);

            Track currentTrack = FitTrack({globalHit1, globalHit2});

            // 初始化击中索引数组
            std::map<int, int> hitIndices;
            hitIndices[seedTrackerId1] = i1;
            hitIndices[seedTrackerId2] = i2;

            bool allValid = true;

            for (size_t i = 0; i < m_trackerIDs.size(); ++i) {

                int tid = m_trackerIDs[i];
                if (tid == seedTrackerId1 || tid == seedTrackerId2) continue;

                auto detector = factory.GetDetector(tid);

                // 预测击中位置
                GlobalHit predictedGlobal = detector->CalcHitFromTrack(currentTrack);

                // 计算预测误差
                auto [sigma_pred_X, sigma_pred_Y] = ComputePredictionError(tid);

                // 寻找最近击中
                const auto& hits_i = event.detectorFramesMap.at(tid)->LocalHits();
                double minDist = numeric_limits<double>::infinity();
                int bestIdx = -1;

                for (size_t j = 0; j < hits_i.size(); ++j) {
                    GlobalHit globalHit = detector->LocalToGlobal(hits_i[j].localPos);
                    double resX = globalHit.X() - predictedGlobal.X();
                    double resY = globalHit.Y() - predictedGlobal.Y();

                    double normDist = sqrt((resX / sigma_pred_X) * (resX / sigma_pred_X) + (resY / sigma_pred_Y) * (resY / sigma_pred_Y));

                    if (normDist < minDist) {
                        minDist = normDist;
                        bestIdx = j;
                    }
                }

                if (minDist > 3.0) {
                    allValid = false;
                    break;
                }

                hitIndices[tid] = bestIdx;
            }

            if (allValid) {

                vector<GlobalHit> allGlobalHits;
                for (size_t i = 0; i < m_trackerIDs.size(); ++i) {
                    int tid = m_trackerIDs[i];
                    int hitIdx = hitIndices[tid];
                    LocalHit localHit = event.detectorFramesMap.at(tid)->LocalHits().at(hitIdx);
                    auto detector = factory.GetDetector(tid);
                    GlobalHit globalHit = detector->LocalToGlobal(localHit.localPos);
                    allGlobalHits.push_back(globalHit);
                }

                Track finalTrack = FitTrack(allGlobalHits);

                return {finalTrack, hitIndices, true};
            }
        }
    }

    return {Track{}, std::map<int, int>(), false};
}

void TrackAnalysisScript::RunTrackerAlign(const vector<Event>& events, TFile* file) {
    auto& factory = DetectorFactory::GetInstance();
    const auto& trackers = factory.GetDetectorsByRole(Detector::Role::Tracker);

    cout << "[TrackerAlign] Perform tracker alignment? (y/n): ";

    char choice;
    cin >> choice;
    cin.ignore();

    if (choice != 'y' && choice != 'Y') return;

    // 1) Collect single-hit events
    vector<Event> single;
    single.reserve(events.size());

    for (const auto& evt : events) {
        bool ok = true;
        for (auto& det : trackers) {
            int id = det->GetID();
            for (int type : det->getConfig().readoutPlaneType) {
                if (evt.detectorFramesMap.at(id)->Clusters(type).size() != 1) {
                    ok = false;
                    break;
                }
            }
        }
        if (ok) single.push_back(evt);
    }
    cout << "[TrackerAlign] Found " << single.size() << " single-hit events." << endl;
    if (single.size() < 20) {
        cerr << "[Alignment] WARNING: Too few single-hit events: " << single.size() << endl;
    }

    // 2) Coarse alignment (residual mean shifts)
    int refID = m_trackerIDs[0];

    map<int, vector<double>> resX, resY;

    for (const auto& evt : single) {
        auto refDet = factory.GetDetector(refID);
        auto ref = refDet->LocalToGlobal(evt.detectorFramesMap.at(refID)->LocalHits().at(0).localPos);

        for (size_t i = 1; i < m_trackerIDs.size(); ++i) {
            int tid = m_trackerIDs[i];
            auto detector = factory.GetDetector(tid);
            auto hit = detector->LocalToGlobal(evt.detectorFramesMap.at(tid)->LocalHits().at(0).localPos);
            resX[tid].push_back(hit.X() - ref.X());
            resY[tid].push_back(hit.Y() - ref.Y());
        }
    }

    file->mkdir("Alignment");
    file->cd("Alignment");

    for (size_t i = 1; i < m_trackerIDs.size(); ++i) {
        int tid = m_trackerIDs[i];

        // X
        auto [xmin, xmax] = GetRange(resX[tid]);
        TH1D* hX = new TH1D(Form("hResX_%d", tid),
                            Form("Tracker %d Residual X;#DeltaX [mm];Events", tid),
                            200, xmin, xmax);
        for (double r : resX[tid]) hX->Fill(r);
        hX->Fit("gaus", "Q");
        hX->Write();

        double dx = hX->GetFunction("gaus")->GetParameter(1);

        // Y
        auto [ymin, ymax] = GetRange(resY[tid]);
        TH1D* hY = new TH1D(Form("hResY_%d", tid),
                            Form("Tracker %d Residual Y;#DeltaY [mm];Events", tid),
                            200, ymin, ymax);
        for (double r : resY[tid]) hY->Fill(r);
        hY->Fit("gaus", "Q");
        hY->Write();

        double dy = hY->GetFunction("gaus")->GetParameter(1);

        // apply coarse shifts
        auto detector = factory.GetDetector(tid);
        detector->SetAlignment(dx, dy, 0, 0, 0, 0);
    }

    AlignTrackers(single);

    // 3) first σ measurement
    auto sigmaMap = ComputeTrackError(single, file);

    // 4) filter events by 3σ
    vector<Event> filtered;
    filtered.reserve(single.size());

    for (auto& e : single) {
        bool good = true;

        vector<GlobalHit> hits;
        for (int tid : m_trackerIDs) {
            auto detector = factory.GetDetector(tid);
            hits.push_back(detector->LocalToGlobal(e.detectorFramesMap.at(tid)->LocalHits().at(0).localPos));
        }

        Track t = FitTrack(hits);

        for (size_t i = 0; i < m_trackerIDs.size(); ++i) {
            int tid = m_trackerIDs[i];
            auto detector = factory.GetDetector(tid);
            auto pred = detector->CalcHitFromTrack(t);
            auto meas = hits[i];

            double dx = meas.X() - pred.X();
            double dy = meas.Y() - pred.Y();

            if (fabs(dx) > 3 * sigmaMap[tid].first || fabs(dy) > 3 * sigmaMap[tid].second) {
                good = false;
                break;
            }
        }
        if (good) filtered.push_back(e);
    }

    // 5) second fine alignment
    AlignTrackers(filtered);

    // 6) final σ measurement
    m_sigmaMap = ComputeTrackError(filtered, file);

    cout << "\n[Tracker Alignment] Final detector positions:" << endl;
    for (int tid : m_trackerIDs) {
        auto det = factory.GetDetector(tid);
        TVector3 pos = det->GetPos();
        TVector3 rot = det->GetRot();
        cout << "  " << det->GetName() << " (ID=" << tid << "):" << endl;
        cout << "    \"position\": [" << fixed << setprecision(5)
             << pos.X() << ", " << pos.Y() << ", " << pos.Z() << "]," << endl;
        cout << "    \"rotation\": ["
             << rot.X() << ", " << rot.Y() << ", " << rot.Z() << "]," << endl;
    }

    cout << "\n[Tracker Alignment] done." << endl;
}

REGISTER_SCRIPT("TrackAnalysis", TrackAnalysisScript);
