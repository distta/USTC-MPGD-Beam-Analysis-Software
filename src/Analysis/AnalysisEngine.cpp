#include "AnalysisEngine.h"

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>

#include "Analysis/TrackAnalysis.h"
#include "DataModel.h"
#include "Detector.h"
#include "Detector/DetectorFactory.h"
#include "Detector/Planar.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>

using namespace std;

// ========== 全局辅助函数 ==========

Track FitTrack(const vector<GlobalHit>& hits) {
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

auto GetRange = [](const std::vector<double>& v) {
    if (v.size() < 3) return std::make_pair(0.0, 1.0);
    const double k = 5;

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
};

// ========== AnalysisEngine实现 ==========

AnalysisEngine::AnalysisEngine(const string& cfg, const string& rawDir,
                               const string& resultDir, const string& runID)
    : m_configFile(cfg), m_rawDir(rawDir), m_resultDir(resultDir), m_runID(runID) {

    // 构造输出目录
    m_outputDir = m_resultDir + "/" + m_runID + "/";
    filesystem::create_directories(m_outputDir);
}

void AnalysisEngine::Initialize() {
    cout << "\n========================================" << endl;
    cout << "Initializing BeamAnalysis" << endl;
    cout << "========================================" << endl;

    // 读取配置
    ifstream in(m_configFile);
    if (!in.is_open()) throw runtime_error("Cannot open config: " + m_configFile);
    in >> m_config;

    if (!m_config.contains("detectors"))
        throw runtime_error("No detectors in config");

    // 使用DetectorFactory创建探测器
    auto& factory = DetectorFactory::GetInstance();
    if (!factory.Initialize(m_config)) {
        throw runtime_error("Failed to initialize DetectorFactory");
    }

    // 构造原始数据文件路径
    string rawFile = m_rawDir + "/run" + m_runID + ".root";

    cout << "Config file: " << m_configFile << endl;
    cout << "Raw file   : " << rawFile << endl;
    cout << "Output dir : " << m_outputDir << endl;
    cout << "Run ID     : " << m_runID << endl;

    // 初始化parser
    m_parser = make_shared<RawDataParser>(rawFile);
    if (!m_parser->Initialize())
        throw runtime_error("Failed to initialize parser");

    cout << "Total events: " << m_parser->GetTotalEvents() << endl;
    cout << "Initialization complete" << endl;
}

void AnalysisEngine::RunTrackAnalysis() {
    auto t0 = chrono::high_resolution_clock::now();
    auto& factory = DetectorFactory::GetInstance();

    cout << "\n========================================" << endl;
    cout << "Track Analysis" << endl;
    cout << "========================================" << endl;

    Long64_t total = m_parser->GetTotalEvents();
    m_events.clear();

    // 获取所有Tracker
    auto trackers = factory.GetDetectorsByRole(Detector::Role::Tracker);

    for (Long64_t i = 0; i < total; ++i) {
        if (i % 10000 == 0 || i == total - 1) {
            cout << "\r[Track Analysis] Processed Event: " << i + 1 << "/" << total << flush;
        }

        for (const auto& det : trackers) det->ClearData();

        if (!m_parser->LoadEvent(i)) continue;

        for (const auto& det : trackers) {
            det->Reconstruct();
        }

        Event evt;
        evt.eventID = i;
        bool valid = true;

        for (const auto& det : trackers) {
            auto hits = det->GetLocalHits();
            if (hits.empty()) {
                valid = false;
                break;
            }
            int tid = det->GetID();
            evt.recLocalHits[tid] = hits;
            evt.recClusters[tid] = det->GetRecClusters();
        }

        if (valid) m_events.push_back(move(evt));
    }
    cout << endl;
    cout << "[Track Analysis] " << m_events.size() << " valid events (all tracker has at least one Hit)" << endl;

    if (m_events.empty()) {
        cerr << "Error: No valid events!" << endl;
        return;
    }

    string trackFile = m_outputDir + "TrackInfo.root";
    TFile* f = new TFile(trackFile.c_str(), "RECREATE");

    // 创建TrackAnalysis实例并执行对齐
    TrackAnalysis trackAnalysis(m_outputDir, m_runID);
    trackAnalysis.RunTrackerAlign(m_events, f);

    f->cd();

    TTree* tTrack = new TTree("Tracks", "Track info");
    Int_t eventID;
    Track track;
    tTrack->Branch("eventID", &eventID);
    tTrack->Branch("track", &track);

    TTree* tVal = new TTree("TrackerValidation", "Tracker QA");
    Int_t detID;
    Double_t hitX, hitY, resX, resY;
    Cluster cluster0, cluster1;
    tVal->Branch("eventID", &eventID);
    tVal->Branch("detID", &detID);
    tVal->Branch("hitX", &hitX);
    tVal->Branch("hitY", &hitY);
    tVal->Branch("resX", &resX);
    tVal->Branch("resY", &resY);
    tVal->Branch("cluster0", &cluster0);
    tVal->Branch("cluster1", &cluster1);

    int saved = 0;
    int totalEvents = m_events.size();

    for (auto& evt : m_events) {

        auto [bestTrack, hitIndices, success] = trackAnalysis.FindBestTrack(evt);

        if (!success) {
            continue;
        }

        eventID = evt.eventID;
        track = bestTrack;
        tTrack->Fill();
        saved++;

        for (const auto& det : trackers) {
            int tid = det->GetID();
            detID = tid;

            int hitIdx = hitIndices[tid];
            LocalHit hit = evt.recLocalHits[tid][hitIdx];

            auto detector = factory.GetDetector(tid);
            GlobalHit pred = detector->GetHitFromTrack(bestTrack);
            LocalHit predL = detector->GlobalToLocal(pred);

            hitX = hit.X();
            hitY = hit.Y();
            resX = hit.X() - predL.X();
            resY = hit.Y() - predL.Y();

            cluster0 = Cluster();
            cluster1 = Cluster();
            if (!evt.recClusters[tid].empty()) {
                auto& rc = evt.recClusters[tid][hitIdx];
                if (rc.size() > 0) cluster0 = rc[0];
                if (rc.size() > 1) cluster1 = rc[1];
            }
            tVal->Fill();
        }
    }

    tTrack->Write();
    tVal->Write();
    f->Close();
    delete f;

    cout << "[Track Analysis] Final Events " << saved << " / " << totalEvents << " events" << endl;

    auto t1 = chrono::high_resolution_clock::now();
    double sec = chrono::duration<double>(t1 - t0).count();

    cout << "\n========================================" << endl;
    cout << "Track Analysis Complete" << endl;
    cout << "Time: " << fixed << setprecision(2) << sec << " seconds" << endl;
    cout << "========================================" << endl;
}

void AnalysisEngine::RunDUTAnalysis() {
    auto t0 = chrono::high_resolution_clock::now();
    auto& factory = DetectorFactory::GetInstance();

    cout << "\n========================================" << endl;
    cout << "DUT Analysis" << endl;
    cout << "========================================" << endl;

    // 加载track信息
    string trackFile = m_outputDir + "TrackInfo.root";
    cout << "\n[Step 1] Loading track info..." << endl;
    cout << "  File: " << trackFile << endl;

    TFile* f = TFile::Open(trackFile.c_str(), "READ");
    if (!f || f->IsZombie()) {
        cerr << "Error: Cannot open " << trackFile << endl;
        cerr << "Please run Track Analysis first!" << endl;
        return;
    }

    TTree* t = (TTree*)f->Get("Tracks");
    if (!t) {
        cerr << "Error: No Tracks tree!" << endl;
        f->Close();
        return;
    }

    Int_t eventID;
    Track* track;
    t->SetBranchAddress("eventID", &eventID);
    t->SetBranchAddress("track", &track);

    map<int, Track> trackMap;
    for (Long64_t i = 0; i < t->GetEntries(); ++i) {
        t->GetEntry(i);
        trackMap[eventID] = *track;
    }
    f->Close();

    cout << "  Loaded " << trackMap.size() << " tracks" << endl;

    // 处理DUT
    cout << "\n[Step 2] Processing DUT data..." << endl;
    m_events.clear();

    int processed = 0;
    for (auto& [evtID, track] : trackMap) {
        if (processed % 10000 == 0)
            cout << "\r  Processing: " << processed << "/" << trackMap.size() << flush;

        auto& allDets = factory.GetAllDetectors();
        for (auto& [id, det] : allDets) det->ClearData();
        if (!m_parser->LoadEvent(evtID)) continue;

        for (auto& [id, det] : allDets) {
            if (det->isDUT()) det->Reconstruct();
        }

        Event evt;
        evt.eventID = evtID;
        evt.track = track;

        for (auto& [id, det] : allDets) {
            if (det->isDUT()) {
                evt.recLocalHits[id] = det->GetLocalHits();
                evt.recClusters[id] = det->GetRecClusters();
            }
        }
        m_events.push_back(move(evt));
        processed++;
    }
    cout << endl;
    cout << "  Processed " << m_events.size() << " DUT events" << endl;

    // Step 2.5: DUT对齐（新增）
    RunDUTAlign();

    // 保存DUT结果
    cout << "\n[Step 3] Saving DUT results..." << endl;

    string dutFile = m_outputDir + "DUTInfo.root";
    TFile* fDut = new TFile(dutFile.c_str(), "RECREATE");
    TTree* tDut = new TTree("DUTTree", "DUT data");

    Int_t dutID;
    Bool_t hasHit;
    Double_t hitX, hitY, resX, resY;
    Cluster cluster0, cluster1;

    tDut->Branch("eventID", &eventID);
    tDut->Branch("dutID", &dutID);
    tDut->Branch("hasHit", &hasHit);
    tDut->Branch("hitX", &hitX);
    tDut->Branch("hitY", &hitY);
    tDut->Branch("resX", &resX);
    tDut->Branch("resY", &resY);
    tDut->Branch("cluster0", &cluster0);
    tDut->Branch("cluster1", &cluster1);

    auto& allDets = factory.GetAllDetectors();
    for (auto& [id, det] : allDets) {
        if (!det->isDUT()) continue;

        dutID = id;

        for (auto& evt : m_events) {
            eventID = evt.eventID;
            hasHit = false;
            hitX = hitY = resX = resY = 0;
            cluster0 = Cluster();
            cluster1 = Cluster();

            GlobalHit predHit = det->GetHitFromTrack(evt.track);
            LocalHit localPredHit = det->GlobalToLocal(predHit);

            // Check if the DUT has a hit in the current event
            auto hitIt = evt.recLocalHits.find(id);
            hasHit = (hitIt != evt.recLocalHits.end() && !hitIt->second.empty());
            int hitId = -1;
            int count = -1;

            if (hasHit) {
                // Find the closest hit (minimize residual)
                double minResidual = std::numeric_limits<double>::infinity();
                const LocalHit* closestHit = nullptr;

                for (const auto& localHit : hitIt->second) {
                    double currentResidualX = localHit.X() == -999 ? 0 : localHit.X() - localPredHit.X();
                    double currentResidualY = localHit.Y() == -999 ? 0 : localHit.Y() - localPredHit.Y();
                    double currentResidual = std::sqrt(currentResidualX * currentResidualX + currentResidualY * currentResidualY);
                    count++;

                    if (currentResidual < minResidual) {
                        minResidual = currentResidual;
                        closestHit = &localHit;
                        resX = currentResidualX;
                        resY = currentResidualY;
                        hitId = count;
                    }
                }

                if (minResidual > 3) continue;

                if (closestHit) {
                    hitX = closestHit->X();
                    hitY = closestHit->Y();
                }
            } else {
                hitX = -999.0;
                hitY = -999.0;
                resX = 0;
                resY = 0;
            }

            if (hitId == -1) continue;
            auto itC = evt.recClusters.find(id);
            if (itC != evt.recClusters.end() && !itC->second.empty()) {
                auto& rc = itC->second[hitId];
                if (rc.size() > 0) cluster0 = rc[0];
                if (rc.size() > 1) cluster1 = rc[1];
            }

            tDut->Fill();
        }
    }

    tDut->Write();
    fDut->Close();
    delete fDut;

    cout << "  Output: " << dutFile << endl;

    auto t1 = chrono::high_resolution_clock::now();
    double sec = chrono::duration<double>(t1 - t0).count();

    cout << "\n========================================" << endl;
    cout << "DUT Analysis Complete" << endl;
    cout << "Time: " << fixed << setprecision(2) << sec << " seconds" << endl;
    cout << "========================================" << endl;
}

void AnalysisEngine::RunDUTAlign() {
    auto& factory = DetectorFactory::GetInstance();

    cout << "\nPerform DUT alignment? (y/n): ";
    char choice;
    cin >> choice;
    cin.ignore();

    if (choice != 'y' && choice != 'Y') {
        cout << "Skipping DUT alignment." << endl;
        return;
    }

    // 识别所有DUT
    vector<int> dutIDs;
    auto& allDets = factory.GetAllDetectors();
    for (auto& [id, det] : allDets) {
        if (det->isDUT()) {
            dutIDs.push_back(id);
        }
    }

    if (dutIDs.empty()) {
        cout << "No DUT detectors found, skipping alignment." << endl;
        return;
    }

    cout << "\n[DUT Alignment] Found " << dutIDs.size() << " DUT(s)" << endl;

    // 对每个DUT执行6参数对齐
    for (int dutID : dutIDs) {
        auto detector = factory.GetDetector(dutID);

        if (m_events.empty()) {
            cerr << "[DUT Alignment] No events to analyze for DUT " << dutID << "!" << endl;
            continue;
        }

        cout << "[DUT " << dutID << "] Aligning (6-parameter)..." << endl;
        auto minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
        minimizer->SetTolerance(0.005);
        minimizer->SetPrintLevel(0);

        UInt_t nPar = 6;

        // 使用lambda捕获this和参数
        auto chi2Func = [this, &detector, dutID](const double* par) -> double {
            return this->DUTChi2Objective(par, this->m_events, detector, dutID);
        };

        ROOT::Math::Functor f(chi2Func, nPar);
        minimizer->SetFunction(f);

        for (UInt_t i = 0; i < nPar; ++i) {
            minimizer->SetVariable(i, Form("p%d", i), 0.0, 0.001);
        }

        minimizer->Minimize();

        // 应用结果
        const double* result = minimizer->X();
        double dx = result[0];
        double dy = result[1];
        double dz = result[2];
        double rotX = result[3];
        double rotY = result[4];
        double rotZ = result[5];

        detector->SetAlignment(dx, dy, dz, rotX, rotY, rotZ);

        TVector3 pos = detector->GetPos();
        TVector3 rot = detector->GetRot();

        cout << "DUT " << dutID << " alignment: "
             << fixed << setprecision(5)
             << "\"position\": [" << pos.X() << "," << pos.Y() << "," << pos.Z() << "],"
             << "\"rotation\": [" << rot.X() << "," << rot.Y() << "," << rot.Z() << "]"
             << endl;
        cout << "DUT " << dutID << " alignment corrections: "
             << "dx=" << dx << ", dy=" << dy << ", dz=" << dz << ", "
             << "rotX=" << rotX << ", rotY=" << rotY << ", rotZ=" << rotZ
             << endl;

        cout << "[DUT " << dutID << "] Final chi2: " << minimizer->MinValue() << endl;

        delete minimizer;
    }

    cout << "[DUT Alignment] Complete" << endl;
}

// ========== DUT对齐私有方法 ==========

bool AnalysisEngine::CalcuDutResidual(
    std::shared_ptr<Detector> detector,
    const std::vector<LocalHit>& hits,
    const Track& track,
    double& hitX, double& hitY,
    double& residualX, double& residualY) {

    if (hits.empty()) {
        return false;
    }

    // 使用track预测DUT位置
    GlobalHit predHit = detector->GetHitFromTrack(track);
    LocalHit localPredHit = detector->GlobalToLocal(predHit);

    // 找到最近的hit（最小化残差）
    double minResidual = std::numeric_limits<double>::infinity();
    const LocalHit* closestHit = nullptr;

    for (const auto& localHit : hits) {
        double currentResidualX = localHit.X() == -999 ? 0 : localHit.X() - localPredHit.X();
        double currentResidualY = localHit.Y() == -999 ? 0 : localHit.Y() - localPredHit.Y();
        double currentResidual = std::sqrt(currentResidualX * currentResidualX + currentResidualY * currentResidualY);

        if (currentResidual < minResidual) {
            minResidual = currentResidual;
            closestHit = &localHit;
            residualX = currentResidualX;
            residualY = currentResidualY;
        }
    }

    if (closestHit) {
        hitX = closestHit->X();
        hitY = closestHit->Y();
        return true;
    }

    return false;
}

double AnalysisEngine::DUTChi2Objective(
    const double* par,
    const std::vector<Event>& events,
    std::shared_ptr<Detector> detector,
    int detID) {

    const double dx = par[0];
    const double dy = par[1];
    const double dz = par[2];
    const double rotX = par[3];
    const double rotY = par[4];
    const double rotZ = par[5];

    detector->SetAlignment(dx, dy, dz, rotX, rotY, rotZ);

    // 对所有事件求平均 χ²
    double chi2 = 0.0;
    int nEvents = 0;

    for (const auto& evt : events) {
        double residualX = 0.0;
        double residualY = 0.0;
        double hitX, hitY;

        auto hitIt = evt.recLocalHits.find(detID);
        bool hasHit = (hitIt == evt.recLocalHits.end()) ? false : CalcuDutResidual(detector, hitIt->second, evt.track, hitX, hitY, residualX, residualY);

        double res = std::sqrt(residualX * residualX + residualY * residualY);
        if (res > 2) continue;

        if (hasHit) {
            chi2 += residualX * residualX + residualY * residualY;
            nEvents++;
        }
    }

    // 返回平均 χ²
    return (nEvents > 0) ? chi2 / nEvents : 1e9;
}
