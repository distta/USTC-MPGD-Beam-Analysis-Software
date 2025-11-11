#include "AnalysisEngine.h"
#include "DataModel.h"

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>

#include "Fit/Fitter.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include <TFile.h>
#include <TMinuit.h>
#include <TTree.h>
#include <TVector3.h>

#include <Eigen/Dense>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <ostream>

using namespace ROOT;
using namespace std;

// -------------------- Assist Function --------------------
constexpr std::array<int, 16> kBoardToRawIndex = {
    0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7};
// MapBoardChannel: adopt original mapping logic (board pairs -> index -> det/type)
std::tuple<int, int, int> AnalysisEngine::MapBoardChannel(unsigned int boardID, unsigned int channelID, unsigned int mm_strip) const {

    int rawDataIndex = (boardID < kBoardToRawIndex.size())
                           ? kBoardToRawIndex[boardID]
                           : static_cast<int>(boardID) / 2;

    int type = (rawDataIndex % 2 == 0) ? 0 : 1;
    int detID = (rawDataIndex / 2) + 1;
    int stripID = static_cast<int>(mm_strip);

    if (boardID == 12)
        stripID = channelID;
    else if (boardID == 13)
        stripID = channelID + 128;
    else if (boardID == 14)
        stripID = 257 - channelID;
    else if (boardID == 15)
        stripID = 129 - channelID;

    return {detID, stripID, type};
}

bool AnalysisEngine::EventFilter() {

    vector<GlobalHit> globalHits;
    for (const auto& [detID, det] : m_dets) {

        if (!det->isTracker()) continue;

        if (det->GetNumOfHits() != 1) return true;

        globalHits.push_back(det->LocalToGlobal(det->GetLocalHits().front()));
    }

    Track t = FitTrack(globalHits);
    if (t.chi2 > 3)
        return true;

    return false;
}

Track AnalysisEngine::FitTrack(const std::vector<GlobalHit>& globalHits) const {
    Track track{};
    const size_t n = globalHits.size();
    if (n < 2) {
        std::cerr << "[FitTrack] Not enough points to fit track." << std::endl;
        return track;
    }

    // --- Step 1: compute mean ---
    double sumX = 0, sumY = 0, sumZ = 0;
    for (const auto& p : globalHits) {
        sumX += p.X();
        sumY += p.Y();
        sumZ += p.Z();
    }
    const double meanX = sumX / n;
    const double meanY = sumY / n;
    const double meanZ = sumZ / n;

    // --- Step 2: linear regression x(z), y(z) ---
    double Szz = 0, Szx = 0, Szy = 0;
    for (const auto& p : globalHits) {
        const double dz = p.Z() - meanZ;
        Szz += dz * dz;
        Szx += dz * (p.X() - meanX);
        Szy += dz * (p.Y() - meanY);
    }

    track.kx = Szx / Szz;
    track.ky = Szy / Szz;
    track.bx = meanX - track.kx * meanZ;
    track.by = meanY - track.ky * meanZ;

    // --- Step 3: compute chi2 ---
    double chi2 = 0;
    for (const auto& p : globalHits) {
        const double x_fit = track.kx * p.Z() + track.bx;
        const double y_fit = track.ky * p.Z() + track.by;
        chi2 += (p.X() - x_fit) * (p.X() - x_fit) + (p.Y() - y_fit) * (p.Y() - y_fit);
    }
    track.chi2 = chi2 / (2 * n - 4);  // reduced χ²

    return track;
}

double AnalysisEngine::TrackChi2(const double* par) {
    // --- 1. 应用对齐参数到每个 tracker ---
    for (size_t i = 0; i < m_trackerIDs.size(); ++i) {
        int detID = m_trackerIDs[i];
        if (m_dets.find(detID) == m_dets.end()) continue;

        const double dx = par[i * 3 + 0];
        const double dy = par[i * 3 + 1];
        const double rotZ = par[i * 3 + 2];

        m_dets[detID]->SetAlignment(dx, dy, 0.0, 0.0, 0.0, rotZ);
    }

    // --- 2. 对所有事件求平均 χ² ---
    double chi2 = 0.0;
    int nEvents = 0;

    for (const auto& evt : m_events) {
        std::vector<GlobalHit> globalHits;

        for (int detID : m_trackerIDs) {
            auto it = evt.recLocalHits.find(detID);
            if (it == evt.recLocalHits.end() || it->second.empty()) continue;

            const LocalHit& lhit = it->second.front();
            GlobalHit ghit = m_dets.at(detID)->LocalToGlobal(lhit);
            globalHits.emplace_back(ghit);
        }

        if (globalHits.size() < 2) continue;  // 不足两点无法拟合

        Track track = FitTrack(globalHits);
        chi2 += track.chi2;
        nEvents++;
    }

    // --- 3. 返回平均 χ² ---
    return (nEvents > 0) ? chi2 / nEvents : 1e9;
}

double AnalysisEngine::DUTChi2(const double* par, int detID) {

    const double dx = par[0];
    const double dy = par[1];
    const double rotZ = par[2];

    auto det = m_dets.at(detID);
    det->SetAlignment(dx, dy, 0.0, 0.0, 0.0, rotZ);

    // --- 2. 对所有事件求平均 χ² ---
    double chi2 = 0.0;
    int nEvents = 0;

    for (const auto& evt : m_events) {
        GlobalHit predHit = det->GetHitFromTrack(evt.track);
        // 转换到局部坐标系
        LocalHit localPredHit = det->GlobalToLocal(predHit);
        double trackPredX = localPredHit.X();
        double trackPredY = localPredHit.Y();

        // 检查DUT是否有击中
        auto hitIt = evt.recLocalHits.find(detID);
        bool hasHit = (hitIt != evt.recLocalHits.end() && !hitIt->second.empty());

        if (hasHit) {
            // 找到最近的击中点
            double residualX = 0.0;
            double residualY = 0.0;

            double minResidual = std::numeric_limits<double>::infinity();

            for (const auto& localHit : hitIt->second) {
                double currentResidualX = localHit.X() - trackPredX;
                double currentResidualY = localHit.Y() - trackPredY;
                double currentResidual = std::sqrt(currentResidualX * currentResidualX + currentResidualY * currentResidualY);

                if (currentResidual < minResidual) {
                    minResidual = currentResidual;
                    residualX = currentResidualX;
                    residualY = currentResidualY;
                }
            }

            chi2 += (residualX * residualX + residualY * residualY);
            nEvents++;
        }
    }

    // --- 3. 返回平均 χ² ---
    return (nEvents > 0) ? chi2 / nEvents : 1e9;
}

double AnalysisEngine::ProfileNLL(const double* par) {
    const size_t nTrackers = m_trackerIDs.size();
    const size_t kDim = 2 * nTrackers;
    double totalNLL = 0.0;

    const double PENALTY_MISS = 1e4;
    const double EPS = 1e-9;

    // --- σ 正则化参数 ---
    const double sigmaPriorMean = 0.05;   // 物理预期分辨率 (mm)
    const double sigmaPriorWidth = 0.02;  // 允许波动范围 (mm)
    const double priorWeight = 1.0;       // 正则强度系数，可调 0.5~5.0

    for (const auto& event : m_events) {
        Eigen::VectorXd obs(kDim);
        Eigen::VectorXd w(kDim);
        bool validEvent = true;

        for (size_t i = 0; i < nTrackers; ++i) {
            int detID = m_trackerIDs[i];
            auto it = event.recLocalHits.find(detID);
            if (it == event.recLocalHits.end() || it->second.empty()) {
                totalNLL += PENALTY_MISS;
                validEvent = false;
                break;
            }

            double sigmaX = par[2 * i];
            double sigmaY = par[2 * i + 1];

            // --- 下界约束防止σ→0 ---
            const double MIN_SIGMA = 1e-3;
            if (sigmaX < MIN_SIGMA) sigmaX = MIN_SIGMA;
            if (sigmaY < MIN_SIGMA) sigmaY = MIN_SIGMA;

            GlobalHit ahit = m_dets.at(detID)->LocalToGlobal(it->second.front());
            obs(2 * i) = ahit.x();
            obs(2 * i + 1) = ahit.y();
            w(2 * i) = 1.0 / (sigmaX * sigmaX);
            w(2 * i + 1) = 1.0 / (sigmaY * sigmaY);
        }
        if (!validEvent) continue;

        // --- 设计矩阵 A ---
        Eigen::MatrixXd A(kDim, 4);
        for (size_t i = 0; i < nTrackers; ++i) {
            double z = m_dets.at(m_trackerIDs[i])->GetPos().z();
            A.row(2 * i) << 1, z, 0, 0;
            A.row(2 * i + 1) << 0, 0, 1, z;
        }

        // --- 加权最小二乘 ---
        Eigen::MatrixXd W = w.asDiagonal();
        Eigen::MatrixXd ATA = A.transpose() * W * A + Eigen::Matrix4d::Identity() * EPS;
        Eigen::VectorXd ATy = A.transpose() * W * obs;
        Eigen::VectorXd theta = ATA.ldlt().solve(ATy);

        // --- 残差 ---
        Eigen::VectorXd res = obs - A * theta;
        double chi2 = 0.0;
        for (size_t i = 0; i < kDim; ++i)
            chi2 += res(i) * res(i) * w(i);

        double logDetSigma = 0.0;
        for (size_t i = 0; i < nTrackers; ++i) {
            double sigmaX = par[2 * i];
            double sigmaY = par[2 * i + 1];
            logDetSigma += std::log(2 * M_PI * sigmaX * sigmaX) + std::log(2 * M_PI * sigmaY * sigmaY);
        }

        totalNLL += 0.5 * (chi2 + logDetSigma);
    }

    // ---------- σ 先验正则项 ----------
    for (size_t i = 0; i < nTrackers; ++i) {
        double sigmaX = par[2 * i];
        double sigmaY = par[2 * i + 1];
        totalNLL += 0.5 * priorWeight * (std::pow((sigmaX - sigmaPriorMean) / sigmaPriorWidth, 2) + std::pow((sigmaY - sigmaPriorMean) / sigmaPriorWidth, 2));
    }

    return totalNLL;
}

// -------------------- Main Analysis --------------------
void AnalysisEngine::Initialize() {

    std::ifstream in(m_configFile);
    if (!in.is_open()) {
        throw std::runtime_error("AnalysisEngine: cannot open config file: " + m_configFile);
    }
    in >> m_config;

    if (!m_config.contains("detectors") || !m_config["detectors"].is_array()) {
        throw std::runtime_error("AnalysisEngine: no detectors defined in config");
    }

    for (const auto& detCfg : m_config["detectors"]) {
        if (!detCfg.contains("id") || !detCfg.contains("name")) {
            throw std::runtime_error("AnalysisEngine: detector missing id or name");
        }
        int detID = detCfg["id"].get<int>();
        std::string name = detCfg["name"].get<std::string>();
        std::string type = detCfg.value("type", "planar");

        if (m_dets.find(detID) != m_dets.end()) {
            throw std::runtime_error("AnalysisEngine: duplicate detector id " + std::to_string(detID));
        }

        // Only 'planar' implemented here (you have Planar class)
        if (type == "planar") {
            // Planar constructor signature: Planar(int id, const std::string& name, const json& cfg)
            m_dets[detID] = std::make_shared<Planar>(detID, name, detCfg);
        } else {
            throw std::runtime_error("AnalysisEngine: unsupported detector type: " + type);
        }
    }

    for (auto& [detID, det] : m_dets) {
        if (det->isTracker()) {
            m_trackerIDs.push_back(detID);
        }
    }

    if (m_trackerIDs.size() < 2) {
        std::cerr << "[RunTrackAnalysis] Need at least 2 trackers!" << std::endl;
        return;
    }
}

void AnalysisEngine::Run(const std::string& rawFile, const std::string& cacheFile, const std::string& outFile) {
    auto t0 = std::chrono::high_resolution_clock::now();

    if (!std::filesystem::exists(cacheFile)) {
        std::cout << "[AnalysisEngine] Cluster file not found. Running Clustering stage to produce: " << cacheFile << std::endl;
        RunClustering(rawFile, cacheFile);
    } else {
        std::cout << "[AnalysisEngine] Cluster file '" << cacheFile << "' already exists. Overwrite? (y/n): ";
        char choice;
        std::cin >> choice;
        if (choice != 'y' && choice != 'Y') {
            std::cout << "[AnalysisEngine] Cluster file exists: " << cacheFile << " -- skipping clustering stage." << std::endl;
        } else {
            RunClustering(rawFile, cacheFile);
        }
    }

    // Light stage always runs (fast, used for parameter scans)
    std::cout << "[AnalysisEngine] Running Analysis stage" << std::endl;
    RunAnalysis(cacheFile, outFile);

    auto t1 = std::chrono::high_resolution_clock::now();
    double sec = std::chrono::duration<double>(t1 - t0).count();
    std::cout << "[AnalysisEngine] Done. Total time: " << sec << " s\n";
}

void AnalysisEngine::RunClustering(const std::string& rawFile, const std::string& cacheFile) {

    TFile* f_in = TFile::Open(rawFile.c_str(), "READ");
    if (!f_in || f_in->IsZombie()) {
        std::cerr << "Failed to open raw file: " << rawFile << std::endl;
        return;
    }
    TTree* t_raw = (TTree*)f_in->Get("raw");
    if (!t_raw) {
        std::cerr << "raw TTree not found!" << std::endl;
        return;
    }

    std::vector<unsigned int>* apv_id = nullptr;
    std::vector<unsigned int>* apv_ch = nullptr;
    std::vector<unsigned int>* mm_strip = nullptr;
    std::vector<std::vector<short>>* apv_q = nullptr;
    unsigned int apv_evt = 0;

    t_raw->SetBranchAddress("apv_id", &apv_id);
    t_raw->SetBranchAddress("apv_ch", &apv_ch);
    t_raw->SetBranchAddress("mm_strip", &mm_strip);
    t_raw->SetBranchAddress("apv_q", &apv_q);
    t_raw->SetBranchAddress("apv_evt", &apv_evt);

    TFile* f_out = new TFile(cacheFile.c_str(), "RECREATE");

    UInt_t eventID;
    ULong64_t clusterIndex;
    UInt_t clusterCount;

    TTree* headerTree = new TTree("Header", "Index for clusters");
    headerTree->Branch("eventID", &eventID);
    headerTree->Branch("clusterIndex", &clusterIndex);
    headerTree->Branch("clusterCount", &clusterCount);

    Int_t detID;
    Cluster cluster;
    TTree* dataTree = new TTree("Clusters", "clusters");
    dataTree->Branch("detID", &detID);
    dataTree->Branch("cluster", &cluster);

    ULong64_t currentIndex = 0;

    // ---------- Event Loop ----------
    Long64_t nEntries = t_raw->GetEntries();
    Long64_t goodEvents = 0;
    Long64_t processedEvents = 0;

    for (Long64_t i = 0; i < nEntries; ++i) {
        t_raw->GetEntry(i);
        processedEvents++;

        if (i % 10000 == 0 || i == nEntries - 1) {
            std::cout << "\r[Clustering]: Event " << processedEvents
                      << "/" << nEntries
                      << " | After Filter: " << goodEvents << std::flush;
        }

        // -------- 清空 Detector 内数据 --------
        for (auto& [id, det] : m_dets) det->ClearData();

        // -------- Step 1: 将 RawData 按 detector 分类 --------
        for (size_t j = 0; j < apv_id->size(); ++j) {
            auto [dID, stripID, type] = MapBoardChannel((*apv_id)[j], (*apv_ch)[j], (*mm_strip)[j]);
            if (m_dets.find(dID) == m_dets.end()) continue;
            m_dets[dID]->AddRawData(RawData{stripID, type, (*apv_q)[j]});
        }

        // -------- Step 2: Tracker Reconstruction --------
        for (auto& [dID, det] : m_dets) {
            if (det->isTracker()) det->Reconstruct();
        }

        // -------- Step 3: Tracker-based Event Filter --------
        if (EventFilter()) continue;
        goodEvents++;

        // -------- Step 4: DUT Reconstruction --------
        for (auto& [dID, det] : m_dets) {
            if (det->isDUT()) det->Reconstruct();
        }

        // -------- Step 5: Clustering Output --------
        clusterIndex = currentIndex;
        clusterCount = 0;

        for (auto& [dID, det] : m_dets) {
            for (const auto& recCluster : det->GetRecClusters()) {
                for (const auto& c : recCluster) {
                    detID = dID;
                    cluster = c;
                    dataTree->Fill();
                    clusterCount++;
                    currentIndex++;
                }
            }
        }

        if (clusterCount > 0) {
            eventID = apv_evt;
            headerTree->Fill();
        }
    }

    headerTree->Write();
    dataTree->Write();
    f_out->Close();
    f_in->Close();

    std::cout << std::endl;
    std::cout << "[Clustering]: Clustering successfully. Total events after clustering: " << goodEvents << std::endl;
}

void AnalysisEngine::RunAnalysis(const std::string& cacheFile, const std::string& outFile) {

    TFile* inFile = TFile::Open(cacheFile.c_str(), "READ");
    TTree* headerTree = (TTree*)inFile->Get("Header");
    TTree* clusterTree = (TTree*)inFile->Get("Clusters");
    if (!headerTree || !clusterTree) {
        std::cerr << "Missing required trees in cache file" << std::endl;
        return;
    }

    ULong64_t clusterIndex;
    UInt_t clusterCount;
    UInt_t eventID;

    headerTree->SetBranchAddress("eventID", &eventID);
    headerTree->SetBranchAddress("clusterIndex", &clusterIndex);
    headerTree->SetBranchAddress("clusterCount", &clusterCount);

    Int_t detID;
    Cluster* cluster = nullptr;

    clusterTree->SetBranchAddress("detID", &detID);
    clusterTree->SetBranchAddress("cluster", &cluster);

    int nEvents = headerTree->GetEntries();
    for (int i = 0; i < nEvents; ++i) {
        headerTree->GetEntry(i);

        Event event;
        event.eventID = eventID;

        RecCluster aRecCluster;
        int matchIDFlag = -1;
        int detIDFlag = -1;

        for (ULong64_t j = clusterIndex; j < clusterIndex + clusterCount; ++j) {
            clusterTree->GetEntry(j);
            // Process cluster
            if (detID != detIDFlag || cluster->matchID != matchIDFlag) {
                if (!aRecCluster.empty()) {
                    event.recLocalHits[detIDFlag].push_back(m_dets[detID]->GetLocalHitFromCluster(aRecCluster));
                    aRecCluster.clear();
                }
                detIDFlag = detID;
                matchIDFlag = cluster->matchID;
            }

            aRecCluster.push_back(*cluster);
        }

        if (!aRecCluster.empty()) {
            event.recLocalHits[detIDFlag].push_back(m_dets[detID]->GetLocalHitFromCluster(aRecCluster));
            aRecCluster.clear();
        }
        m_events.push_back(std::move(event));
    }

    TFile* f_out = new TFile(outFile.c_str(), "RECREATE");

    RunTrackAnalysis(f_out);

    RunDUTAnalysis(f_out);

    f_out->Close();
}

void AnalysisEngine::RunTrackAnalysis(TFile* outFile) {
    if (m_events.empty()) {
        std::cerr << "[RunTrackAnalysis] No events to analyze!" << std::endl;
        return;
    }

    std::cout << "[RunTrackAnalysis] Alignment? (y/n): ";
    char choice;
    std::cin >> choice;
    if (choice == 'y' || choice == 'Y') {

        auto minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
        minimizer->SetTolerance(0.01);
        minimizer->SetPrintLevel(0);

        UInt_t nPar = m_trackerIDs.size() * 3;

        auto chi2Func = [this](const double* par) -> double {
            return this->TrackChi2(const_cast<double*>(par));
        };

        ROOT::Math::Functor f(chi2Func, nPar);
        minimizer->SetFunction(f);

        for (UInt_t i = 0; i < nPar; ++i) {
            minimizer->SetVariable(i, Form("p%d", i), 0.0, 0.001);
        }

        minimizer->Minimize();

        // 应用结果
        for (size_t i = 0; i < m_trackerIDs.size(); ++i) {
            int detIdx = m_trackerIDs[i];
            const double* result = minimizer->X();

            double dx = result[i * 3];
            double dy = result[i * 3 + 1];
            double rotZ = result[i * 3 + 2];

            auto& det = m_dets[detIdx];
            det->SetAlignment(dx, dy, 0, 0, 0, rotZ);

            TVector3 pos = det->GetPos();
            TVector3 rot = det->GetRot();

            std::cout << "Tracker " << detIdx << " alignment: "
                      << std::fixed << std::setprecision(5)  // 控制浮点数格式
                      << "\"position\": [" << pos.X() << "," << pos.Y() << "," << pos.Z() << "],"
                      << "\"rotation\": [" << rot.X() << "," << rot.Y() << "," << rot.Z() << "]"
                      << std::endl;
        }

        std::cout << "[TrackAnalysis] Final chi2: " << minimizer->MinValue() << std::endl;
    }

    // std::cout << "[RunResolutionFit] Starting profile likelihood fit for tracker resolution...\n";

    // UInt_t nPar = m_trackerIDs.size() * 2;  // 每个tracker有 sigma_x, sigma_y

    // auto minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    // minimizer->SetPrintLevel(0);
    // minimizer->SetTolerance(1e-3);

    // auto nllFunc = [this](const double* par) -> double {
    //     return this->ProfileNLL(par);
    // };
    // ROOT::Math::Functor f(nllFunc, nPar);
    // minimizer->SetFunction(f);

    // // 初值设置
    // for (UInt_t i = 0; i < nPar / 2; ++i) {
    //     minimizer->SetLimitedVariable(2 * i, Form("sigmaX_det%d", i), 0.1, 0.01, 1e-4, 10.0);
    //     minimizer->SetLimitedVariable(2 * i + 1, Form("sigmaY_det%d", i), 0.1, 0.01, 1e-4, 10.0);
    // }

    // minimizer->Minimize();

    // const double* bestFit = minimizer->X();
    // std::cout << "\n[Result] Tracker resolutions:\n";
    // for (size_t i = 0; i < m_trackerIDs.size(); ++i) {
    //     std::cout << "Tracker " << m_trackerIDs[i]
    //               << " sigmaX=" << bestFit[2 * i]
    //               << " sigmaY=" << bestFit[2 * i + 1] << std::endl;
    // }

    TTree* residualTree = new TTree("TrackResiduals", "track residuals");
    int eventID;
    int detID;
    double localX, localY;
    double globalX, globalY;
    double resX, resY;
    double trackSlopeX, trackInterceptX;
    double trackSlopeY, trackInterceptY;
    double chi2;

    residualTree->Branch("eventID", &eventID);
    residualTree->Branch("detID", &detID);
    residualTree->Branch("localX", &localX);
    residualTree->Branch("localY", &localY);
    residualTree->Branch("globalX", &globalX);
    residualTree->Branch("globalY", &globalY);
    residualTree->Branch("resX", &resX);
    residualTree->Branch("resY", &resY);
    residualTree->Branch("trackSlopeX", &trackSlopeX);
    residualTree->Branch("trackInterceptX", &trackInterceptX);
    residualTree->Branch("trackSlopeY", &trackSlopeY);
    residualTree->Branch("trackInterceptY", &trackInterceptY);
    residualTree->Branch("chi2", &chi2);

    for (auto& event : m_events) {
        eventID = event.eventID;
        vector<GlobalHit> allHits;
        GlobalHit targetHit;

        for (size_t i = 0; i < m_trackerIDs.size(); i++) {

            allHits.clear();
            int excludedDetID = m_trackerIDs[i];
            detID = excludedDetID;

            for (const auto& [dtID, hits] : event.recLocalHits) {
                if (dtID == excludedDetID) {
                    if (!hits.empty()) {
                        localX = hits[0].X();
                        localY = hits[0].Y();
                        targetHit = m_dets[dtID]->LocalToGlobal(hits[0]);
                    }
                    continue;
                }
                if (m_dets[dtID]->isDUT()) continue;
                if (!hits.empty()) {
                    allHits.push_back(m_dets[dtID]->LocalToGlobal(hits[0]));
                }
            }

            if (allHits.size() < 2) continue;

            Track track = FitTrack(allHits);

            globalX = targetHit.X();
            globalY = targetHit.Y();
            trackSlopeX = track.kx;
            trackInterceptX = track.bx;
            trackSlopeY = track.ky;
            trackInterceptY = track.by;
            chi2 = event.track.chi2;

            double predX = trackSlopeX * m_dets[excludedDetID]->GetPos().Z() + trackInterceptX;
            double predY = trackSlopeY * m_dets[excludedDetID]->GetPos().Z() + trackInterceptY;

            resX = targetHit.X() - predX;
            resY = targetHit.Y() - predY;

            residualTree->Fill();
        }
        allHits.push_back(targetHit);
        event.track = FitTrack(allHits);
    }

    outFile->cd();
    residualTree->Write();
}

void AnalysisEngine::RunDUTAnalysis(TFile* outFile) {
    if (m_events.empty()) {
        std::cerr << "[RunDUTAnalysis] No events to analyze!" << std::endl;
        return;
    }

    // First, ask if alignment is needed
    std::cout << "[RunDUTAnalysis] DUT alignment? (y/n): ";
    char choice;
    std::cin >> choice;

    if (choice == 'y' || choice == 'Y') {
        // Loop over each DUT to align them individually
        for (auto& [detID, det] : m_dets) {
            if (!det->isDUT()) continue;  // Skip non-DUT detectors

            std::cout << "[DUT " << detID << "] Aligning..." << std::endl;
            auto minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
            minimizer->SetTolerance(0.01);
            minimizer->SetPrintLevel(0);

            // Number of parameters: dx, dy, rotZ for each DUT
            UInt_t nPar = 3;

            // Define the chi2 function for alignment minimization
            auto chi2Func = [this, detID](const double* par) -> double {
                return this->DUTChi2(const_cast<double*>(par), detID);  // Pass detID to chi2 function
            };

            ROOT::Math::Functor f(chi2Func, nPar);
            minimizer->SetFunction(f);

            // Set initial guesses for alignment (dx, dy, rotZ)
            for (UInt_t i = 0; i < nPar; ++i) {
                minimizer->SetVariable(i, Form("p%d", i), 0.0, 0.001);  // Initial values for dx, dy, rotZ
            }

            minimizer->Minimize();

            // Apply the alignment results
            const double* result = minimizer->X();
            double dx = result[0];
            double dy = result[1];
            double rotZ = result[2];

            det->SetAlignment(dx, dy, 0, 0, 0, rotZ);

            TVector3 pos = det->GetPos();
            TVector3 rot = det->GetRot();

            std::cout << "DUT " << detID << " alignment: "
                      << std::fixed << std::setprecision(5)
                      << "\"position\": [" << pos.X() << "," << pos.Y() << "," << pos.Z() << "],"
                      << "\"rotation\": [" << rot.X() << "," << rot.Y() << "," << rot.Z() << "]"
                      << std::endl;

            std::cout << "[DUT " << detID << "] Final chi2: " << minimizer->MinValue() << std::endl;
        }
    }

    // Create DUT efficiency analysis tree
    TTree* dutTree = new TTree("DUTEfficiency", "DUT efficiency analysis");

    int eventID, dutDetID;
    bool hasTrack, hasHit;
    double trackPredX, trackPredY;
    double hitX, hitY;
    double residualX, residualY;

    dutTree->Branch("eventID", &eventID);
    dutTree->Branch("dutDetID", &dutDetID);
    dutTree->Branch("hasTrack", &hasTrack);
    dutTree->Branch("hasHit", &hasHit);
    dutTree->Branch("trackPredX", &trackPredX);
    dutTree->Branch("trackPredY", &trackPredY);
    dutTree->Branch("hitX", &hitX);
    dutTree->Branch("hitY", &hitY);
    dutTree->Branch("resX", &residualX);
    dutTree->Branch("resY", &residualY);

    // Now loop through each DUT and analyze its efficiency with respect to the events
    for (auto& [detID, det] : m_dets) {
        if (!det->isDUT()) continue;  // Only process DUT detectors

        dutDetID = detID;

        // Loop through events for each DUT
        for (const auto& event : m_events) {
            eventID = event.eventID;
            const Track& track = event.track;

            hasTrack = (track.chi2 >= 0 && track.chi2 < 1e6);  // Valid track check

            // Get predicted hit on DUT plane from the track
            GlobalHit predHit = det->GetHitFromTrack(track);
            LocalHit localPredHit = det->GlobalToLocal(predHit);
            trackPredX = localPredHit.X();
            trackPredY = localPredHit.Y();

            // Check if the DUT has a hit in the current event
            auto hitIt = event.recLocalHits.find(detID);
            hasHit = (hitIt != event.recLocalHits.end() && !hitIt->second.empty());

            if (hasHit) {
                // Find the closest hit (minimize residual)
                double minResidual = std::numeric_limits<double>::infinity();
                const LocalHit* closestHit = nullptr;

                for (const auto& localHit : hitIt->second) {
                    double currentResidualX = localHit.X() - trackPredX;
                    double currentResidualY = localHit.Y() - trackPredY;
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
                }
            } else {
                hitX = -999.0;
                hitY = -999.0;
                residualX = -999.0;
                residualY = -999.0;
            }

            // Fill tree with data for this event and DUT
            dutTree->Fill();
        }
    }

    outFile->cd();
    dutTree->Write();
}
