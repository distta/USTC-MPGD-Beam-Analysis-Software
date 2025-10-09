#include "Pipeline.h"
#include "DataModel.h"

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>

#include "Detector/Cylinder.h"
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

constexpr std::array<int, 14> kBoardToRawIndex = {
    0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6};
// MapBoardChannel: adopt original mapping logic (board pairs -> index -> det/type)
std::tuple<int, int, int> Pipeline::MapBoardChannel(unsigned int boardID, unsigned int channelID) const {

    int rawDataIndex = (boardID < kBoardToRawIndex.size())
                           ? kBoardToRawIndex[boardID]
                           : static_cast<int>(boardID) / 2;

    int type = (rawDataIndex % 2 == 0) ? 0 : 1;
    int detID = (rawDataIndex / 2) + 1;
    int stripID = static_cast<int>(channelID);
    return {detID, stripID, type};
}

bool Pipeline::EventFilter(const std::unordered_map<int, std::vector<RecHit>>& recHits) {

    vector<GlobalHit> globalHits;
    for (const auto& [detID, det] : m_dets) {

        if (det->isDUT()) continue;
        if (recHits.find(detID) == recHits.end())
            return true;

        if (recHits.at(detID).size() != 1)
            return true;

        globalHits.push_back(m_dets[detID]->LocalToGlobal(recHits.at(detID).front()));
    }

    Track t = FitTrack(globalHits);
    if (t.chi2 > 10)
        return true;

    return false;
}

Track Pipeline::FitTrack(const std::vector<GlobalHit>& globalHits) const {
    Track t;
    std::vector<double> zs, xs, ys;
    // collect first hit per detector (if exist). You may change to use centroid or best hit.
    for (const auto& hit : globalHits) {
        zs.push_back(hit.z);
        xs.push_back(hit.x);
        ys.push_back(hit.y);
    }

    const size_t n = zs.size();
    if (n < 2) {
        t.slope_x = t.slope_y = 0.0;
        t.intercept_x = t.intercept_y = 0.0;
        t.chi2 = 1e9;
        return t;
    }

    // Fit x(z)
    double sum_z = 0, sum_z2 = 0, sum_x = 0, sum_zx = 0;
    for (size_t i = 0; i < n; ++i) {
        sum_z += zs[i];
        sum_z2 += zs[i] * zs[i];
        sum_x += xs[i];
        sum_zx += zs[i] * xs[i];
    }
    double denom = n * sum_z2 - sum_z * sum_z;
    if (std::abs(denom) < 1e-12) {
        t.slope_x = 0;
        t.intercept_x = sum_x / n;
    } else {
        t.slope_x = (n * sum_zx - sum_z * sum_x) / denom;
        t.intercept_x = (sum_x - t.slope_x * sum_z) / n;
    }

    // Fit y(z)
    double sum_y = 0, sum_zy = 0;
    for (size_t i = 0; i < n; ++i) {
        sum_y += ys[i];
        sum_zy += zs[i] * ys[i];
    }
    if (std::abs(denom) < 1e-12) {
        t.slope_y = 0;
        t.intercept_y = sum_y / n;
    } else {
        t.slope_y = (n * sum_zy - sum_z * sum_y) / denom;
        t.intercept_y = (sum_y - t.slope_y * sum_z) / n;
    }

    // Compute chi2 as sum of squared residuals
    double chi2 = 0;
    for (size_t i = 0; i < n; ++i) {
        double dx = xs[i] - (t.intercept_x + t.slope_x * zs[i]);
        double dy = ys[i] - (t.intercept_y + t.slope_y * zs[i]);
        chi2 += dx * dx + dy * dy;
    }
    t.chi2 = chi2 / std::max(1.0, static_cast<double>(n - 2));
    return t;
}

double Pipeline::TrackChi2(const double* par) {
    for (size_t i = 0; i < m_trackerIDs.size(); ++i) {
        int detIdx = m_trackerIDs[i];
        int base = i * 3;

        double dx = par[base + 0];
        double dy = par[base + 1];
        double rotZ = par[base + 2];

        m_dets[detIdx]->Alignment(dx, dy, rotZ);
    }

    double chi2 = 0.0;
    int nEvents = 0;

    for (const auto& evt : m_events) {
        std::vector<GlobalHit> points;
        for (int detID : m_trackerIDs) {
            auto it = evt.recLocalHits.find(detID);
            if (it == evt.recLocalHits.end() || it->second.empty()) continue;

            const RecHit& hit = it->second.front();
            points.push_back(m_dets[detID]->LocalToGlobal(hit));
        }
        if (points.size() < 2) continue;

        Track track = FitTrack(points);

        for (auto& point : points) {
            pair<double, double> pos = track.getPositionAtZ(point.z);
            chi2 += std::pow(point.x - pos.first, 2) + std::pow(point.y - pos.second, 2);
        }
        nEvents++;
    }

    return (nEvents > 0) ? chi2 / nEvents : 1e9;
}

double Pipeline::ProfileNLL(const double* par) {
    const size_t nTrackers = m_trackerIDs.size();
    const size_t kDim = 2 * nTrackers;  // 观测维数
    double totalNLL = 0.0;

    /* 惩罚系数：比典型 NLL 大 2~3 个量级即可，但保持可导 */
    const double PENALTY_SIGMA = 1e6;
    const double PENALTY_MISS = 1e5;
    const double EPS = 1e-12;  // 数值正则

    for (const auto& event : m_events) {
        // ---------- 构造权重矩阵 W 和观测向量 obs ----------
        Eigen::VectorXd obs(kDim);
        Eigen::MatrixXd W = Eigen::MatrixXd::Zero(kDim, kDim);
        double logDetSigma = 0.0;  // log|Σ|

        for (size_t i = 0; i < nTrackers; ++i) {
            const double sigmaX = par[2 * i];
            const double sigmaY = par[2 * i + 1];

            /* 平滑惩罚：sigma<=0 */
            if (sigmaX <= 0.0 || sigmaY <= 0.0)
                return PENALTY_SIGMA + PENALTY_SIGMA * (std::abs(sigmaX) + std::abs(sigmaY));

            int detID = m_trackerIDs[i];
            auto it = event.recLocalHits.find(detID);
            if (it == event.recLocalHits.end() || it->second.empty())
                return PENALTY_MISS;

            GlobalHit ahit = m_dets[detID]->LocalToGlobal(it->second.front());
            obs(2 * i) = ahit.x;
            obs(2 * i + 1) = ahit.y;

            const double wx = 1.0 / (sigmaX * sigmaX);
            const double wy = 1.0 / (sigmaY * sigmaY);
            W(2 * i, 2 * i) = wx;
            W(2 * i + 1, 2 * i + 1) = wy;

            logDetSigma += std::log(sigmaX * sigmaX) + std::log(sigmaY * sigmaY);
        }

        // ---------- 设计矩阵 A ----------
        Eigen::MatrixXd A(kDim, 4);
        for (size_t i = 0; i < nTrackers; ++i) {
            double z = m_dets[m_trackerIDs[i]]->GetPosZ();
            A.row(2 * i) << 1, z, 0, 0;
            A.row(2 * i + 1) << 0, 0, 1, z;
        }

        // ---------- 加权最小二乘解析解 ----------
        Eigen::MatrixXd ATW = A.transpose() * W;
        Eigen::MatrixXd ATA = ATW * A;
        Eigen::VectorXd ATy = ATW * obs;

        ATA += Eigen::Matrix4d::Identity() * EPS;  // 正则
        Eigen::VectorXd theta = ATA.ldlt().solve(ATy);

        Eigen::VectorXd res = obs - A * theta;
        double chi2 = res.transpose() * W * res;

        // ---------- 负对数似然 ----------
        totalNLL += 0.5 * (chi2 + logDetSigma + kDim * std::log(2 * M_PI));
    }
    return totalNLL;
}

void Pipeline::Initialize() {

    std::ifstream in(m_configFile);
    if (!in.is_open()) {
        throw std::runtime_error("Pipeline: cannot open config file: " + m_configFile);
    }
    in >> m_config;

    if (!m_config.contains("detectors") || !m_config["detectors"].is_array()) {
        throw std::runtime_error("Pipeline: no detectors defined in config");
    }

    for (const auto& detCfg : m_config["detectors"]) {
        if (!detCfg.contains("id") || !detCfg.contains("name")) {
            throw std::runtime_error("Pipeline: detector missing id or name");
        }
        int detID = detCfg["id"].get<int>();
        std::string name = detCfg["name"].get<std::string>();
        std::string type = detCfg.value("type", "planar");

        if (m_dets.find(detID) != m_dets.end()) {
            throw std::runtime_error("Pipeline: duplicate detector id " + std::to_string(detID));
        }

        // Only 'planar' implemented here (you have Planar class)
        if (type == "planar") {
            // Planar constructor signature: Planar(int id, const std::string& name, const json& cfg)
            m_dets[detID] = std::make_shared<Planar>(detID, name, detCfg);
        } else if (type == "Cylinder") {
            m_dets[detID] = std::make_shared<Cylinder>(detID, name, detCfg);
        } else {
            throw std::runtime_error("Pipeline: unsupported detector type: " + type);
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

void Pipeline::Run(const std::string& rawFile, const std::string& cacheFile, const std::string& outFile) {
    auto t0 = std::chrono::high_resolution_clock::now();

    if (!std::filesystem::exists(cacheFile)) {
        std::cout << "[Pipeline] Cache not found. Running Clustering stage to produce: " << cacheFile << std::endl;
        RunClustering(rawFile, cacheFile);
    } else {
        std::cout << "[Pipeline] Cache file '" << cacheFile << "' already exists. Overwrite? (y/n): ";
        char choice;
        std::cin >> choice;
        if (choice != 'y' && choice != 'Y') {
            std::cout << "[Pipeline] Cache exists: " << cacheFile << " -- skipping clustering stage." << std::endl;
        } else {
            RunClustering(rawFile, cacheFile);
        }
    }

    // Light stage always runs (fast, used for parameter scans)
    std::cout << "[Pipeline] Running light stage (matching, geometry, track fitting) ..." << std::endl;
    RunAnalysis(cacheFile, outFile);

    auto t1 = std::chrono::high_resolution_clock::now();
    double sec = std::chrono::duration<double>(t1 - t0).count();
    std::cout << "[Pipeline] Done. Total time: " << sec << " s\n";
}

void Pipeline::RunClustering(const std::string& rawFile, const std::string& cacheFile) {

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
    std::vector<std::vector<short>>* apv_q = nullptr;
    unsigned int apv_evt = 0;

    // raw data branches
    t_raw->SetBranchAddress("apv_id", &apv_id);
    t_raw->SetBranchAddress("mm_strip", &apv_ch);
    t_raw->SetBranchAddress("apv_q", &apv_q);
    t_raw->SetBranchAddress("apv_evt", &apv_evt);

    TFile* f_out = new TFile(cacheFile.c_str(), "RECREATE");

    UInt_t eventID;
    ULong64_t clusterIndex;
    UInt_t clusterCount;

    // Header tree: eventID, clusterIndex (start), clusterCount
    TTree* headerTree = new TTree("Header", "Index for clusters");
    headerTree->Branch("eventID", &eventID);
    headerTree->Branch("clusterIndex", &clusterIndex);
    headerTree->Branch("clusterCount", &clusterCount);

    Int_t detID;
    RecCluster cluster;

    // Data tree: detID, cluster
    TTree* dataTree = new TTree("Clusters", "clusters");
    dataTree->Branch("detID", &detID);
    dataTree->Branch("cluster", &cluster);

    ULong64_t currentIndex = 0;

    // -------- Event loop --------
    Long64_t nEntries = t_raw->GetEntries();
    dataTree->SetCacheSize(nEntries * m_dets.size() * sizeof(RecCluster) * 2);

    Long64_t goodEvents = 0;
    Long64_t processedEvents = 0;

    for (Long64_t i = 0; i < nEntries; ++i) {
        t_raw->GetEntry(i);
        processedEvents++;
        if (i % 10000 == 0 || i == nEntries - 1) {
            std::cout << "\r[Clustering]: Event " << processedEvents << "/" << nEntries << " | After Filter: " << goodEvents << std::flush;
        }

        // -------- Build RawHits processing --------
        std::unordered_map<int, std::vector<RawData>> rawHits;
        for (size_t j = 0; j < apv_id->size(); ++j) {
            auto [dID, stripID, type] = MapBoardChannel((*apv_id)[j], (*apv_ch)[j]);
            if (m_dets.find(dID) != m_dets.end()) {
                rawHits[dID].push_back(RawData{stripID, type, (*apv_q)[j]});
            }
        }

        //-------- Reconstruction --------
        std::unordered_map<int, std::vector<RecHit>> recHits;
        for (auto& [dID, raws] : rawHits) {
            recHits[dID] = m_dets[dID]->Reconstruction(raws);
        }

        //-------- Event selection --------
        if (EventFilter(recHits)) continue;
        goodEvents++;

        //-------- Clustering output --------
        clusterIndex = currentIndex;
        clusterCount = 0;

        for (const auto& [dID, hits] : recHits) {
            for (const auto& h : hits) {
                for (const auto& c : h.cluster) {
                    detID = dID;
                    cluster = std::move(c);
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

void Pipeline::RunAnalysis(const std::string& cacheFile, const std::string& outFile) {

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
    RecCluster* cluster = nullptr;

    clusterTree->SetBranchAddress("detID", &detID);
    clusterTree->SetBranchAddress("cluster", &cluster);

    int nEvents = headerTree->GetEntries();
    for (int i = 0; i < nEvents; ++i) {
        headerTree->GetEntry(i);

        Event event;
        event.eventID = eventID;

        RecHit aHit;
        int matchIDFlag = -1;
        int detIDFlag = -1;

        for (ULong64_t j = clusterIndex; j < clusterIndex + clusterCount; ++j) {
            clusterTree->GetEntry(j);
            // Process cluster
            if (detID != detIDFlag || cluster->matchID != matchIDFlag) {
                if (!aHit.cluster.empty()) {
                    event.recLocalHits[detIDFlag].push_back(aHit);
                    aHit.cluster.clear();
                }
                detIDFlag = detID;
                matchIDFlag = cluster->matchID;
            }

            aHit.cluster.push_back(*cluster);
        }

        if (!aHit.cluster.empty()) {
            event.recLocalHits[detIDFlag].push_back(aHit);
            aHit.cluster.clear();
        }
        m_events.push_back(std::move(event));
    }

    TFile* f_out = new TFile(outFile.c_str(), "RECREATE");

    RunTrackAnalysis(f_out);

    RunDUTAnalysis(f_out);

    f_out->Close();
}

void Pipeline::RunTrackAnalysis(TFile* outFile) {
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
            det->Alignment(dx, dy, rotZ);

            double xPos = det->GetPosX();
            double yPos = det->GetPosY();
            double zPos = det->GetPosZ();
            double xRot = det->GetRotX();
            double yRot = det->GetRotY();
            double zRot = det->GetRotZ();
            std::cout << "Tracker " << detIdx << " alignment: "
                      << std::fixed << std::setprecision(5)  // 控制浮点数格式
                      << "position: [" << xPos << "," << yPos << "," << zPos << "],"
                      << "rotation: [" << xRot << "," << yRot << "," << zRot << "]"
                      << std::endl;
        }

        std::cout << "[TrackAnalysis] Final chi2: " << minimizer->MinValue() << std::endl;
    }

    // std::cout << "[RunResolutionFit] Starting profile likelihood fit for tracker resolution...\n";

    // UInt_t nPar = m_trackerIDs.size() * 2;  // 每个tracker有 sigma_x, sigma_y

    // auto minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    // minimizer->SetPrintLevel(2);
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

    {
        TTree* residualTree = new TTree("TrackResiduals", "track residuals");
        int eventID;
        int detID;
        double localX, localY;
        double globalX, globalY;
        double resX, resY;
        double trackSlopeX, trackInterceptX;
        double trackSlopeY, trackInterceptY;

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

        for (auto& event : m_events) {
            eventID = event.eventID;
            for (int i = 0; i < m_trackerIDs.size(); i++) {

                int excludedDetID = m_trackerIDs[i];
                detID = excludedDetID;
                vector<GlobalHit> allHits;
                GlobalHit targetHit;

                for (auto& [dtID, hits] : event.recLocalHits) {
                    if (dtID == excludedDetID) {
                        localX = hits[0].cluster[0].pos;
                        localY = hits[0].cluster.size() > 1 ? hits[0].cluster[1].pos : 0.0;
                        targetHit = m_dets[dtID]->LocalToGlobal(hits[0]);
                        continue;
                    }
                    if (m_dets[dtID]->isDUT()) continue;
                    allHits.push_back(m_dets[dtID]->LocalToGlobal(hits[0]));
                }

                if (allHits.size() < 2) continue;

                Track track = FitTrack(allHits);
                event.track = track;

                globalX = targetHit.x;
                globalY = targetHit.y;
                trackSlopeX = track.slope_x;
                trackInterceptX = track.intercept_x;
                trackSlopeY = track.slope_y;
                trackInterceptY = track.intercept_y;

                double predX = trackSlopeX * m_dets[excludedDetID]->GetPosZ() + trackInterceptX;
                double predY = trackSlopeY * m_dets[excludedDetID]->GetPosZ() + trackInterceptY;

                resX = targetHit.x - predX;
                resY = targetHit.y - predY;

                residualTree->Fill();
            }
        }

        outFile->cd();
        residualTree->Write();
    }
}

void Pipeline::RunDUTAnalysis(TFile* outFile) {
    // 创建DUT效率分析树
    TTree* dutTree = new TTree("DUTEfficiency", "DUT efficiency analysis");

    int eventID, dutDetID;
    bool hasTrack, hasHit;
    double trackPredX, trackPredY;
    double hitX, hitY;
    double residualX, residualY;

    dutTree->Branch("eventID", &eventID);
    dutTree->Branch("dutDetID", &dutDetID);
    dutTree->Branch("trackPredX", &trackPredX);
    dutTree->Branch("hitX", &hitX);
    dutTree->Branch("residualX", &residualX);

    // 处理每个事件
    for (auto& event : m_events) {
        eventID = event.eventID;
        Track& track = event.track;

        for (const auto& [detID, det] : m_dets) {
            if (!det->isDUT()) continue;

            dutDetID = detID;
            trackPredX = det->GetLocalHit(track, 0);

            double minDistance = 1e9;

            for (const auto& hit : event.recLocalHits[detID]) {
                double hx = hit.cluster[0].pos;

                double dist = fabs(hx - trackPredX);

                if (dist < minDistance) {
                    minDistance = dist;
                    hitX = hx;
                    residualX = hx - trackPredX;
                }
            }

            dutTree->Fill();
        }
    }

    outFile->cd();
    dutTree->Write();
}