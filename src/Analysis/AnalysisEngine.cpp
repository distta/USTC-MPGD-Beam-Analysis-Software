#include "AnalysisEngine.h"

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TTree.h>

#include "Analysis/TrackAnalysis.h"
#include "DataModel.h"
#include "Detector.h"
#include "Detector/DetectorFactory.h"
#include "Detector/Planar.h"
#include "Event/DetectorFrame.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>

using namespace std;

// ========== 全局辅助函数 ==========

Track FitTrack(const vector<TVector3>& hits) {
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

std::pair<double, double> GetRange(const std::vector<double>& v) {
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

// ========== DUTStatistics 实现 ==========

double DUTStatistics::CalculateMean(const std::vector<double>& values) {
    if (values.empty()) return 0.0;
    
    double sum = 0.0;
    for (double val : values) {
        sum += val;
    }
    return sum / values.size();
}

double DUTStatistics::CalculateRMS(const std::vector<double>& values) {
    if (values.empty()) return 0.0;
    
    double mean = CalculateMean(values);
    double sumSq = 0.0;
    for (double val : values) {
        sumSq += (val - mean) * (val - mean);
    }
    return std::sqrt(sumSq / values.size());
}

std::pair<int, int> DUTStatistics::GetBinIndices(double predX, double predY) const {
    const auto& binning = m_config.binning;
    
    if (predX < binning.predX_min || predX > binning.predX_max || 
        predY < binning.predY_min || predY > binning.predY_max) {
        return {-1, -1};
    }
    
    int binX = static_cast<int>((predX - binning.predX_min) / 
                                 (binning.predX_max - binning.predX_min) * binning.nBinsX);
    int binY = static_cast<int>((predY - binning.predY_min) / 
                                 (binning.predY_max - binning.predY_min) * binning.nBinsY);
    
    // 处理边界情况
    if (binX >= binning.nBinsX) binX = binning.nBinsX - 1;
    if (binY >= binning.nBinsY) binY = binning.nBinsY - 1;
    
    return {binX, binY};
}

void DUTStatistics::AddBinData(int dutID, int binX, int binY, 
                                bool hasValidHit, double resX, double resY) {
    if (binX < 0 || binY < 0) return;
    
    auto& binData = m_binDataMap[dutID][{binX, binY}];
    binData.totalEvents++;
    
    if (hasValidHit) {
        binData.hitEvents++;
        binData.resX_values.push_back(resX);
        binData.resY_values.push_back(resY);
    }
}

// ========== AnalysisEngine辅助函数实现 ==========

// 创建无效Cluster对象
Cluster AnalysisEngine::CreateInvalidCluster(int type) {
    Cluster invalidCluster;
    invalidCluster.type = type;
    invalidCluster.size = DUTAnalysisConfig::kInvalidSize;
    invalidCluster.range = DUTAnalysisConfig::kInvalidSize;
    invalidCluster.charge = DUTAnalysisConfig::kInvalidValue;
    invalidCluster.maxAmp = DUTAnalysisConfig::kInvalidValue;
    invalidCluster.time = DUTAnalysisConfig::kInvalidValue;
    invalidCluster.centroid = DUTAnalysisConfig::kInvalidValue;
    invalidCluster.pos = DUTAnalysisConfig::kInvalidValue;
    invalidCluster.stripHitIndices.clear();
    return invalidCluster;
}

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
    using Clock = std::chrono::high_resolution_clock;
    auto t1 = Clock::now();

    std::cout << "\n========================================\n"
              << " Track Analysis\n"
              << "========================================\n";

    auto& factory = DetectorFactory::GetInstance();
    const auto trackers = factory.GetDetectorsByRole(Detector::Role::Tracker);

    const Long64_t total = m_parser->GetTotalEvents();
    m_events.clear();
    m_events.reserve(total);

    // =======================================
    // Event Loop
    // =======================================
    for (Long64_t i = 0; i < total; ++i) {

        if ((i % 10000) == 0 || i == total - 1) std::cout << "\r[Track Analysis] Processed " << (i + 1) << "/" << total << std::flush;

        auto rawHits = m_parser->LoadEvent(i);
        if (rawHits.empty())
            continue;

        Event evt{.eventID = int(i)};
        bool validEvent = true;

        for (auto& det : trackers) {
            auto detEvt = std::make_shared<DetectorFrame>(*det);

            detEvt->SetRawData(rawHits[det->GetID()]);

            if (!detEvt->Process()) {
                validEvent = false;
                break;
            }
            const int id = det->GetID();
            evt.detectorFramesMap[id] = std::move(detEvt);
        }

        if (validEvent)
            m_events.push_back(std::move(evt));
    }

    std::cout << "\n[Track Analysis] Valid events: " << m_events.size() << "\n";

    if (m_events.empty()) {
        std::cerr << "Error: No valid events!" << std::endl;
        return;
    }

    string trackFile = m_outputDir + "TrackInfo.root";
    TFile* f = new TFile(trackFile.c_str(), "RECREATE");

    TrackAnalysis trackAnalysis(m_outputDir);
    trackAnalysis.RunTrackerAlign(m_events, f);

    f->cd();

    TTree* tTrack = new TTree("Tracks", "Track info");
    Int_t eventID;
    Track track;
    double t0;
    tTrack->Branch("eventID", &eventID);
    tTrack->Branch("track", &track);
    tTrack->Branch("t0", &t0);

    TTree* tVal = new TTree("TrackerValidation", "Tracker QA");
    Int_t detID;
    Double_t resX, resY;
    Double_t hitX, hitY;
    std::vector<Int_t> clusterIndices;
    std::vector<StripHit> stripHits;
    std::vector<Cluster> clusters;
    tVal->Branch("eventID", &eventID);
    tVal->Branch("detID", &detID);
    tVal->Branch("resX", &resX);
    tVal->Branch("resY", &resY);
    tVal->Branch("hitX", &hitX);
    tVal->Branch("hitY", &hitY);
    tVal->Branch("clusterIndices", &clusterIndices);
    tVal->Branch("stripHits", &stripHits);
    tVal->Branch("clusters", &clusters);

    int saved = 0;
    int totalEvents = m_events.size();

    for (Event& evt : m_events) {

        vector<double> t0Vec;
        auto [bestTrack, hitIndices, success] = trackAnalysis.FindBestTrack(evt);

        if (!success)
            continue;

        for (const auto& det : trackers) {
            int tid = det->GetID();
            detID = tid;

            int hitIdx = hitIndices[tid];
            const auto& detFrame = evt.detectorFramesMap[tid];

            // 填充完整数据
            stripHits = detFrame->StripHits();
            clusters = detFrame->Clusters();
            const LocalHit& localHit = detFrame->LocalHits()[hitIdx];

            for(auto cluster:clusters)
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

        eventID = evt.eventID;
        track = bestTrack;
        t0 = std::accumulate(t0Vec.begin(), t0Vec.end(), 0.0) / t0Vec.size();
        tTrack->Fill();
        saved++;
    }

    tTrack->Write();
    tVal->Write();
    f->Close();

    cout << "[Track Analysis] Final Events " << saved << " / " << totalEvents << " events" << endl;

    auto t2 = chrono::high_resolution_clock::now();
    double sec = chrono::duration<double>(t2 - t1).count();

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
    cout << " Loading track info..." << endl;
    cout << " File: " << trackFile << endl;

    TFile* f = TFile::Open(trackFile.c_str(), "READ");
    if (!f || f->IsZombie()) {
        cerr << "Error: Cannot open " << trackFile << endl;
        cerr << "Please run Track Analysis first!" << endl;
        return;
    }

    TTree* trackTree = (TTree*)f->Get("Tracks");
    if (!trackTree) {
        cerr << "Error: No Tracks tree!" << endl;
        f->Close();
        return;
    }

    Int_t eventID;
    Track* track = nullptr;
    trackTree->SetBranchAddress("eventID", &eventID);
    trackTree->SetBranchAddress("track", &track);

    cout << "\nProcessing DUT data..." << endl;
    m_events.clear();

    int processed = 0;
    Long64_t nEntries = trackTree->GetEntries();
    const auto& duts = factory.GetDetectorsByRole(Detector::Role::DUT);

    for (Long64_t i = 0; i < nEntries; ++i) {
        trackTree->GetEntry(i);

        if (processed % 1000 == 0) cout << "\r  Processing: " << processed << "/" << nEntries << flush;
        auto rawHits = m_parser->LoadEvent(eventID);
        if (rawHits.empty()) continue;

        Event evt{.eventID = int(eventID), .track = *track};

        for (auto& det : duts) {
            const int id = det->GetID();
            auto detEvt = std::make_shared<DetectorFrame>(*det);
            detEvt->SetRawData(rawHits[det->GetID()]);

            detEvt->Process();

            evt.detectorFramesMap[id] = std::move(detEvt);
        }

        m_events.push_back(std::move(evt));

        processed++;
    }
    cout << endl;
    cout << "  Processed " << m_events.size() << " DUT events" << endl;

    RunDUTAlign();

    // 初始化配置和统计类
    DUTAnalysisConfig analysisConfig;
    DUTStatistics statistics(analysisConfig);

    // ======== 不再进行全局残差计算，直接在遍历中收集每个bin的残差 ========

    // ======== 第二遍遍历：填充DUTTree并累积统计 ========
    cout << "\nSaving DUT results..." << endl;

    string dutFile = m_outputDir + "DUTInfo.root";
    TFile* fDut = new TFile(dutFile.c_str(), "RECREATE");
    TTree* tDut = new TTree("DUTTree", "DUT data");

    // 定义原有分支
    Int_t dutID;
    Double_t resX, resY, predX, predY;
    Double_t hitX, hitY;
    std::vector<Int_t> clusterIndices;
    std::vector<StripHit> stripHits;
    std::vector<Cluster> clusters;
    
    // 定义新增分支
    Int_t hitFlag;
    Cluster clusterX, clusterY;
    std::vector<StripHit> stripHitsX, stripHitsY;

    // 注册原有分支
    tDut->Branch("eventID", &eventID);
    tDut->Branch("dutID", &dutID);
    tDut->Branch("predX", &predX);
    tDut->Branch("predY", &predY);
    tDut->Branch("resX", &resX);
    tDut->Branch("resY", &resY);
    tDut->Branch("hitX", &hitX);
    tDut->Branch("hitY", &hitY);
    tDut->Branch("clusterIndices", &clusterIndices);
    tDut->Branch("stripHits", &stripHits);
    tDut->Branch("clusters", &clusters);
    
    // 注册新增分支
    tDut->Branch("hitFlag", &hitFlag);
    tDut->Branch("clusterX", &clusterX);
    tDut->Branch("clusterY", &clusterY);
    tDut->Branch("stripHitsX", &stripHitsX);
    tDut->Branch("stripHitsY", &stripHitsY);

    for (auto& det : duts) {
        int id = det->GetID();
        for (auto& evt : m_events) {
            eventID = evt.eventID;
            dutID = id;

            TVector3 predG = det->CalcHitFromTrack(evt.track);
            TVector3 predL = det->GlobalToLocal(predG);
            predX = predL.X();
            predY = predL.Y();

            if (!analysisConfig.IsInValidRange(predX, predY)) continue;

            // 获取bin索引
            auto [binX, binY] = statistics.GetBinIndices(predX, predY);
            
            auto frameIt = evt.detectorFramesMap.find(id);
            if (frameIt != evt.detectorFramesMap.end()) {
                const auto& detFrame = frameIt->second;

                // 填充完整数据（兼容性）
                stripHits = detFrame->StripHits();
                clusters = detFrame->Clusters();

                if (!clusters.empty()) {
                    // 复用CalcuDutResidual查找最佳匹配cluster并计算残差
                    LocalHit localHit = CalcuDutResidual(det, clusters, predL, resX, resY);
                    const auto& clusterIdx = localHit.clusterIndices;
                    
                    // 提取clusterX和clusterY
                    int idxX = clusterIdx[0];
                    int idxY = clusterIdx[1];
                    
                    clusterX = (idxX >= 0) ? clusters[idxX] : CreateInvalidCluster(DUTAnalysisConfig::kTypeX);
                    clusterY = (idxY >= 0) ? clusters[idxY] : CreateInvalidCluster(DUTAnalysisConfig::kTypeY);
                    
                    // 计算hitX和hitY
                    hitX = localHit.localPos.X();
                    hitY = localHit.localPos.Y();
                    
                    // 设置hitFlag
                    bool hasX = (idxX >= 0);
                    bool hasY = (idxY >= 0);
                    if (hasX && hasY) hitFlag = 3;
                    else if (hasX) hitFlag = 1;
                    else if (hasY) hitFlag = 2;
                    else hitFlag = 0;
                    
                    // 提取stripHitsX
                    stripHitsX.clear();
                    if (hasX) {
                        for (int idx : clusterX.stripHitIndices) {
                            if (idx >= 0 && idx < static_cast<int>(stripHits.size())) {
                                stripHitsX.push_back(stripHits[idx]);
                            }
                        }
                    }
                    
                    // 提取stripHitsY
                    stripHitsY.clear();
                    if (hasY) {
                        for (int idx : clusterY.stripHitIndices) {
                            if (idx >= 0 && idx < static_cast<int>(stripHits.size())) {
                                stripHitsY.push_back(stripHits[idx]);
                            }
                        }
                    }
                    
                    // 填充clusterIndices（兼容性）
                    clusterIndices = {idxX, idxY};
                    
                    // 累积分区统计（不进行5sigma筛选，直接收集残差）
                    bool hasValidHit = hasX && hasY;  // 仅判断是否有X和Y cluster
                    statistics.AddBinData(dutID, binX, binY, hasValidHit, resX, resY);
                    
                } else {
                    // 无cluster时设置默认值
                    clusterX = CreateInvalidCluster(DUTAnalysisConfig::kTypeX);
                    clusterY = CreateInvalidCluster(DUTAnalysisConfig::kTypeY);
                    stripHitsX.clear();
                    stripHitsY.clear();
                    hitX = DUTAnalysisConfig::kInvalidValue;
                    hitY = DUTAnalysisConfig::kInvalidValue;
                    resX = resY = 0;
                    clusterIndices = {-1, -1};
                    hitFlag = 0;
                    
                    // 累积总事件数
                    statistics.AddBinData(dutID, binX, binY, false, 0, 0);
                }
            } else {
                // 无数据帧时清空
                stripHits.clear();
                clusters.clear();
                clusterX = CreateInvalidCluster(DUTAnalysisConfig::kTypeX);
                clusterY = CreateInvalidCluster(DUTAnalysisConfig::kTypeY);
                stripHitsX.clear();
                stripHitsY.clear();
                hitX = DUTAnalysisConfig::kInvalidValue;
                hitY = DUTAnalysisConfig::kInvalidValue;
                resX = resY = 0;
                clusterIndices = {-1, -1};
                hitFlag = 0;
            }

            tDut->Fill();
        }
    }

    tDut->Write();
    
    // ======== 创建EfficiencyTree ========
    cout << "\nCalculating efficiency statistics..." << endl;
    TTree* tEff = new TTree("EfficiencyTree", "Efficiency data");
    
    const auto& binning = analysisConfig.binning;
    Int_t eff_dutID, eff_binX, eff_binY;
    Double_t eff_predX_center, eff_predY_center;
    Int_t eff_totalEvents, eff_hitEvents;
    Double_t efficiency, efficiency_error;
    
    tEff->Branch("dutID", &eff_dutID);
    tEff->Branch("binX", &eff_binX);
    tEff->Branch("binY", &eff_binY);
    tEff->Branch("predX_center", &eff_predX_center);
    tEff->Branch("predY_center", &eff_predY_center);
    tEff->Branch("totalEvents", &eff_totalEvents);
    tEff->Branch("hitEvents", &eff_hitEvents);
    tEff->Branch("efficiency", &efficiency);
    tEff->Branch("efficiency_error", &efficiency_error);
    
    TH2D* hEff2D = new TH2D("Efficiency2D", "2D Efficiency;predX (mm);predY (mm);Efficiency", 
                             binning.nBinsX, binning.predX_min, binning.predX_max, 
                             binning.nBinsY, binning.predY_min, binning.predY_max);
    TH2D* hEvents2D = new TH2D("Events2D", "2D Events;predX (mm);predY (mm);Events",
                               binning.nBinsX, binning.predX_min, binning.predX_max,
                               binning.nBinsY, binning.predY_min, binning.predY_max);
    
    const auto& binDataMap = statistics.GetBinDataMap();
    double binWidthX = (binning.predX_max - binning.predX_min) / binning.nBinsX;
    double binWidthY = (binning.predY_max - binning.predY_min) / binning.nBinsY;
    
    struct BinQuality {
        int dutID, binX, binY;
        int events;
        double efficiency, resolution;
        double score;
    };
    std::vector<BinQuality> qualityList;
    
    for (const auto& [dutID_key, binMap] : binDataMap) {
        for (const auto& [binPair, binData] : binMap) {
            eff_dutID = dutID_key;
            eff_binX = binPair.first;
            eff_binY = binPair.second;
            eff_totalEvents = binData.totalEvents;
            eff_predX_center = binning.predX_min + (eff_binX + 0.5) * binWidthX;
            eff_predY_center = binning.predY_min + (eff_binY + 0.5) * binWidthY;
            
            eff_hitEvents = 0;
            double resolutionXY = 0;
            if (!binData.resX_values.empty()) {
                auto [xMin, xMax] = GetRange(binData.resX_values);
                auto [yMin, yMax] = GetRange(binData.resY_values);
                
                TH1D hTempX("hTempX", "", 100, xMin, xMax);
                TH1D hTempY("hTempY", "", 100, yMin, yMax);
                for (double val : binData.resX_values) hTempX.Fill(val);
                for (double val : binData.resY_values) hTempY.Fill(val);
                
                hTempX.Fit("gaus", "Q0");
                hTempY.Fit("gaus", "Q0");
                TF1* fitX = hTempX.GetFunction("gaus");
                TF1* fitY = hTempY.GetFunction("gaus");
                
                double meanX = fitX ? fitX->GetParameter(1) : (xMin + xMax) / 2;
                double sigmaX = fitX ? fitX->GetParameter(2) : (xMax - xMin) / 10;
                double meanY = fitY ? fitY->GetParameter(1) : (yMin + yMax) / 2;
                double sigmaY = fitY ? fitY->GetParameter(2) : (yMax - yMin) / 10;
                
                double cutX = DUTAnalysisConfig::kSigmaFactor * sigmaX;
                double cutY = DUTAnalysisConfig::kSigmaFactor * sigmaY;
                resolutionXY = std::sqrt(sigmaX * sigmaX + sigmaY * sigmaY);
                
                for (size_t i = 0; i < binData.resX_values.size(); ++i) {
                    double resX = binData.resX_values[i], resY = binData.resY_values[i];
                    if (std::abs(resX - meanX) < cutX && std::abs(resY - meanY) < cutY) eff_hitEvents++;
                }
            } else {
                eff_hitEvents = binData.hitEvents;
            }
            
            efficiency = eff_totalEvents > 0 ? static_cast<double>(eff_hitEvents) / eff_totalEvents : 0;
            efficiency_error = eff_totalEvents > 0 ? std::sqrt(efficiency * (1 - efficiency) / eff_totalEvents) : 0;
            
            tEff->Fill();
            hEff2D->SetBinContent(eff_binX + 1, eff_binY + 1, efficiency);
            hEvents2D->SetBinContent(eff_binX + 1, eff_binY + 1, eff_totalEvents);
            
            if (eff_totalEvents > 50 && efficiency > 0.8 && resolutionXY > 0 && resolutionXY < 0.5) {
                double score = eff_totalEvents * efficiency / resolutionXY;
                qualityList.push_back({eff_dutID, eff_binX, eff_binY, eff_totalEvents, efficiency, resolutionXY, score});
            }
        }
    }
    
    std::sort(qualityList.begin(), qualityList.end(), 
              [](const BinQuality& a, const BinQuality& b) { return a.score > b.score; });
    
    int nOutput = std::min(10, (int)qualityList.size());
    cout << "Saving " << nOutput << " best residual histograms..." << endl;
    
    for (int i = 0; i < nOutput; ++i) {
        const auto& bq = qualityList[i];
        auto it1 = binDataMap.find(bq.dutID);
        if (it1 == binDataMap.end()) continue;
        auto it2 = it1->second.find({bq.binX, bq.binY});
        if (it2 == it1->second.end()) continue;
        
        const auto& binData = it2->second;
        auto [xMin, xMax] = GetRange(binData.resX_values);
        auto [yMin, yMax] = GetRange(binData.resY_values);
        
        TH1D hResX(Form("Residual_DUT%d_X_Bin%02d%02d", bq.dutID, bq.binX, bq.binY),
                  Form("DUT%d Bin[%d,%d] ResX (Evt=%d Eff=%.2f);Residual X (mm);Entries", 
                       bq.dutID, bq.binX, bq.binY, bq.events, bq.efficiency),
                  100, xMin, xMax);
        TH1D hResY(Form("Residual_DUT%d_Y_Bin%02d%02d", bq.dutID, bq.binX, bq.binY),
                  Form("DUT%d Bin[%d,%d] ResY (Evt=%d Eff=%.2f);Residual Y (mm);Entries",
                       bq.dutID, bq.binX, bq.binY, bq.events, bq.efficiency),
                  100, yMin, yMax);
        
        for (double val : binData.resX_values) hResX.Fill(val);
        for (double val : binData.resY_values) hResY.Fill(val);
        
        hResX.Fit("gaus", "Q");
        hResY.Fit("gaus", "Q");
        
        hResX.Write();
        hResY.Write();
    }
    
    tEff->Write();
    hEff2D->Write();
    hEvents2D->Write();
    
    cout << "Calculating resolution statistics..." << endl;
    TTree* tRes = new TTree("ResolutionTree", "Resolution data");
    Int_t res_dutID, res_binX, res_binY, res_nHits;
    Double_t res_predX_center, res_predY_center, resX_mean, resX_RMS, resY_mean, resY_RMS;
    
    tRes->Branch("dutID", &res_dutID);
    tRes->Branch("binX", &res_binX);
    tRes->Branch("binY", &res_binY);
    tRes->Branch("predX_center", &res_predX_center);
    tRes->Branch("predY_center", &res_predY_center);
    tRes->Branch("nHits", &res_nHits);
    tRes->Branch("resX_mean", &resX_mean);
    tRes->Branch("resX_RMS", &resX_RMS);
    tRes->Branch("resY_mean", &resY_mean);
    tRes->Branch("resY_RMS", &resY_RMS);
    
    TH2D* hResX2D = new TH2D("ResolutionX2D", "X Resolution 2D;predX (mm);predY (mm);Resolution (mm)",
                             binning.nBinsX, binning.predX_min, binning.predX_max,
                             binning.nBinsY, binning.predY_min, binning.predY_max);
    TH2D* hResY2D = new TH2D("ResolutionY2D", "Y Resolution 2D;predX (mm);predY (mm);Resolution (mm)",
                             binning.nBinsX, binning.predX_min, binning.predX_max,
                             binning.nBinsY, binning.predY_min, binning.predY_max);
    
    for (const auto& [dutID_key, binMap] : binDataMap) {
        for (const auto& [binPair, binData] : binMap) {
            if (binData.resX_values.empty()) continue;
            
            res_dutID = dutID_key;
            res_binX = binPair.first;
            res_binY = binPair.second;
            res_nHits = binData.resX_values.size();
            res_predX_center = binning.predX_min + (res_binX + 0.5) * binWidthX;
            res_predY_center = binning.predY_min + (res_binY + 0.5) * binWidthY;
            
            auto [xMin, xMax] = GetRange(binData.resX_values);
            auto [yMin, yMax] = GetRange(binData.resY_values);
            
            TH1D hTempX("hTempX", "", 100, xMin, xMax);
            TH1D hTempY("hTempY", "", 100, yMin, yMax);
            for (double val : binData.resX_values) hTempX.Fill(val);
            for (double val : binData.resY_values) hTempY.Fill(val);
            
            hTempX.Fit("gaus", "Q0");
            hTempY.Fit("gaus", "Q0");
            TF1* fitX = hTempX.GetFunction("gaus");
            TF1* fitY = hTempY.GetFunction("gaus");
            
            resX_mean = fitX ? fitX->GetParameter(1) : (xMin + xMax) / 2;
            resX_RMS = fitX ? fitX->GetParameter(2) : (xMax - xMin) / 10;
            resY_mean = fitY ? fitY->GetParameter(1) : (yMin + yMax) / 2;
            resY_RMS = fitY ? fitY->GetParameter(2) : (yMax - yMin) / 10;
            
            tRes->Fill();
            hResX2D->SetBinContent(res_binX + 1, res_binY + 1, resX_RMS);
            hResY2D->SetBinContent(res_binX + 1, res_binY + 1, resY_RMS);
        }
    }
    
    tRes->Write();
    hResX2D->Write();
    hResY2D->Write();
    
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

    vector<int> dutIDs = factory.GetDetectorIDsByRole(Detector::Role::DUT);

    if (dutIDs.empty()) {
        cout << "No DUT detectors found, skipping alignment." << endl;
        return;
    }

    cout << "\n[DUT Alignment] Found " << dutIDs.size() << " DUT(s)" << endl;

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

LocalHit AnalysisEngine::CalcuDutResidual(std::shared_ptr<Detector> detector, const std::vector<Cluster>& clusters, const TVector3& predL, double& residualX, double& residualY) {

    double predX = predL.X();
    double predY = predL.Y();

    // 从探测器配置读取参数，消除魔法数字
    const auto& config = detector->getConfig();
    const int typeX = DUTAnalysisConfig::kTypeX;
    const int typeY = DUTAnalysisConfig::kTypeY;
    
    // 从配置读取pitch值
    double pitchX = config.readoutPlanePitch.at(typeX);
    double pitchY = config.readoutPlanePitch.at(typeY);

    // X方向处理：找到最优cluster
    int bestClusterXIndex = -1;
    double minResX = std::numeric_limits<double>::infinity();
    double bestPosX = DUTAnalysisConfig::kInvalidValue;

    for (size_t i = 0; i < clusters.size(); ++i) {
        if (clusters[i].type == typeX) {
            double currentResX = std::abs(clusters[i].pos * pitchX - predX);
            if (currentResX < minResX) {
                minResX = currentResX;
                bestClusterXIndex = static_cast<int>(i);
                bestPosX = clusters[i].pos * pitchX;
            }
        }
    }

    // Y方向处理：找到最优cluster
    int bestClusterYIndex = -1;
    double minResY = std::numeric_limits<double>::infinity();
    double bestPosY = DUTAnalysisConfig::kInvalidValue;

    for (size_t i = 0; i < clusters.size(); ++i) {
        if (clusters[i].type == typeY) {
            double currentResY = std::abs(clusters[i].pos * pitchY - predY);
            if (currentResY < minResY) {
                minResY = currentResY;
                bestClusterYIndex = static_cast<int>(i);
                bestPosY = clusters[i].pos * pitchY;
            }
        }
    }

    // 构建LocalHit
    LocalHit localHit;

    if (bestClusterXIndex != -1) {
        residualX = bestPosX - predX;
    } else {
        bestPosX = DUTAnalysisConfig::kInvalidValue;
        residualX = 0;
    }

    if (bestClusterYIndex != -1) {
        residualY = bestPosY - predY;
    } else {
        bestPosY = DUTAnalysisConfig::kInvalidValue;
        residualY = 0;
    }

    localHit.localPos.SetXYZ(bestPosX, bestPosY, 0);
    localHit.clusterIndices = {bestClusterXIndex, bestClusterYIndex};

    return localHit;
}

double AnalysisEngine::DUTChi2Objective(const double* par, const std::vector<Event>& events, std::shared_ptr<Detector> detector, int detID) {

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

        // 使用detectorFramesMap获取Clusters
        auto frameIt = evt.detectorFramesMap.find(detID);
        if (frameIt == evt.detectorFramesMap.end() || frameIt->second->Clusters().empty()) {
            continue;
        }

        const auto& clusters = frameIt->second->Clusters();
        TVector3 predG = detector->CalcHitFromTrack(evt.track);
        TVector3 predL = detector->GlobalToLocal(predG);
        // 调用新的CalcuDutResidual方法
        LocalHit localHit = CalcuDutResidual(detector, clusters, predL, residualX, residualY);

        double res = residualX * residualX + residualY * residualY;
        if (res > 3) continue;
        chi2 += res;
        nEvents++;
    }

    // 返回平均 χ²
    return (nEvents > 0) ? chi2 / nEvents : 1e9;
}