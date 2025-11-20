#include "DataModel.h"
#include "EventDisplayManager.h"

#include "TMultiGraph.h"
#include <TCanvas.h>
#include <TFile.h>
#include <TGaxis.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TMarker.h>
#include <TPad.h>
#include <TTree.h>

#include <algorithm>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>

EventDisplayManager::EventDisplayManager(const std::string& rawDir,
                                         const std::string& resultDir,
                                         const std::string& runID)
    : m_rawDir(rawDir), m_resultDir(resultDir), m_runID(runID) {
    m_rawFilePath = m_rawDir + "/run" + m_runID + ".root";
    m_trackFilePath = m_resultDir + "/" + m_runID + "/TrackInfo.root";
    m_outBaseDir = m_resultDir + "/" + m_runID + "/EventDisplay/";
}

EventDisplayManager::~EventDisplayManager() {
    for (auto* c : m_canvases)
        if (c) delete c;
    m_canvases.clear();
}

// 初始化 RawDataParser 并加载 track entries
bool EventDisplayManager::Initialize() {
    // 初始化 parser
    m_parser = std::make_unique<RawDataParser>(m_rawFilePath);
    if (!m_parser->Initialize()) {
        std::cerr << "[EventDisplay] RawDataParser 初始化失败: " << m_rawFilePath << std::endl;
        return false;
    }

    // Load track entries from TrackInfo.root
    if (!LoadTrackEntries()) {
        std::cerr << "[EventDisplay] LoadTrackEntries 失败: " << m_trackFilePath << std::endl;
        return false;
    }

    // 创建输出目录（存在 AnalysisEngine 时可复用其 output dir）
    std::filesystem::create_directories(m_outBaseDir);
    std::cout << "[EventDisplay] ready. track entries: " << m_trackEntries.size() << std::endl;
    return true;
}

bool EventDisplayManager::LoadTrackEntries() {
    m_trackEntries.clear();
    TFile* f = TFile::Open(m_trackFilePath.c_str(), "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "[EventDisplay] 无法打开 TrackInfo: " << m_trackFilePath << std::endl;
        return false;
    }

    TTree* t = (TTree*)f->Get("Tracks");
    if (!t) {
        std::cerr << "[EventDisplay] Tracks tree 不存在 in " << m_trackFilePath << std::endl;
        f->Close();
        return false;
    }

    Int_t eventID;
    Track* trackPtr = nullptr;
    t->SetBranchAddress("eventID", &eventID);
    t->SetBranchAddress("track", &trackPtr);

    for (Long64_t i = 0; i < t->GetEntries(); ++i) {
        t->GetEntry(i);
        if (!trackPtr) continue;
        TrackEntry te;
        te.eventID = eventID;
        te.track = *trackPtr;
        m_trackEntries.push_back(te);
    }
    f->Close();
    return true;
}

void EventDisplayManager::RunInteractive() {
    if (m_trackEntries.empty()) {
        std::cerr << "[EventDisplay] 没有 track 条目，先运行 RunTrackAnalysis / RunDUTAnalysis。" << std::endl;
        return;
    }

    while (true) {
        std::cout << "\nAvailable track entries: " << m_trackEntries.size() << std::endl;
        std::cout << "请输入要显示的 event id, 或 -1 退出: ";
        int idx;
        std::cin >> idx;
        if (idx == -1) break;

        if (idx < 0 || idx >= (int)m_trackEntries.back().eventID) {
            std::cerr << "选择超出范围" << std::endl;
            continue;
        }

        TrackEntry te;
        for (const auto evt : m_trackEntries) {
            if (evt.eventID == idx) {
                te = evt;
                break;
            }
        }

        if (!ProcessEntry(te)) {
            std::cerr << "[EventDisplay] ProcessEntry 失败 for event " << te.eventID << std::endl;
            continue;
        }
        std::cout << "[EventDisplay] 已生成 DUT overview for event " << te.eventID << std::endl;

        // 进入波形查询循环
        auto& factory = DetectorFactory::GetInstance();
        auto allDets = factory.GetAllDetectors();

        while (true) {
            std::cout << "\n查询波形: 输入 DUT ID, type, target stripID (空格分隔), 或 -1 退出: ";
            int dutID, type, targetStrip;
            std::cin >> dutID;
            if (dutID == -1) break;
            std::cin >> type >> targetStrip;

            auto it = allDets.find(dutID);
            if (it == allDets.end() || !it->second->isDUT()) {
                std::cerr << "无效 DUT ID" << std::endl;
                break;
            }

            auto& rawData = it->second->GetRawData();

            // 找到目标 strip 的相邻 strip（前后各一个，如果存在）
            std::vector<int> nearbyStrips;
            int target = targetStrip;

            if (FindRawForStrip(rawData, targetStrip, type) != nullptr)
                nearbyStrips.push_back(targetStrip);

            // 向左寻找
            int dis = 0;
            int strip = targetStrip - 1;
            while (dis < 2 && strip >= 0) {
                if (FindRawForStrip(rawData, strip, type) != nullptr) {
                    nearbyStrips.push_back(strip);
                    dis = 0;  // reset consecutive missing count
                } else {
                    dis++;
                }
                strip--;
            }

            // 向右寻找
            dis = 0;
            strip = targetStrip + 1;
            while (dis < 2 && strip < 512) {  // 假设最大 strip ID=511
                if (FindRawForStrip(rawData, strip, type) != nullptr) {
                    nearbyStrips.push_back(strip);
                    dis = 0;
                } else {
                    dis++;
                }
                strip++;
            }

            std::sort(nearbyStrips.begin(), nearbyStrips.end());

            if (nearbyStrips.empty()) {
                std::cerr << "未找到目标 strip 及相邻 strip" << std::endl;
                continue;
            }

            // 绘制多条波形到同一 canvas
            TCanvas* c = new TCanvas(Form("wave_evt%d_dut%d_type%d", te.eventID, dutID, type),
                                     Form("Event %d DUT %d type %d", te.eventID, dutID, type),
                                     1000, 600);
            m_canvases.push_back(c);
            gPad->SetGrid(1, 1);

            TMultiGraph* mg = new TMultiGraph();
            std::vector<int> colors = {kBlue, kGreen + 2, kOrange, kMagenta, kCyan, kGray + 1};

            TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
            legend->SetBorderSize(0);
            legend->SetFillStyle(0);

            int colorIdx = 0;
            for (int sid : nearbyStrips) {
                const RawData* r = FindRawForStrip(rawData, sid, type);
                if (!r) continue;

                TGraph* g = new TGraph();
                for (size_t i = 0; i < r->adc.size(); ++i)
                    g->SetPoint(i, (double)i, (double)r->adc[i]);

                g->SetLineWidth(2);
                g->SetMarkerSize(2);
                g->SetLineColor(colors[colorIdx % colors.size()]);
                mg->Add(g, "APL");

                legend->AddEntry(g, Form("strip %d", sid), "l");

                colorIdx++;
            }

            mg->Draw("A L");
            mg->SetTitle(Form("Event %d DUT %d Type %d;Sample;ADC", te.eventID, dutID, type));
            legend->Draw();

            // 保存输出
            std::string outDir = m_outBaseDir + "Event_" + std::to_string(te.eventID) + "/";
            std::filesystem::create_directories(outDir);
            c->SaveAs((outDir + Form("wave_nearby_dut%d_type%d_strip%d.png", dutID, type, targetStrip)).c_str());

            std::cout << "[EventDisplay] 波形已保存: target strip " << targetStrip << " 与相邻 strip\n";
        }
    }
}

// 处理单个 entry：LoadEvent -> Reconstruct -> 对每个 DUT 绘图并写入 ROOT
bool EventDisplayManager::ProcessEntry(const TrackEntry& te) {
    auto& factory = DetectorFactory::GetInstance();
    auto allDets = factory.GetAllDetectors();

    // 1) clear detectors, load raw event (m_parser::LoadEvent 会把 RawData 分发到 Detectors via DetectorFactory)
    for (auto& [id, det] : allDets) det->ClearData();
    if (!m_parser->LoadEvent(te.eventID)) {
        std::cerr << "[EventDisplay] RawDataParser::LoadEvent failed for event " << te.eventID << std::endl;
        return false;
    }

    // 2) reconstruct DUTs
    for (auto& [id, det] : allDets) {
        if (det->isDUT()) det->Reconstruct();
    }

    // 3) prepare output ROOT file for this event

    std::string outDir = m_outBaseDir + "Event_" + std::to_string(te.eventID) + "/";
    std::filesystem::create_directories(outDir);

    std::string outFile = outDir + "Event_" + std::to_string(te.eventID) + ".root";
    TFile* fout = TFile::Open(outFile.c_str(), "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "[EventDisplay] 无法创建输出文件: " << outFile << std::endl;
        return false;
    }

    // 4) for each DUT: draw overview and save into its directory
    for (auto& [id, det] : allDets) {
        if (!det->isDUT()) continue;

        DrawDUTOverview(te.eventID, det, te.track);
    }

    fout->Write();
    fout->Close();
    delete fout;
    return true;
}

void EventDisplayManager::DrawDUTOverview(int eventID, std::shared_ptr<Detector> det, const Track& track) {

    const int DUT_ID = det->GetID();
    const std::string detName = det->GetName();

    auto stripHitsMap = det->GetStripHits();  // map<int, vector<StripHit>>

    if (stripHitsMap.empty()) {
        TLatex note;
        note.SetTextSize(0.04);
        note.DrawLatexNDC(0.15, 0.5, Form("DUT %s (ID=%d) : no strip hits for event %d", detName.c_str(), DUT_ID, eventID));
        TH1F placeholder(Form("placeholder_evt%d_dut%d", eventID, DUT_ID), "no data", 1, 0, 1);
        placeholder.Write();
        return;
    }

    std::vector<int> types;
    for (const auto& kv : stripHitsMap) types.push_back(kv.first);

    int nTypes = static_cast<int>(types.size());
    TCanvas* c = new TCanvas(Form("DUT_%d_evt%d_overview", DUT_ID, eventID),
                             Form("DUT %d - event %d overview", DUT_ID, eventID),
                             1200, std::max(500, nTypes * 300));
    m_canvases.push_back(c);
    c->Divide(1, nTypes, 0.005, 0.005);

    const int STRIP_XMIN = 0;
    const int STRIP_XMAX = 512;
    const double TIME_YMIN = 0;
    const double TIME_YMAX = 5;

    GlobalHit globalHit = det->GetHitFromTrack(track);
    LocalHit localHit = det->GlobalToLocal(globalHit);
    for (auto hit : det->GetLocalHits()) {
        std::cout << "DUT " << DUT_ID << " Local Hit: (" << hit.X() << ", " << hit.Y() << ")" << std::endl;
    }
    std::cout << localHit.X() << "," << localHit.Y() << std::endl;

    for (int iType = 0; iType < nTypes; ++iType) {
        int type = types[iType];
        const auto& hits = stripHitsMap[type];

        c->cd(iType + 1);
        gPad->SetGrid(1, 1);

        // Collect strip amplitude and time
        std::map<int, double> stripAmp, stripTime;
        for (const auto& sh : hits) {
            stripAmp[sh.stripID] = sh.amp;
            stripTime[sh.stripID] = (sh.time - 80) * 0.02;
        }

        // --- Amplitude histogram (left Y-axis) ---
        TH1F* hAmp = new TH1F(Form("hAmp_evt%d_dut%d_type%d", eventID, DUT_ID, type),
                              "", STRIP_XMAX - STRIP_XMIN, STRIP_XMIN, STRIP_XMAX);
        for (const auto& [id, amp] : stripAmp)
            if (id >= STRIP_XMIN && id < STRIP_XMAX)
                hAmp->SetBinContent(id + 1, amp);

        hAmp->SetStats(0);
        hAmp->SetLineColor(kBlue);
        hAmp->SetFillColor(kBlue);
        hAmp->SetFillStyle(3004);
        hAmp->GetXaxis()->SetTitle("Strip ID");
        hAmp->GetYaxis()->SetTitle("Amplitude [ADC]");
        hAmp->GetYaxis()->SetTitleColor(kBlue);
        hAmp->GetYaxis()->SetLabelColor(kBlue);
        hAmp->Draw("HIST");

        // --- Time points (right Y-axis) ---
        // Create a TGaxis for right Y
        TGaxis* axis = new TGaxis(STRIP_XMAX, hAmp->GetMinimum(), STRIP_XMAX, hAmp->GetMaximum(),
                                  TIME_YMIN, TIME_YMAX, 510, "+L");
        axis->SetLineColor(kRed);
        axis->SetLabelColor(kRed);
        axis->SetTitle("Time [ns]");
        axis->SetTitleColor(kRed);
        axis->Draw();

        TGraph* gTime = new TGraph();
        for (const auto& [id, time] : stripTime)
            if (id >= STRIP_XMIN && id < STRIP_XMAX)
                gTime->SetPoint(gTime->GetN(), id, time * hAmp->GetMaximum() / TIME_YMAX);  // scale to amplitude

        gTime->SetMarkerColor(kRed);
        gTime->SetMarkerStyle(20);
        gTime->SetMarkerSize(0.8);
        gTime->Draw("P SAME");

        // --- Predicted hit marker ---
        double predHit = 0;
        if (type == 0)
            predHit = localHit.x() / det->getConfig().readoutPlanePitch.at(0);
        else if (type == 1)
            predHit = localHit.y() / det->getConfig().readoutPlanePitch.at(1);

        TMarker* mHit = new TMarker(predHit, hAmp->GetMaximum() * 0.5, 29);
        mHit->SetMarkerColor(kMagenta);
        mHit->SetMarkerSize(2.0);
        mHit->Draw();

        // --- Title ---
        TLatex latex;
        latex.SetTextSize(0.06);
        latex.SetTextColor(kBlack);
        latex.SetNDC();
        latex.DrawLatex(0.15, 0.94, Form("%s type=%d", detName.c_str(), type));
    }

    c->Update();
    c->Write(Form("canvas_evt%d_dut%d", eventID, DUT_ID));
}
// 为给定 cluster 绘制所有参与 strip 的波形（每条 strip 一条曲线）并保存到 dir 下
void EventDisplayManager::DrawClusterWaveforms(int eventID, std::shared_ptr<Detector> det, const Cluster& cluster, TDirectory* dir) {
    if (!det) return;
    const auto& raw = det->GetRawData();

    TCanvas* c = new TCanvas(Form("wave_cluster_evt%d_dut%d_clust%d", eventID, det->GetID(), cluster.matchID),
                             Form("Waveforms evt%d dut%d cluster%d", eventID, det->GetID(), cluster.matchID),
                             1200, 800);
    m_canvases.push_back(c);
    gPad->SetGrid(1, 1);
    int nStrips = (int)cluster.strips.size();
    int nCols = std::min(4, std::max(1, nStrips));
    int nRows = (nStrips + nCols - 1) / nCols;
    c->Divide(nCols, nRows);

    int pad = 1;
    for (const auto& sh : cluster.strips) {
        c->cd(pad++);
        gPad->SetGrid(1, 1);
        gPad->SetLeftMargin(0.12);
        gPad->SetRightMargin(0.05);

        // find rawdata
        const RawData* r = FindRawForStrip(raw, sh.stripID, sh.type);
        if (!r) {
            TLatex no;
            no.SetTextSize(0.04);
            no.DrawLatexNDC(0.2, 0.5, "No raw waveform");
            continue;
        }
        // draw all ADC points into TGraph
        TGraph* g = new TGraph();
        for (size_t i = 0; i < r->adc.size(); ++i) g->SetPoint((int)i, (double)i, (double)r->adc[i]);
        g->SetLineWidth(2);
        // color by isValid
        if (sh.isValid)
            g->SetLineColor(kBlue);
        else
            g->SetLineColor(kGray + 1);
        g->Draw("AL");
        g->GetXaxis()->SetTitle("Sample");
        g->GetYaxis()->SetTitle("ADC");
        // mark peak
        if (!r->adc.empty()) {
            auto it = std::max_element(r->adc.begin(), r->adc.end());
            int idx = (int)std::distance(r->adc.begin(), it);
            double val = *it;
            TMarker* m = new TMarker(idx, val, 29);
            m->SetMarkerColor(kRed);
            m->SetMarkerSize(1.6);
            m->Draw();
        }
        // write graph to dir
        dir->cd();
        g->Write(Form("wave_evt%d_dut%d_strip%d", eventID, det->GetID(), sh.stripID));
    }

    dir->cd();
    c->Write(Form("wave_canvas_evt%d_dut%d_cluster%d", eventID, det->GetID(), cluster.matchID));
}

// 在 det->GetRawData() 中查找对应 strip 的 RawData（按 stripID 和 type 匹配）
const RawData* EventDisplayManager::FindRawForStrip(const std::vector<RawData>& raw, int stripID, int type) const {
    for (const auto& r : raw) {
        if (r.stripID == stripID && r.type == type) return &r;
    }
    return nullptr;
}