#include "DataModel.h"
#include "DetectorFrame.h"
#include "EventDisplayManager.h"

#include "TMultiGraph.h"
#include <TCanvas.h>
#include <TF1.h>
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

EventDisplayManager::EventDisplayManager(
    std::shared_ptr<RawDataParser> parser,
    const std::string& outputDirectory)
    : m_parser(std::move(parser)) {
    const auto output = std::filesystem::path(outputDirectory);
    m_trackFilePath = (output / "TrackInfo.root").string();
    m_outBaseDir = (output / "EventDisplay").string() + "/";
}

EventDisplayManager::~EventDisplayManager() {
    for (auto* c : m_canvases)
        if (c) delete c;
    m_canvases.clear();
}

// 初始化 RawDataParser 并加载 track entries
bool EventDisplayManager::Initialize() {
    if (!m_parser) {
        std::cerr << "[EventDisplay] RawDataParser 未初始化" << std::endl;
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

        std::vector<const TrackEntry*> matches;
        for (const auto& entry : m_trackEntries)
            if (entry.eventID == idx) matches.push_back(&entry);
        if (matches.empty()) {
            std::cerr << "未找到 event ID " << idx << std::endl;
            continue;
        }
        size_t trackIndex = 0;
        if (matches.size() > 1) {
            std::cout << "该事件有 " << matches.size() << " 条轨迹，请选择 [0-"
                      << matches.size() - 1 << "]: ";
            std::cin >> trackIndex;
            if (trackIndex >= matches.size()) {
                std::cerr << "无效的轨迹编号" << std::endl;
                continue;
            }
        }
        const auto& selected = *matches[trackIndex];
        if (!ProcessEntry(selected)) {
            std::cerr << "[EventDisplay] ProcessEntry 失败 for event " << selected.eventID << std::endl;
            continue;
        }
        std::cout << "[EventDisplay] 已生成 DUT overview for event " << selected.eventID << std::endl;

        // 波形查询循环
        QueryWaveforms(selected);
    }
}

void EventDisplayManager::QueryWaveforms(const TrackEntry& te) {
    auto& factory = DetectorFactory::GetInstance();
    auto eventRawHits = m_parser->LoadEvent(te.eventID);

    while (true) {
        std::cout << "\n查询波形: 输入 DUT ID, type, target stripID (空格分隔), 或 -1 退出: ";
        int dutID, type, targetStrip;
        std::cin >> dutID;
        if (dutID == -1) break;
        std::cin >> type >> targetStrip;

        auto det = factory.GetDetector(dutID);
        if (!det || !det->isDUT()) {
            std::cerr << "无效 DUT ID" << std::endl;
            continue;
        }

        auto rawIt = eventRawHits.find(dutID);
        if (rawIt == eventRawHits.end()) {
            std::cerr << "DUT " << dutID << " 没有数据" << std::endl;
            continue;
        }

        DrawWaveforms(te.eventID, dutID, type, targetStrip, rawIt->second);
    }
}

void EventDisplayManager::DrawWaveforms(int eventID, int dutID, int type, int targetStrip,
                                        const std::vector<RawData>& rawData) {
    // 查找目标strip及相邻strips
    std::vector<int> nearbyStrips = FindNearbyStrips(rawData, targetStrip, type);

    if (nearbyStrips.empty()) {
        std::cerr << "未找到目标 strip 及相邻 strip" << std::endl;
        return;
    }

    // 绘制波形
    TCanvas* c = new TCanvas(Form("wave_evt%d_dut%d_type%d", eventID, dutID, type),
                             Form("Event %d DUT %d type %d", eventID, dutID, type),
                             1000, 600);
    m_canvases.push_back(c);
    gPad->SetGrid(1, 1);

    TMultiGraph* mg = new TMultiGraph();
    const std::vector<int> colors = {kBlue, kGreen + 2, kOrange, kMagenta, kCyan, kGray + 1};
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);

    for (size_t i = 0; i < nearbyStrips.size(); ++i) {
        const RawData* r = FindRawForStrip(rawData, nearbyStrips[i], type);
        if (!r) continue;

        TGraph* g = new TGraph();
        for (size_t j = 0; j < r->adc.size(); ++j)
            g->SetPoint(j, (double)j, (double)r->adc[j]);

        g->SetLineWidth(2);
        g->SetLineColor(colors[i % colors.size()]);
        mg->Add(g, "L");
        legend->AddEntry(g, Form("strip %d", nearbyStrips[i]), "l");
    }

    mg->Draw("A");
    mg->SetTitle(Form("Event %d DUT %d Type %d;Sample;ADC", eventID, dutID, type));
    legend->Draw();

    // 保存输出
    std::string outDir = m_outBaseDir + "Event_" + std::to_string(eventID) + "/";
    std::filesystem::create_directories(outDir);
    c->SaveAs((outDir + Form("wave_nearby_dut%d_type%d_strip%d.png", dutID, type, targetStrip)).c_str());

    std::cout << "[EventDisplay] 波形已保存: target strip " << targetStrip << " 与相邻 strip\n";
}

std::vector<int> EventDisplayManager::FindNearbyStrips(const std::vector<RawData>& rawData,
                                                       int targetStrip, int type) const {
    std::vector<int> nearbyStrips;

    // 添加目标strip
    if (FindRawForStrip(rawData, targetStrip, type))
        nearbyStrips.push_back(targetStrip);

    // 向两侧查找,最多跳过1个缺失的strip
    auto findDirection = [&](int start, int step, int limit) {
        int missingCount = 0;
        for (int strip = start; (step > 0 ? strip < limit : strip >= limit); strip += step) {
            if (FindRawForStrip(rawData, strip, type)) {
                nearbyStrips.push_back(strip);
                missingCount = 0;
            } else if (++missingCount >= 2) {
                break;
            }
        }
    };

    findDirection(targetStrip - 1, -1, 0);   // 向左
    findDirection(targetStrip + 1, 1, 512);  // 向右

    std::sort(nearbyStrips.begin(), nearbyStrips.end());
    return nearbyStrips;
}

// 处理单个 entry：LoadEvent -> Reconstruct -> 绘图并保存
bool EventDisplayManager::ProcessEntry(const TrackEntry& te) {
    auto& factory = DetectorFactory::GetInstance();
    auto rawHits = m_parser->LoadEvent(te.eventID);
    if (rawHits.empty()) {
        std::cerr << "[EventDisplay] 无法加载 event " << te.eventID << " 的数据" << std::endl;
        return false;
    }

    // 准备输出ROOT文件
    std::string outDir = m_outBaseDir + "Event_" + std::to_string(te.eventID) + "/";
    std::filesystem::create_directories(outDir);
    std::string outFile = outDir + "Event_" + std::to_string(te.eventID) + ".root";

    TFile* fout = TFile::Open(outFile.c_str(), "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "[EventDisplay] 无法创建输出文件: " << outFile << std::endl;
        return false;
    }

    // 遍历所有DUT，创建 DetectorFrame 并绘图
    for (auto& [id, det] : factory.GetAllDetectors()) {
        if (!det->isDUT()) continue;

        auto detFrame = std::make_shared<DetectorFrame>(*det);
        auto it = rawHits.find(id);
        if (it != rawHits.end()) {
            detFrame->SetRawData(it->second);
            detFrame->Process();
        }

        DrawDUTOverview(te.eventID, det, detFrame, te.track);
    }

    fout->Write();
    fout->Close();
    delete fout;
    return true;
}

void EventDisplayManager::DrawDUTOverview(int eventID, std::shared_ptr<Detector> det,
                                          std::shared_ptr<DetectorFrame> detFrame, const Track& track) {
    const int DUT_ID = det->GetID();
    const std::string detName = det->GetName();

    // 从 DetectorFrame 获取 ChannelHits 并按 type 分组
    const auto& allChannelHits = detFrame->ChannelHits();
    const auto* planarConfig = det->GetPlanarConfig();
    const auto* padConfig = det->GetPlanarPadConfig();
    if (padConfig) {
        const double xMin = -0.5 * padConfig->pitchX;
        const double xMax =
            (padConfig->columns - 0.5) * padConfig->pitchX;
        const double yMin = -0.5 * padConfig->pitchY;
        const double yMax =
            (padConfig->rows - 0.5) * padConfig->pitchY;
        auto* canvas = new TCanvas(Form("DUT_%d_evt%d_pad_overview", DUT_ID, eventID),
                                   Form("DUT %d - event %d pad overview", DUT_ID, eventID),
                                   900, 800);
        m_canvases.push_back(canvas);
        auto* occupancy = new TH2F(Form("hPadAmp_evt%d_dut%d", eventID, DUT_ID),
                                   Form("%s;Local X [mm];Local Y [mm]", detName.c_str()),
                                   padConfig->columns, xMin, xMax,
                                   padConfig->rows, yMin, yMax);
        for (const auto& hit : allChannelHits) {
            if (!hit.isValid || hit.type != padConfig->planeType ||
                !hit.HasID1())
                continue;
            const int row = hit.id1;
            const int column = hit.id0;
            if (row < 0 || row >= padConfig->rows ||
                column < 0 || column >= padConfig->columns) continue;
            occupancy->SetBinContent(column + 1, row + 1, hit.amp);
        }
        occupancy->SetStats(0);
        occupancy->Draw("COLZ TEXT");
        const auto prediction = det->GlobalToLocal(det->CalcHitFromTrack(track));
        auto* marker = new TMarker(prediction.X(), prediction.Y(), 29);
        marker->SetMarkerColor(kMagenta);
        marker->SetMarkerSize(2.0);
        marker->Draw();
        canvas->Write(Form("canvas_evt%d_dut%d", eventID, DUT_ID));
        return;
    }
    if (!planarConfig) {
        std::cerr << "[EventDisplay] DUT overview currently supports planar strip detectors only; detector "
                  << det->GetName() << " uses another readout type" << std::endl;
        return;
    }
    std::map<int, std::vector<ChannelHit>> channelHitsMap;
    for (const auto& sh : allChannelHits) {
        channelHitsMap[sh.type].push_back(sh);
    }

    if (channelHitsMap.empty()) {
        TLatex note;
        note.SetTextSize(0.04);
        note.DrawLatexNDC(0.15, 0.5, Form("DUT %s (ID=%d) : no strip hits for event %d", detName.c_str(), DUT_ID, eventID));
        TH1F placeholder(Form("placeholder_evt%d_dut%d", eventID, DUT_ID), "no data", 1, 0, 1);
        placeholder.Write();
        return;
    }

    int nTypes = channelHitsMap.size();
    TCanvas* c = new TCanvas(Form("DUT_%d_evt%d_overview", DUT_ID, eventID),
                             Form("DUT %d - event %d overview", DUT_ID, eventID),
                             1200, std::max(500, nTypes * 300));
    m_canvases.push_back(c);
    c->Divide(1, nTypes, 0.005, 0.005);

    // 常量
    double TIME_YMIN = 0.0;
    double TIME_YMAX = 5.0;
    double TIME_SCALE = 0.021;  // ns

    // 计算预测击中位置
    GlobalHit globalHit = det->CalcHitFromTrack(track);
    TVector3 localPos = det->GlobalToLocal(globalHit);

    int padIdx = 1;
    for (const auto& [type, hits] : channelHitsMap) {

        int STRIP_XMIN = 0;
        int STRIP_XMAX = planarConfig->readoutPlaneStripNumber.at(type);

        c->cd(padIdx++);
        gPad->SetGrid(1, 1);

        // 收集strip幅度和时间
        std::map<int, double> stripAmp, stripTime;
        for (const auto& sh : hits) {
            stripAmp[sh.id0] = sh.amp;
            stripTime[sh.id0] = (sh.time - 81) * TIME_SCALE;
        }

        // 幅度直方图 (左Y轴)
        TH1F* hAmp = new TH1F(Form("hAmp_evt%d_dut%d_type%d", eventID, DUT_ID, type),
                              "", STRIP_XMAX - STRIP_XMIN, STRIP_XMIN, STRIP_XMAX);
        for (const auto& [id, amp] : stripAmp) {
            if (id >= STRIP_XMIN && id < STRIP_XMAX)
                hAmp->SetBinContent(id + 1, amp);
        }

        hAmp->SetStats(0);
        hAmp->SetLineColor(kBlue);
        hAmp->SetFillColor(kBlue);
        hAmp->SetFillStyle(3004);
        hAmp->GetXaxis()->SetTitle("Strip ID");
        hAmp->GetYaxis()->SetTitle("Amplitude [ADC]");
        hAmp->GetYaxis()->SetTitleColor(kBlue);
        hAmp->GetYaxis()->SetLabelColor(kBlue);
        hAmp->Draw("HIST");

        // 时间点 (右Y轴)
        TGaxis* axis = new TGaxis(STRIP_XMAX, hAmp->GetMinimum(), STRIP_XMAX, hAmp->GetMaximum(),
                                  TIME_YMIN, TIME_YMAX, 510, "+L");
        axis->SetLineColor(kRed);
        axis->SetLabelColor(kRed);
        axis->SetTitle("Time [ns]");
        axis->SetTitleColor(kRed);
        axis->Draw();

        TGraph* gTime = new TGraph();
        for (const auto& [id, time] : stripTime) {
            if (id >= STRIP_XMIN && id < STRIP_XMAX)
                gTime->SetPoint(gTime->GetN(), id, time * hAmp->GetMaximum() / TIME_YMAX);
        }
        gTime->SetMarkerColor(kRed);
        gTime->SetMarkerStyle(20);
        gTime->SetMarkerSize(0.8);
        gTime->Draw("P SAME");

        // 预测击中标记

        double predHit = (type == 0) ? localPos.X() / planarConfig->readoutPlanePitch.at(type)
                                     : localPos.Y() / planarConfig->readoutPlanePitch.at(type);

        TMarker* mHit = new TMarker(predHit, hAmp->GetMaximum() * 0.5, 29);
        mHit->SetMarkerColor(kMagenta);
        mHit->SetMarkerSize(2.0);
        mHit->Draw();

        // 标题
        TLatex latex;
        latex.SetTextSize(0.06);
        latex.SetTextColor(kBlack);
        latex.SetNDC();
        latex.DrawLatex(0.15, 0.94, Form("%s type=%d", detName.c_str(), type));
    }

    c->Update();
    c->Write(Form("canvas_evt%d_dut%d", eventID, DUT_ID));
}

// 在 det->GetRawData() 中查找对应 strip 的 RawData（按 id0 和 type 匹配）
const RawData* EventDisplayManager::FindRawForStrip(const std::vector<RawData>& raw, int stripID, int type) const {
    for (const auto& r : raw) {
        if (r.id0 == stripID && r.type == type) return &r;
    }
    return nullptr;
}
