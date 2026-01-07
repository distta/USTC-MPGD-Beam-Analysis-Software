#include "Script/ChannelNoiseScript.h"
#include "Script/Base/RawDataParser.h"
#include "Script/Base/ScriptFactory.h"
#include "Detector/DetectorFactory.h"
#include "Algorithm/AnalysisUtils.h"

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TStyle.h>

#include <algorithm>
#include <iostream>

void ChannelNoiseScript::LoadConfig(const json& config) {
    m_maxEvents = config.value("maxEvents", 0);  // 0 = 全部事例
    m_progressInterval = config.value("progressInterval", 1000);
    m_noisyThreshold = config.value("noisyThreshold", 10.0);  // ADC
    m_deadThreshold = config.value("deadThreshold", 1.0);     // ADC
}

void ChannelNoiseScript::Print() const {
    std::cout << "ChannelNoiseScript Configuration:" << std::endl;
    std::cout << "  Max Events: " << (m_maxEvents == 0 ? "All" : std::to_string(m_maxEvents)) << std::endl;
    std::cout << "  Progress Interval: " << m_progressInterval << std::endl;
    std::cout << "  Noisy Threshold: " << m_noisyThreshold << " ADC" << std::endl;
    std::cout << "  Dead Threshold: " << m_deadThreshold << " ADC" << std::endl;
}

bool ChannelNoiseScript::Execute() {
    std::cout << "\n========================================" << std::endl;
    std::cout << " Channel Noise Analysis Script" << std::endl;
    std::cout << "========================================" << std::endl;

    auto parser = GetParser();
    if (!parser) {
        std::cerr << "Error: Parser not set!" << std::endl;
        return false;
    }

    // 分析通道噪声
    AnalyzeChannelNoise();

    // 保存结果
    if (!SaveResults()) {
        std::cerr << "Error: Failed to save results!" << std::endl;
        return false;
    }

    std::cout << "\n========================================" << std::endl;
    std::cout << " Channel Noise Analysis Completed" << std::endl;
    std::cout << "========================================" << std::endl;

    return true;
}

void ChannelNoiseScript::AnalyzeChannelNoise() {
    auto parser = GetParser();
    Long64_t totalEvents = parser->GetTotalEvents();
    Long64_t eventsToProcess = (m_maxEvents > 0 && m_maxEvents < totalEvents) ? m_maxEvents : totalEvents;

    std::cout << "\nAnalyzing " << eventsToProcess << " events..." << std::endl;

    // 临时存储：[detID][type][stripID] -> baseline_list
    std::map<int, std::map<int, std::map<int, std::vector<double>>>> tempData;

    // 遍历所有事例
    for (Long64_t iEvt = 0; iEvt < eventsToProcess; ++iEvt) {
        if ((iEvt % m_progressInterval) == 0 || iEvt == eventsToProcess - 1) {
            std::cout << "\r[Channel Noise] Processed " << (iEvt + 1) << "/" << eventsToProcess << std::flush;
        }

        auto rawHitsMap = parser->LoadEvent(iEvt);
        if (rawHitsMap.empty())
            continue;

        // 遍历所有探测器
        for (const auto& [detID, rawHits] : rawHitsMap) {
            // 遍历该探测器的所有hit
            for (const auto& raw : rawHits) {
                int stripID = raw.stripID;
                int type = raw.type;
                const auto& waveform = raw.adc;

                if (waveform.empty())
                    continue;

                // 只取第一个采样点
                double baseline = static_cast<double>(waveform[0]);
                tempData[detID][type][stripID].push_back(baseline);
            }
        }
    }

    std::cout << "\n\nCalculating channel statistics..." << std::endl;

    // 计算每个通道的统计量
    m_channelNoiseMap.clear();
    int totalChannels = 0;
    int noisyChannels = 0;
    int deadChannels = 0;

    for (const auto& [detID, typeMap] : tempData) {
        for (const auto& [type, stripMap] : typeMap) {
            for (const auto& [stripID, baselineList] : stripMap) {
                if (baselineList.empty())
                    continue;

                ChannelNoise cn;
                cn.detectorID = detID;
                cn.stripID = stripID;
                cn.type = type;
                cn.nSamples = baselineList.size();
                cn.baselineList = baselineList;  // 保存原始数据用于绘制直方图

                // 使用 AnalysisUtils 计算平均基线和RMS
                cn.meanBaseline = AnalysisUtils::CalculateMean(baselineList);
                cn.rmsBaseline = AnalysisUtils::CalculateRMS(baselineList);

                // meanNoise 和 maxNoise 在这个实现中设为基线的RMS
                cn.meanNoise = cn.rmsBaseline;
                cn.maxNoise = cn.rmsBaseline;

                // 判断通道状态
                cn.isNoisy = (cn.rmsBaseline > m_noisyThreshold);
                cn.isDead = (cn.rmsBaseline < m_deadThreshold);

                m_channelNoiseMap[detID][type][stripID] = cn;

                totalChannels++;
                if (cn.isNoisy) noisyChannels++;
                if (cn.isDead) deadChannels++;
            }
        }
    }

    std::cout << "\nChannel Statistics:" << std::endl;
    std::cout << "  Total Channels: " << totalChannels << std::endl;
    std::cout << "  Noisy Channels: " << noisyChannels << " (" 
              << (100.0 * noisyChannels / totalChannels) << "%)" << std::endl;
    std::cout << "  Dead Channels: " << deadChannels << " (" 
              << (100.0 * deadChannels / totalChannels) << "%)" << std::endl;
}

bool ChannelNoiseScript::SaveResults() {
    std::string outputFile = GetOutputDir() + "ChannelNoise.root";
    std::cout << "\nSaving results to: " << outputFile << std::endl;

    TFile* file = new TFile(outputFile.c_str(), "RECREATE");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot create output file: " << outputFile << std::endl;
        return false;
    }

    // 创建TTree保存通道噪声统计
    TTree* tree = new TTree("ChannelNoise", "Channel Noise Statistics");
    
    Int_t detectorID, stripID, type, nSamples;
    Double_t meanBaseline, rmsBaseline, meanNoise, maxNoise;
    Bool_t isNoisy, isDead;

    tree->Branch("detectorID", &detectorID);
    tree->Branch("stripID", &stripID);
    tree->Branch("type", &type);
    tree->Branch("meanBaseline", &meanBaseline);
    tree->Branch("rmsBaseline", &rmsBaseline);
    tree->Branch("meanNoise", &meanNoise);
    tree->Branch("maxNoise", &maxNoise);
    tree->Branch("nSamples", &nSamples);
    tree->Branch("isNoisy", &isNoisy);
    tree->Branch("isDead", &isDead);

    // 填充TTree
    for (const auto& [detID, typeMap] : m_channelNoiseMap) {
        for (const auto& [tp, stripMap] : typeMap) {
            for (const auto& [strip, cn] : stripMap) {
                detectorID = cn.detectorID;
                stripID = cn.stripID;
                type = cn.type;
                meanBaseline = cn.meanBaseline;
                rmsBaseline = cn.rmsBaseline;
                meanNoise = cn.meanNoise;
                maxNoise = cn.maxNoise;
                nSamples = cn.nSamples;
                isNoisy = cn.isNoisy;
                isDead = cn.isDead;
                
                tree->Fill();
            }
        }
    }

    tree->Write();

    // 创建每个探测器的通道噪声分布直方图
    auto& factory = DetectorFactory::GetInstance();
    const auto& detectors = factory.GetAllDetectors();

    for (const auto& [detID, det] : detectors) {
        if (m_channelNoiseMap.find(detID) == m_channelNoiseMap.end())
            continue;

        file->mkdir(Form("Detector_%d", detID));
        file->cd(Form("Detector_%d", detID));

        // 为X和Y方向分别创建直方图
        for (int type = 0; type < 2; ++type) {
            if (m_channelNoiseMap[detID].find(type) == m_channelNoiseMap[detID].end())
                continue;

            const auto& stripMap = m_channelNoiseMap[detID][type];
            
            // 找出strip范围
            int minStrip = 9999, maxStrip = -1;
            for (const auto& [strip, cn] : stripMap) {
                if (strip < minStrip) minStrip = strip;
                if (strip > maxStrip) maxStrip = strip;
            }

            int nBins = maxStrip - minStrip + 1;
            std::string typeName = (type == 0) ? "X" : "Y";

            // 平均噪声 vs 通道号
            TH1D* hNoise = new TH1D(Form("hNoise_%s", typeName.c_str()),
                                    Form("Detector %d - %s Channels Mean Noise;Strip ID;Mean Noise (ADC)", detID, typeName.c_str()),
                                    nBins, minStrip - 0.5, maxStrip + 0.5);

            // 平均基线 vs 通道号
            TH1D* hBaseline = new TH1D(Form("hBaseline_%s", typeName.c_str()),
                                       Form("Detector %d - %s Channels Mean Baseline;Strip ID;Mean Baseline (ADC)", detID, typeName.c_str()),
                                       nBins, minStrip - 0.5, maxStrip + 0.5);

            // 基线RMS vs 通道号
            TH1D* hBaselineRMS = new TH1D(Form("hBaselineRMS_%s", typeName.c_str()),
                                          Form("Detector %d - %s Channels Baseline RMS;Strip ID;Baseline RMS (ADC)", detID, typeName.c_str()),
                                          nBins, minStrip - 0.5, maxStrip + 0.5);

            // 填充直方图
            for (const auto& [strip, cn] : stripMap) {
                hNoise->SetBinContent(strip - minStrip + 1, cn.meanNoise);
                hBaseline->SetBinContent(strip - minStrip + 1, cn.meanBaseline);
                hBaselineRMS->SetBinContent(strip - minStrip + 1, cn.rmsBaseline);
            }

            hNoise->Write();
            hBaseline->Write();
            hBaselineRMS->Write();

            delete hNoise;
            delete hBaseline;
            delete hBaselineRMS;
            
            // 为噪声通道创建详细的基线分布直方图
            file->mkdir(Form("Detector_%d/NoisyChannels_%s", detID, typeName.c_str()));
            file->cd(Form("Detector_%d/NoisyChannels_%s", detID, typeName.c_str()));
            
            for (const auto& [strip, cn] : stripMap) {
                if (cn.isNoisy) {
                    // 为每个噪声通道创建基线分布直方图
                    TH1D* hChannelBaseline = new TH1D(
                        Form("hBaseline_Det%d_%s_Strip%d", detID, typeName.c_str(), strip),
                        Form("Detector %d - %s Strip %d Baseline Distribution;Baseline (ADC);Entries", 
                             detID, typeName.c_str(), strip),
                        100,-100,100
                    );
                    hChannelBaseline->StatOverflows(kTRUE);
                    
                    // 填充原始基线数据
                    for (double baseline : cn.baselineList) {
                        hChannelBaseline->Fill(baseline);
                    }
                    std::cout<<AnalysisUtils::CalculateMean(cn.baselineList)<<" "<<hChannelBaseline->GetMean()<<std::endl;
                    std::cout << cn.baselineList.size() << " " << AnalysisUtils::CalculateRMS(cn.baselineList) <<" "<< hChannelBaseline->GetStdDev() << std::endl;
                    
                    std::cout << "\n  Created histogram for noisy channel: Det" << detID 
                              << " " << typeName << " Strip " << strip 
                              << " (RMS=" << cn.rmsBaseline << ")" << std::flush;

                    hChannelBaseline->Write();
                    delete hChannelBaseline;
                }
            }
        }

        file->cd();
    }

    // 创建全局噪声分布直方图
    file->cd();
    TH1D* hNoiseDistribution = new TH1D("hNoiseDistribution", 
                                        "Global Channel Noise Distribution;Mean Noise (ADC);Channels",
                                        100, 0, 20);
    
    TH1D* hBaselineDistribution = new TH1D("hBaselineDistribution",
                                           "Global Channel Baseline Distribution;Mean Baseline (ADC);Channels",
                                           100, 0, 500);

    for (const auto& [detID, typeMap] : m_channelNoiseMap) {
        for (const auto& [tp, stripMap] : typeMap) {
            for (const auto& [strip, cn] : stripMap) {
                hNoiseDistribution->Fill(cn.meanNoise);
                hBaselineDistribution->Fill(cn.meanBaseline);
            }
        }
    }

    hNoiseDistribution->Write();
    hBaselineDistribution->Write();

    delete hNoiseDistribution;
    delete hBaselineDistribution;

    file->Close();
    delete file;

    std::cout << "Results saved successfully!" << std::endl;
    return true;
}

// 注册脚本
REGISTER_SCRIPT("ChannelNoiseScript", ChannelNoiseScript);
