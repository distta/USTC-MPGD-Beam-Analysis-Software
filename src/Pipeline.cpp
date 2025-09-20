#include "Pipeline.h"
#include "DataModel.h"
#include "Detector/Detector.h"
#include "Detector/Planar.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <nlohmann/json.hpp>
#include <ostream>

using json = nlohmann::json;
using namespace std;

Pipeline::~Pipeline() {
}

void Pipeline::Finalize() {
    std::cout << "Pipeline finalization completed." << std::endl;
    std::cout << "Total processed events: " << m_events.size() << std::endl;
}

void Pipeline::Initialize(const std::string& configFile) {
    std::ifstream in(configFile);
    if (!in.is_open()) {
        throw std::runtime_error("Failed to open config file: " + configFile);
    }

    json config;
    try {
        in >> config;
    } catch (const json::exception& e) {
        throw std::runtime_error("Failed to parse JSON config: " + std::string(e.what()));
    }

    if (!config.contains("detectors") || !config["detectors"].is_array() || config["detectors"].empty()) {
        throw std::runtime_error("No valid detectors found in configuration");
    }

    const auto& detectors = config["detectors"];
    for (const auto& detConfig : detectors) {

        if (!detConfig.contains("id") || !detConfig.contains("name")) {
            throw std::runtime_error("Detector missing required fields: id, name");
        }

        const int detID = detConfig["id"].get<int>();
        const std::string name = detConfig["name"].get<std::string>();

        if (m_dets.find(detID) != m_dets.end()) {
            throw std::runtime_error("Duplicate detector ID: " + std::to_string(detID));
        }

        const std::string geometryType = detConfig.value("type", "planar");

        std::shared_ptr<Detector> detector;
        if (geometryType == "planar") {
            detector = std::make_shared<Planar>(detID, name, detConfig);
        } else {
            throw std::runtime_error("Unsupported geometry type: " + geometryType);
        }

        AddDetector(detector);
    }
}

void Pipeline::SetRawDataFile(const std::string& dataFile) {
    m_rawDataFileName = dataFile;
}

void Pipeline::SetOutputDirectory(const std::string& outputDir) {
    m_outputDirectory = outputDir;
    // 更新 cache 文件路径以使用新的输出目录
    m_cacheFileName = m_outputDirectory + "/cache.root";
}

void Pipeline::GenerateCache() {

    if (std::filesystem::exists(m_cacheFileName)) {
        std::cout << "Cache file '" << m_cacheFileName << "' already exists. Overwrite? (y/n): ";
        char choice;
        std::cin >> choice;
        if (choice != 'y' && choice != 'Y') {
            std::cout << "Using existing cache file.\n";
            return;
        }
        std::cout << "Overwriting existing cache file.\n";
    }

    TFile outFile(m_cacheFileName.c_str(), "RECREATE");
    if (outFile.IsZombie()) {
        throw std::runtime_error("Failed to create cache ROOT file: " + m_cacheFileName);
    }

    TTree cacheTree("cluster", "cluster");
    std::vector<RecCluster> clusterBuffer;
    cacheTree.Branch("eventID", &m_eventID, "eventID/I");
    cacheTree.Branch("detID", &m_detID, "detID/I");
    cacheTree.Branch("cluster", &clusterBuffer);

    // 原始数据输入
    TChain rawChain("raw");
    if (rawChain.Add(m_rawDataFileName.c_str()) == 0) {
        throw std::runtime_error("Failed to load raw data file: " + m_rawDataFileName);
    }

    unsigned int apv_evt_ = 0;
    std::vector<unsigned int>* apv_id_ = nullptr;
    std::vector<unsigned int>* apv_ch_ = nullptr;
    std::vector<std::vector<short>>* apv_q_ = nullptr;

    rawChain.SetBranchAddress("apv_evt", &apv_evt_);
    rawChain.SetBranchAddress("apv_id", &apv_id_);
    rawChain.SetBranchAddress("apv_ch", &apv_ch_);
    rawChain.SetBranchAddress("apv_q", &apv_q_);

    int totalEntries = rawChain.GetEntries();
    std::map<int, std::vector<RawData>> rawDataBuffer;

    std::cout << "[INFO]: Raw File contains entries: " << totalEntries << std::endl;

    for (int i = 0; i < 10000; ++i) {
        rawChain.GetEntry(i);
        m_eventID = apv_evt_;

        for (size_t j = 0; j < apv_id_->size(); ++j) {
            auto [detID, stripID, type] = ElectronicMap((*apv_id_)[j], (*apv_ch_)[j]);
            rawDataBuffer[detID].push_back({stripID, type, (*apv_q_)[j]});
        }

        for (auto& [detID, raws] : rawDataBuffer) {
            if (m_dets.find(detID) == m_dets.end()) continue;
            clusterBuffer = m_dets[detID]->BuildClusters(raws);
            m_detID = detID;

            if (clusterBuffer.size() > 0)
                cacheTree.Fill();
        }

        rawDataBuffer.clear();

        if ((i % 100 == 0) || i == totalEntries - 1) {
            std::cout << "\r"
                      << "[INFO]: Processed " << (i + 1) << "/" << totalEntries << std::flush;
        }
    }
    std::cout << std::endl;
    cacheTree.Write();
    outFile.Close();
}

void Pipeline::Run() {

    auto start = std::chrono::high_resolution_clock::now();

    GenerateCache();

    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> duration = end - start;

    std::cout << "[INFO]: Generate Cache File execution time: "
              << duration.count() << " seconds" << std::endl;

    TFile* m_cacheFile = TFile::Open(m_cacheFileName.c_str(), "READ");
    if (!m_cacheFile || m_cacheFile->IsZombie()) {
        throw std::runtime_error("Failed to open cache file: " + m_cacheFileName);
    }

    TTree* m_cacheTree = (TTree*)m_cacheFile->Get("cluster");
    if (!m_cacheTree) {
        throw std::runtime_error("Cache file does not contain cluster tree!");
    }

    m_cacheTree->SetBranchAddress("eventID", &m_eventID);
    m_cacheTree->SetBranchAddress("detID", &m_detID);
    m_cacheTree->SetBranchAddress("cluster", &m_clusterBuffer);

    const int nEntries = m_cacheTree->GetEntries();
    if (nEntries == 0) {
        std::cerr << "[Warning]: No entries found in cache" << std::endl;
        return;
    }

    int currentEventID = -1;
    Event event;

    for (int i = 0; i < nEntries; ++i) {
        m_cacheTree->GetEntry(i);

        if (m_eventID != currentEventID) {
            if (currentEventID != -1) {
                event.eventID = currentEventID;

                CreateGlobalHits(event);

                if (EventFilter(event)) {
                    m_events.push_back(event);
                }
            }

            event = Event();
            currentEventID = m_eventID;
        }

        event.recClusters[m_detID] = *m_clusterBuffer;
    }

    std::cout << "Total valid events after filtering: " << m_events.size() << " / " << nEntries << std::endl;
}

bool Pipeline::EventFilter(Event& event) {
    std::vector<double> z_vals;
    std::vector<double> x_vals;
    std::vector<double> y_vals;

    // 1. 检查每个Tracker是否恰好只有一个GlobalHit
    for (auto& [detID, det] : m_dets) {
        if (!det->isTracker()) continue;

        auto it = event.recGlobalHits.find(detID);
        if (it == event.recGlobalHits.end() || it->second.size() != 1) {
            return false;
        }

        const GlobalHit& hit = it->second[0];
        z_vals.push_back(hit.z);
        x_vals.push_back(hit.x);
        y_vals.push_back(hit.y);
    }

    if (z_vals.size() < 2) {
        return false;
    }

    // 计算z、x、y的均值
    auto mean = [](const std::vector<double>& v) {
        double sum = 0.0;
        for (auto& val : v) sum += val;
        return sum / v.size();
    };

    double z_mean = mean(z_vals);
    double x_mean = mean(x_vals);
    double y_mean = mean(y_vals);

    // 计算斜率和截距
    auto calcSlopeIntercept = [&](const std::vector<double>& z, const std::vector<double>& val, double mean_z, double mean_val) {
        double num = 0.0;
        double den = 0.0;
        for (size_t i = 0; i < z.size(); ++i) {
            double dz = z[i] - mean_z;
            num += dz * (val[i] - mean_val);
            den += dz * dz;
        }
        double slope = num / den;
        double intercept = mean_val - slope * mean_z;
        return std::make_pair(slope, intercept);
    };

    auto [slope_x, intercept_x] = calcSlopeIntercept(z_vals, x_vals, z_mean, x_mean);
    auto [slope_y, intercept_y] = calcSlopeIntercept(z_vals, y_vals, z_mean, y_mean);

    // 计算chi2和ndf
    double chi2 = 0.0;
    for (size_t i = 0; i < z_vals.size(); ++i) {
        double dx = x_vals[i] - (intercept_x + slope_x * z_vals[i]);
        double dy = y_vals[i] - (intercept_y + slope_y * z_vals[i]);
        chi2 += dx * dx + dy * dy;
    }
    int ndf = static_cast<int>(2 * z_vals.size() - 4);

    // 5. 保存结果到event.track
    event.track.slope_x = slope_x;
    event.track.intercept_x = intercept_x;
    event.track.slope_y = slope_y;
    event.track.intercept_y = intercept_y;
    event.track.chi2 = chi2 / ndf;

    return true;
}

void Pipeline::CreateGlobalHits(Event& event) {
    for (auto& [detID, recClusters] : event.recClusters) {
        event.recLocalHits[detID] = m_dets[detID]->MatchCluster(recClusters);

        if (m_dets[detID]->isDUT()) continue;
        for (auto& hit : event.recLocalHits[detID]) {
            event.recGlobalHits[detID].push_back(m_dets[detID]->LocalToGlobal(hit));
        }
    }
}

std::tuple<int, int, int> Pipeline::ElectronicMap(int boardID, int channelID) {
    int rawDataIndex = -1;

    if (boardID == 0 || boardID == 1) {
        rawDataIndex = 0;
    } else if (boardID == 2 || boardID == 3) {
        rawDataIndex = 1;
    } else if (boardID == 4 || boardID == 5) {
        rawDataIndex = 2;
    } else if (boardID == 6 || boardID == 7) {
        rawDataIndex = 3;
    } else if (boardID == 12 || boardID == 13) {
        rawDataIndex = 4;
    } else if (boardID == 14 || boardID == 15) {
        rawDataIndex = 5;
    } else if (boardID == 9) {
        rawDataIndex = 6;
    }

    int type = rawDataIndex % 2 == 0 ? 0 : 1;
    int detID = int(rawDataIndex / 2);
    int stripID = channelID;

    return {detID, stripID, type};
}