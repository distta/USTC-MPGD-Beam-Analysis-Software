#include "Pipeline.h"
#include "DataModel.h"

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>

#include <TFile.h>
#include <TTree.h>

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>

using namespace ROOT;
using namespace std;

Pipeline::Pipeline(const std::string& configFile)
    : m_configFile(configFile) {
    std::ifstream in(configFile);
    if (!in.is_open()) {
        throw std::runtime_error("Pipeline: cannot open config file: " + configFile);
    }
    in >> m_config;
    InitializeDetectors();
}

void Pipeline::InitializeDetectors() {
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
        } else {
            throw std::runtime_error("Pipeline: unsupported detector type: " + type);
        }
    }
}

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

void Pipeline::Run(const std::string& rawFile, const std::string& cacheFile, const std::string& outFile) {
    auto t0 = std::chrono::high_resolution_clock::now();

    // enable ROOT implicit MT (auto parallel)
    // ROOT::EnableImplicitMT();

    if (!std::filesystem::exists(cacheFile)) {
        std::cout << "[Pipeline] Cache not found. Running Clustering stage to produce: " << cacheFile << std::endl;
        RunClustering(rawFile, cacheFile);
    } else {
        std::cout << "Cache file '" << cacheFile << "' already exists. Overwrite? (y/n): ";
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
    RunTracking(cacheFile, outFile);

    auto t1 = std::chrono::high_resolution_clock::now();
    double sec = std::chrono::duration<double>(t1 - t0).count();
    std::cout << "[Pipeline] Done. Total time: " << sec << " s\n";
}

void Pipeline::RunClustering(const std::string& rawFile, const std::string& cacheFile) {
    auto start = std::chrono::high_resolution_clock::now();
    ROOT::RDataFrame df("raw", rawFile);

    // --------------------------
    // Step 1: 构建 RawData
   auto df_result = df.Define("ClusterData", 
        [this](const std::vector<unsigned int>& apv_id, 
               const std::vector<unsigned int>& apv_ch, 
               std::vector<std::vector<short>>& apv_q,
               unsigned int apv_evt) {
            
            // 构建RawData
            std::unordered_map<int, std::vector<RawData>> rawHits;
            for (size_t i = 0; i < apv_id.size(); ++i) {
                auto [detID, stripID, type] = MapBoardChannel(apv_id[i], apv_ch[i]);
                if (m_dets.find(detID) != m_dets.end()) {
                    rawHits[detID].emplace_back(RawData{stripID, type, apv_q[i]});
                }
            }
            
             
            int eventID=apv_evt;
            std::vector<int> detIDs;  
            std::vector<RecCluster> clusters;
            
            for (auto& [dID, raws] : rawHits) {
                auto recHits = m_dets[dID]->Reconstruction(raws);
                
                for (const auto& hit : recHits) {
                    for (const auto& cluster : hit.cluster) {
                        detIDs.push_back(dID);
                        clusters.push_back(cluster);
                    }
                }
            }
            
            return std::make_tuple(eventID, detIDs, clusters);
        },
        {"apv_id", "apv_ch", "apv_q", "apv_evt"})
        
        .Define("eventID", [](const std::tuple<int, std::vector<int>, std::vector<RecCluster>>& data) {
            return std::get<0>(data);
        }, {"ClusterData"})
        .Define("detID", [](const std::tuple<int, std::vector<int>, std::vector<RecCluster>>& data) {
            return std::get<1>(data);
        }, {"ClusterData"}) 
        .Define("cluster", [](const std::tuple<int, std::vector<int>, std::vector<RecCluster>>& data) {
            return std::get<2>(data);
        }, {"ClusterData"});
    
    std::cout << "[Pipeline] Writing cache to " << cacheFile << " ..." << std::endl;
    
    df_result.Range(10000).Snapshot("clusters", cacheFile, {"eventID", "detID", "cluster"});
    
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    
    std::cout << "[Pipeline] Cache written successfully. Cost Time: " 
              << duration.count() << " s" << std::endl;
}

void Pipeline::RunTracking(const std::string& cacheFile, const std::string& outFile) {
}

Track Pipeline::FitTrack(const std::map<int, std::vector<GlobalHit>>& globalHits) const {
    Track t;
    std::vector<double> zs, xs, ys;
    // collect first hit per detector (if exist). You may change to use centroid or best hit.
    for (const auto& [detID, vec] : globalHits) {
        if (vec.empty()) continue;
        zs.push_back(vec.front().z);
        xs.push_back(vec.front().x);
        ys.push_back(vec.front().y);
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

bool Pipeline::EventFilter(const std::vector<RecHit>& recHits) {

    std::map<int, int> detectorHitCount;
    for (const auto& hit : recHits) {
        detectorHitCount[hit.detID]++;
    }

    for (const auto& pair : detectorHitCount) {
        auto it = m_dets.find(pair.first);
        if (it == m_dets.end()) {
            std::cerr << "[Error] Unknown detID=" << pair.first << std::endl;
            return false;
        }
        if (it->second->isDUT()) continue;
        if (pair.second != 1) {
            return false;
        }
    }

    return true;
}