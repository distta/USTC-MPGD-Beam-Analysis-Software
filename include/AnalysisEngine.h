

#pragma once

#include <map>
#include <memory>
#include <string>

#include <TFile.h>

#include "DataModel.h"
#include "Detector/Detector.h"
#include "Detector/Planar.h"
#include <nlohmann/json.hpp>

using json = nlohmann::json;

class AnalysisEngine {
   public:
    explicit AnalysisEngine(const std::string& configFile) : m_configFile(configFile) {
    }
    ~AnalysisEngine() = default;

    void Initialize();

    void Run(const std::string& rawFile, const std::string& cacheFile, const std::string& outFile);

   private:
    void RunClustering(const std::string& rawFile, const std::string& cacheFile);

    void RunAnalysis(const std::string& cacheFile, const std::string& outFile);

    void RunTrackAnalysis(TFile* outFile);

    void RunDUTAnalysis(TFile* outFile);

   private:
    // 把硬件 (boardID, channelID) 映射到 (detID, stripID, type)
    std::tuple<int, int, int> MapBoardChannel(unsigned int boardID, unsigned int channelID, unsigned int mm_strip) const;

    bool EventFilter();

    Track FitTrack(const std::vector<GlobalHit>& globalHits) const;

    // Compute track chi2
    double TrackChi2(const double* par);

    // Compute DUT chi2
    double DUTChi2(const double* par, int detID);

    // profile-likelihood
    double ProfileNLL(const double* par);

   private:
    std::string m_configFile;
    std::string m_rawDataFileName;

    json m_config;

    // detID -> shared_ptr<Detector>
    std::map<int, std::shared_ptr<Detector>> m_dets;
    std::vector<int> m_trackerIDs;
    std::vector<Event> m_events;
};
