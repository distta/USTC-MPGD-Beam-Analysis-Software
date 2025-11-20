#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "DataModel.h"
#include "Detector/Detector.h"
#include "RawDataParser.h"
#include <nlohmann/json.hpp>

using json = nlohmann::json;

Track FitTrack(const std::vector<GlobalHit>& globalHits);

class AnalysisEngine {
   public:
    AnalysisEngine(const std::string& configFile,
                   const std::string& rawDir,
                   const std::string& resultDir,
                   const std::string& runID);

    void Initialize();

    void RunTrackAnalysis();

    void RunDUTAnalysis();
    void RunDUTAlign();

   private:
    std::shared_ptr<RawDataParser> m_parser;

    std::string m_configFile;
    std::string m_rawDir;
    std::string m_resultDir;
    std::string m_runID;
    std::string m_outputDir;  // result/runID/

    json m_config;
    std::vector<Event> m_events;

    // DUT对齐相关私有方法
    double DUTChi2Objective(
        const double* par,
        const std::vector<Event>& events,
        std::shared_ptr<Detector> detector,
        int detID);

    bool CalcuDutResidual(
        std::shared_ptr<Detector> detector,
        const std::vector<LocalHit>& hits,
        const Track& track,
        double& hitX, double& hitY,
        double& residualX, double& residualY);
};
