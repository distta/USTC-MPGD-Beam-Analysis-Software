#pragma once

#include "DataModel.h"
#include "Detector/Detector.h"
#include <TChain.h>
#include <map>
#include <memory>
#include <string>

class Pipeline {
   public:
    Pipeline() = default;
    ~Pipeline();

    void Initialize(const std::string& configFile);
    void SetRawDataFile(const std::string& dataFile);
    void SetOutputDirectory(const std::string& outputDir);
    void Run();
    void Finalize();

   private:
    void AddDetector(std::shared_ptr<Detector> det) {
        m_dets[det->GetID()] = det;
    }
    std::tuple<int, int, int> ElectronicMap(int boardID, int channelID);
    bool EventFilter(const std::map<int, std::vector<std::vector<RawData>>>& preCluster);

    void GenerateCache();
    void CreateGlobalHits(Event& event);

   private:
    std::string m_cacheFileName;
    std::string m_rawDataFileName;
    std::string m_outputDirectory;

    int m_eventID;
    int m_detID;
    RecCluster* m_clusterBuffer;

    std::map<int, std::shared_ptr<Detector>> m_dets;
    std::vector<Event> m_events;
};