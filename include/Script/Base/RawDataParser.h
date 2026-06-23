/**
 * @file RawDataParser.h
 * @brief BeamAnalysis标准ROOT读取器
 * @author Huang Qixuan
 */

#pragma once

#include "DataModel.h"
#include <TFile.h>
#include <TTree.h>
#include <map>
#include <string>
#include <vector>

class RawDataParser {
   public:
    RawDataParser(const std::string& rawFile);
    ~RawDataParser();

    bool Initialize();

    std::unordered_map<int, std::vector<RawData>> LoadEvent(int eventID);

    Long64_t GetTotalEvents() const { return m_numOfEvents; };
    ULong64_t GetCurrentEventID() const { return m_eventID; }

    // 获取指定通道的sigma，返回-1表示未找到
    double GetSigma(int detID, int type, int stripID) const;

    bool WriteDebugRoot(const std::string& outputFile);

   private:
    struct ChannelCalibration {
        double pedestal{0.0};
        double gain{1.0};
        int polarity{1};
    };

    std::string m_rawFile;

    // channel calibration: [detID][type][stripID]
    std::map<int, std::map<int, std::map<int, double>>> m_pedSigmaMap;
    std::map<int, std::map<int, std::map<int, ChannelCalibration>>> m_calibrationMap;

    bool LoadCanonicalChannelData();

    Long64_t m_numOfEvents{0};

    TFile* m_file{nullptr};
    TTree* m_tree{nullptr};

    ULong64_t m_eventID{0};
    ULong64_t m_timestamp{0};
    std::vector<int>* m_detectorIDs{nullptr};
    std::vector<int>* m_planeTypes{nullptr};
    std::vector<int>* m_stripIDs{nullptr};
    std::vector<std::vector<short>>* m_waveforms{nullptr};

};
