
#include "RawDataParser.h"
#include "Detector/DetectorFactory.h"
#include <array>
#include <iostream>

RawDataParser::RawDataParser(const std::string& rawFile)
    : m_rawFile(rawFile) {
}

RawDataParser::~RawDataParser() {
    if (m_file) {
        m_file->Close();
        delete m_file;
        m_file = nullptr;
    }
}

bool RawDataParser::Initialize() {
    // 打开ROOT文件
    m_file = TFile::Open(m_rawFile.c_str(), "READ");
    if (!m_file || m_file->IsZombie()) {
        std::cerr << "[RawDataParser] Failed to open file: " << m_rawFile << std::endl;
        return false;
    }

    // 获取TTree
    m_tree = (TTree*)m_file->Get("raw");
    if (!m_tree) {
        std::cerr << "[RawDataParser] raw TTree not found in " << m_rawFile << std::endl;
        return false;
    }

    // 设置分支地址
    if (m_tree) {
        m_tree->SetBranchAddress("apv_id", &m_apv_id);
        m_tree->SetBranchAddress("apv_ch", &m_apv_ch);
        m_tree->SetBranchAddress("mm_strip", &m_mm_strip);
        m_tree->SetBranchAddress("apv_q", &m_apv_q);
        m_tree->SetBranchAddress("apv_evt", &m_apv_evt);
        m_numOfEvents = m_tree->GetEntries();
    }

    std::cout << "[RawDataParser] Initialized with " << m_file->GetName() << " containing " << m_numOfEvents << " events" << std::endl;
    return true;
}

bool RawDataParser::LoadEvent(int eventID) {
    if (!m_tree) {
        std::cerr << "[RawDataParser] TTree not initialized" << std::endl;
        return false;
    }

    if (eventID < 0 || eventID >= GetTotalEvents()) {
        std::cerr << "[RawDataParser] Invalid event index: " << eventID << " is out of range [0, " << GetTotalEvents() - 1 << "]" << std::endl;
        return false;
    }

    m_tree->GetEntry(eventID);

    auto& factory = DetectorFactory::GetInstance();
    auto& detectors = factory.GetAllDetectors();

    for (size_t j = 0; j < m_apv_id->size(); ++j) {
        auto [detID, stripID, type] = MapBoardChannel(
            (*m_apv_id)[j], (*m_apv_ch)[j], (*m_mm_strip)[j]);

        if (detectors.find(detID) == detectors.end()) continue;

        RawData raw{stripID, type, (*m_apv_q)[j]};
        detectors.at(detID)->AddRawData(raw);
    }

    return true;
}

// 硬件映射常量
constexpr std::array<int, 16> kBoardToRawIndex = {
    0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7};

std::tuple<int, int, int> RawDataParser::MapBoardChannel(unsigned int boardID, unsigned int channelID, unsigned int mm_strip) const {

    int rawDataIndex = (boardID < kBoardToRawIndex.size())
                           ? kBoardToRawIndex[boardID]
                           : static_cast<int>(boardID) / 2;

    int type = (rawDataIndex % 2 == 0) ? 0 : 1;
    int detID = (rawDataIndex / 2) + 1;
    int stripID = static_cast<int>(mm_strip);

    // 1726
    if (boardID == 12)
        stripID = channelID + 1;
    else if (boardID == 13)
        stripID = 129 + channelID;
    else if (boardID == 14)
        stripID = 256 - channelID;
    else if (boardID == 15)
        stripID = 128 - channelID;

    // 1774配置的特殊映射
    // if (boardID == 12)
    //     stripID = 384 - channelID;
    // else if (boardID == 13)
    //     stripID = 128 - channelID;
    // else if (boardID == 14)
    //     stripID = 129 + channelID;
    // else if (boardID == 15)
    //     stripID = 385 + channelID;

    // 1813
    // if (boardID == 12)
    //     stripID = 384 - channelID;
    // else if (boardID == 13)
    //     stripID = 128 - channelID;
    // else if (boardID == 14)
    //     stripID = 384 - channelID;
    // else if (boardID == 15)
    //     stripID = 128 - channelID;

    return {detID, stripID, type};
}
