#include "Input/SRSRootConverter.h"

#include <TFile.h>
#include <TTree.h>

#include <array>
#include <filesystem>
#include <iostream>
#include <map>
#include <tuple>
#include <vector>

namespace fs = std::filesystem;

namespace {
struct SRSChannelInfo {
    int boardID{0};
    int chipID{0};
    int channelID{0};
    int detectorID{0};
    int planeType{0};
    int stripID{0};
    double noiseSigma{-1.0};
};

using ChannelKey = std::tuple<int, int, int>;

std::tuple<int, int, int> MapSRSChannel(unsigned int boardID,
                                        unsigned int channelID,
                                        unsigned int mmStrip) {
    constexpr std::array<int, 16> boardToRawIndex = {
        0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7};
    const int rawDataIndex = boardID < boardToRawIndex.size()
                                 ? boardToRawIndex[boardID]
                                 : static_cast<int>(boardID) / 2;
    const int planeType = rawDataIndex % 2 == 0 ? 0 : 1;
    const int detectorID = rawDataIndex / 2 + 1;
    int stripID = static_cast<int>(mmStrip);

    if (boardID == 12)
        stripID = 256 - static_cast<int>(channelID);
    else if (boardID == 13)
        stripID = 128 - static_cast<int>(channelID);

    if (boardID == 14)
        stripID = 256 - static_cast<int>(channelID);
    else if (boardID == 15)
        stripID = 128 - static_cast<int>(channelID);

    return {detectorID, stripID, planeType};
}
}  // namespace

bool SRSRootConverter::Convert(const json& config, const std::string& outputPath) {
    const std::string inputPath = config.value("input", "");
    if (inputPath.empty()) {
        std::cerr << "[SRSRootConverter] input is required\n";
        return false;
    }

    TFile inputFile(inputPath.c_str(), "READ");
    if (inputFile.IsZombie()) {
        std::cerr << "[SRSRootConverter] Cannot open " << inputPath << '\n';
        return false;
    }
    TTree* rawTree = static_cast<TTree*>(inputFile.Get("raw"));
    if (!rawTree) {
        std::cerr << "[SRSRootConverter] raw tree not found\n";
        return false;
    }

    const char* required[] = {"apv_evt", "apv_id", "apv_ch", "mm_strip", "apv_q"};
    for (const char* branch : required) {
        if (!rawTree->GetBranch(branch)) {
            std::cerr << "[SRSRootConverter] Missing raw branch: " << branch << '\n';
            return false;
        }
    }

    const bool hasTime = rawTree->GetBranch("time_s") && rawTree->GetBranch("time_us");
    unsigned int apvEvent = 0;
    int timeSeconds = 0, timeMicroseconds = 0;
    std::vector<unsigned int>* apvIDs = nullptr;
    std::vector<unsigned int>* apvChannels = nullptr;
    std::vector<unsigned int>* strips = nullptr;
    std::vector<std::vector<short>>* waveformsRaw = nullptr;

    rawTree->SetBranchAddress("apv_evt", &apvEvent);
    rawTree->SetBranchAddress("apv_id", &apvIDs);
    rawTree->SetBranchAddress("apv_ch", &apvChannels);
    rawTree->SetBranchAddress("mm_strip", &strips);
    rawTree->SetBranchAddress("apv_q", &waveformsRaw);
    if (hasTime) {
        rawTree->SetBranchAddress("time_s", &timeSeconds);
        rawTree->SetBranchAddress("time_us", &timeMicroseconds);
    }

    fs::path output(outputPath);
    if (output.has_parent_path()) fs::create_directories(output.parent_path());
    TFile outputFile(outputPath.c_str(), "RECREATE");
    if (outputFile.IsZombie()) {
        std::cerr << "[SRSRootConverter] Cannot create " << outputPath << '\n';
        return false;
    }

    ULong64_t eventID = 0, timestamp = 0;
    std::vector<int> detectorIDs, planeTypes, stripIDs;
    std::vector<int> boardIDs, chipIDs, channelIDs;
    std::vector<std::vector<short>> waveforms;
    std::vector<unsigned int> channelFlags;

    TTree eventTree("Events", "BeamAnalysis canonical raw events");
    eventTree.Branch("event_id", &eventID);
    eventTree.Branch("timestamp", &timestamp);
    eventTree.Branch("detector_id", &detectorIDs);
    eventTree.Branch("plane_type", &planeTypes);
    eventTree.Branch("strip_id", &stripIDs);
    eventTree.Branch("board_id", &boardIDs);
    eventTree.Branch("chip_id", &chipIDs);
    eventTree.Branch("channel_id", &channelIDs);
    eventTree.Branch("waveform", &waveforms);
    eventTree.Branch("channel_flags", &channelFlags);

    std::map<ChannelKey, SRSChannelInfo> observedChannels;
    const auto entries = rawTree->GetEntries();
    for (Long64_t entry = 0; entry < entries; ++entry) {
        rawTree->GetEntry(entry);
        const size_t count = apvIDs ? apvIDs->size() : 0;
        if (!apvChannels || !strips || !waveformsRaw ||
            apvChannels->size() != count || strips->size() != count || waveformsRaw->size() != count) {
            std::cerr << "[SRSRootConverter] Inconsistent vectors at entry " << entry << '\n';
            return false;
        }
        eventID = apvEvent;
        timestamp = hasTime
                        ? static_cast<ULong64_t>(timeSeconds) * 1000000ULL +
                              static_cast<ULong64_t>(timeMicroseconds)
                        : static_cast<ULong64_t>(entry);
        detectorIDs.clear();
        planeTypes.clear();
        stripIDs.clear();
        boardIDs.clear();
        chipIDs.clear();
        channelIDs.clear();
        waveforms.clear();
        channelFlags.clear();

        for (size_t i = 0; i < count; ++i) {
            const int boardID = static_cast<int>((*apvIDs)[i]);
            const int channelID = static_cast<int>((*apvChannels)[i]);
            auto [detectorID, stripID, planeType] =
                MapSRSChannel((*apvIDs)[i], (*apvChannels)[i], (*strips)[i]);
            detectorIDs.push_back(detectorID);
            planeTypes.push_back(planeType);
            stripIDs.push_back(stripID);
            boardIDs.push_back(boardID);
            chipIDs.push_back(0);
            channelIDs.push_back(channelID);
            waveforms.push_back((*waveformsRaw)[i]);
            channelFlags.push_back(0);
            observedChannels[{boardID, channelID, stripID}] =
                {boardID, 0, channelID, detectorID, planeType, stripID, -1.0};
        }
        eventTree.Fill();
    }

    TTree* pedestalTree = static_cast<TTree*>(inputFile.Get("pedestals"));
    if (pedestalTree && pedestalTree->GetEntries() > 0) {
        std::vector<unsigned int>* pedAPVIDs = nullptr;
        std::vector<unsigned int>* pedChannels = nullptr;
        std::vector<unsigned int>* pedStrips = nullptr;
        std::vector<double>* pedStd = nullptr;
        pedestalTree->SetBranchAddress("apv_id", &pedAPVIDs);
        pedestalTree->SetBranchAddress("apv_ch", &pedChannels);
        pedestalTree->SetBranchAddress("mm_strip", &pedStrips);
        pedestalTree->SetBranchAddress("apv_pedstd", &pedStd);
        pedestalTree->GetEntry(0);
        if (pedAPVIDs && pedChannels && pedStrips && pedStd) {
            for (size_t i = 0; i < pedAPVIDs->size(); ++i) {
                const int boardID = static_cast<int>((*pedAPVIDs)[i]);
                const int channelID = static_cast<int>((*pedChannels)[i]);
                auto [detectorID, stripID, planeType] =
                    MapSRSChannel((*pedAPVIDs)[i], (*pedChannels)[i], (*pedStrips)[i]);
                auto& info = observedChannels[{boardID, channelID, stripID}];
                info.boardID = boardID;
                info.channelID = channelID;
                info.stripID = stripID;
                info.noiseSigma = (*pedStd)[i];
                info.detectorID = detectorID;
                info.planeType = planeType;
            }
        }
    }

    int mapBoardID, mapChipID, mapChannelID;
    int detectorID, planeType, stripID, polarity = 1;
    double pedestal = 0.0, noiseSigma, gain = 1.0;
    unsigned int status = 0;
    TTree channelTree("Channels", "BeamAnalysis channel map and calibration");
    channelTree.Branch("board_id", &mapBoardID);
    channelTree.Branch("chip_id", &mapChipID);
    channelTree.Branch("channel_id", &mapChannelID);
    channelTree.Branch("detector_id", &detectorID);
    channelTree.Branch("plane_type", &planeType);
    channelTree.Branch("strip_id", &stripID);
    channelTree.Branch("pedestal", &pedestal);
    channelTree.Branch("noise_sigma", &noiseSigma);
    channelTree.Branch("gain", &gain);
    channelTree.Branch("polarity", &polarity);
    channelTree.Branch("status", &status);
    for (const auto& [key, info] : observedChannels) {
        mapBoardID = info.boardID;
        mapChipID = info.chipID;
        mapChannelID = info.channelID;
        detectorID = info.detectorID;
        planeType = info.planeType;
        stripID = info.stripID;
        noiseSigma = info.noiseSigma;
        channelTree.Fill();
    }

    int schemaVersion = 1;
    std::string sourceFormat = "SRS_ROOT";
    std::string converterVersion = "1.1.0";
    TTree metadataTree("Metadata", "BeamAnalysis canonical format metadata");
    metadataTree.Branch("schema_version", &schemaVersion);
    metadataTree.Branch("source_format", &sourceFormat);
    metadataTree.Branch("converter_version", &converterVersion);
    metadataTree.Fill();

    eventTree.Write();
    channelTree.Write();
    metadataTree.Write();
    outputFile.Close();
    inputFile.Close();
    std::cout << "[SRSRootConverter] Wrote " << entries << " events to " << outputPath << '\n';
    return true;
}
