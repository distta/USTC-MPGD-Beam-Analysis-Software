#include "Input/BTAPVDatConverter.h"

#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace fs = std::filesystem;

namespace {

using HardwareKey = std::tuple<int, int, int>;

struct RunningStat {
    size_t count{0};
    double mean{0.0};
    double m2{0.0};

    void Fill(double value) {
        ++count;
        const double delta = value - mean;
        mean += delta / static_cast<double>(count);
        m2 += delta * (value - mean);
    }

    double Sigma() const {
        return count > 1
                   ? std::sqrt(m2 / static_cast<double>(count - 1))
                   : 0.0;
    }
};

std::vector<std::string> CollectDatFiles(const json& config,
                                         const char* inputKey,
                                         const char* filesKey) {
    std::vector<std::string> files;

    if (config.contains(filesKey)) {
        for (const auto& item : config[filesKey]) {
            files.push_back(item.get<std::string>());
        }
    }

    if (config.contains(inputKey)) {
        const fs::path input = config[inputKey].get<std::string>();

        if (fs::is_directory(input)) {
            for (const auto& entry : fs::directory_iterator(input)) {
                if (entry.is_regular_file() && entry.path().extension() == ".dat") {
                    files.push_back(entry.path().string());
                }
            }
        } else if (fs::is_regular_file(input)) {
            files.push_back(input.string());
        }
    }

    std::sort(files.begin(), files.end());
    files.erase(std::unique(files.begin(), files.end()), files.end());

    return files;
}

}  // namespace

bool BTAPVDatConverter::ReadWord(std::istream& input, uint16_t& value) {
    unsigned char bytes[2];

    if (!input.read(reinterpret_cast<char*>(bytes), sizeof(bytes))) {
        return false;
    }

    value = static_cast<uint16_t>(bytes[0]) << 8 | bytes[1];
    return true;
}

bool BTAPVDatConverter::DecodeFile(const std::string& path,
                                   int samplesPerWaveform,
                                   bool decoderDebug,
                                   size_t decoderDebugLimit,
                                   std::map<uint64_t, DecodedEvent>& events) const {
    std::ifstream input(path, std::ios::binary);

    if (!input) {
        std::cerr << "[BTAPVDatConverter] Cannot open " << path << '\n';
        return false;
    }

    uint16_t word = 0;
    bool synchronized = false;

    while (ReadWord(input, word)) {
        if (word == 0xcaca) {
            synchronized = true;
            break;
        }
    }

    if (!synchronized) {
        std::cerr << "[BTAPVDatConverter] Sync marker not found in "
                  << path << '\n';
        return false;
    }

    bool haveLocalEvent = false;
    uint16_t previousLocalEvent = 0;
    uint64_t localEpoch = 0;
    size_t blockCount = 0;

    while (ReadWord(input, word)) {
        if (word == 0xcaca) {
            continue;
        }

        if (word != 0xeeee) {
            continue;
        }

        uint16_t chipWord;
        uint16_t timeHigh;
        uint16_t timeLow;
        uint16_t timeTailWord;
        uint16_t eventWord;
        uint16_t channelCount;

        if (!ReadWord(input, chipWord) ||
            !ReadWord(input, timeHigh) ||
            !ReadWord(input, timeLow) ||
            !ReadWord(input, timeTailWord) ||
            !ReadWord(input, eventWord) ||
            !ReadWord(input, channelCount)) {
            break;
        }

        const int realChip = chipWord & 0xff;
        const int chipID = realChip % 4;
        const int mappedBoardID = realChip / 4;

        const uint64_t timestamp =
            (static_cast<uint64_t>(timeHigh) << 24) |
            (static_cast<uint64_t>(timeLow) << 8) |
            ((timeTailWord >> 8) & 0xff);

        // Extend the counter only for an actual 16-bit wrap. Small backwards
        // movements can occur when chip blocks arrive out of order and must
        // not create a spurious +65536 jump.
        if (haveLocalEvent && previousLocalEvent > 0xf000 && eventWord < 0x1000) {
            localEpoch += 65536ULL;
        }

        const uint64_t extendedEvent = localEpoch + eventWord;
        haveLocalEvent = true;
        previousLocalEvent = eventWord;

        uint16_t reservedWords[3];

        if (!ReadWord(input, reservedWords[0]) ||
            !ReadWord(input, reservedWords[1]) ||
            !ReadWord(input, reservedWords[2])) {
            break;
        }

        bool complete = true;
        std::vector<DecodedChannel> blockChannels;
        blockChannels.reserve(channelCount);

        for (uint16_t channelIndex = 0;
             channelIndex < channelCount;
             ++channelIndex) {
            uint16_t channelWord;

            if (!ReadWord(input, channelWord)) {
                complete = false;
                break;
            }

            DecodedChannel decoded;
            decoded.hardware = {
                mappedBoardID,
                chipID,
                channelWord & 0xff};

            decoded.waveform.reserve(samplesPerWaveform);

            for (int sample = 0; sample < samplesPerWaveform; ++sample) {
                uint16_t adcWord;

                if (!ReadWord(input, adcWord)) {
                    complete = false;
                    break;
                }

                short adc = static_cast<short>(adcWord & 0x0fff);

                if (adc >= 0x0800) {
                    adc = static_cast<short>(adc - 0x1000);
                }

                decoded.waveform.push_back(static_cast<short>(adc));
            }

            if (!complete) {
                break;
            }

            blockChannels.push_back(std::move(decoded));
        }

        if (!complete) {
            break;
        }

        if (decoderDebug &&
            (decoderDebugLimit == 0 || blockCount < decoderDebugLimit)) {
            const int firstChannel = blockChannels.empty()
                                         ? -1
                                         : blockChannels.front().hardware.channelID;
            const size_t waveformSamples = blockChannels.empty()
                                               ? 0
                                               : blockChannels.front().waveform.size();

            std::ostringstream debug;
            debug << "[BTAPVDatConverter][block] index=" << blockCount << '\n'
                  << " event_raw=" << eventWord
                  << " chip_word=0x" << std::hex << std::setw(4)
                  << std::setfill('0') << chipWord
                  << std::dec << std::setfill(' ')
                  << " real_chip=" << realChip
                  << " board=" << mappedBoardID
                  << " chip=" << chipID
                  << " timestamp=" << timestamp
                  << " channels_num=" << blockChannels.size()
                  << " waveform_samples=" << waveformSamples
                  << " time_words=0x" << std::hex << std::setw(4)
                  << std::setfill('0') << timeHigh
                  << ':' << std::setw(4) << timeLow
                  << ':' << std::setw(4) << timeTailWord;
            std::cout << debug.str() << '\n';
        } else if (decoderDebug && blockCount == decoderDebugLimit) {
            std::cout << "[BTAPVDatConverter][block] debug output limited to "
                      << decoderDebugLimit << " blocks for " << path << '\n';
        }

        // extendedEvent is only an internal map key that prevents two events
        // around a 16-bit wrap from being merged. The persisted event ID is
        // exactly the raw eventWord from the DAQ block.
        auto& event = events[extendedEvent];
        event.eventID = static_cast<uint64_t>(eventWord);
        if (!event.hasTimestamp) {
            event.timestamp = timestamp;
            event.hasTimestamp = true;
        } else if (decoderDebug && event.timestamp != timestamp &&
                   (decoderDebugLimit == 0 ||
                    blockCount < decoderDebugLimit)) {
            std::cout << "[BTAPVDatConverter][timestamp] warning: event="
                      << eventWord << " has timestamps "
                      << event.timestamp << " and " << timestamp << '\n';
        }
        event.channels.insert(event.channels.end(),
                              std::make_move_iterator(blockChannels.begin()),
                              std::make_move_iterator(blockChannels.end()));

        // The DAQ writes one padding byte after each block.
        input.ignore(1);
        ++blockCount;
    }

    std::cout << "[BTAPVDatConverter] " << path
              << ": " << blockCount << " blocks\n";

    return blockCount > 0;
}

bool BTAPVDatConverter::Convert(const json& config,
                                const std::string& outputPath) {
    const auto files = CollectDatFiles(config, "input", "files");
    const auto pedestalFiles = CollectDatFiles(
        config, "pedestal_input", "pedestal_files");

    if (files.empty()) {
        std::cerr << "[BTAPVDatConverter] No input .dat files\n";
        return false;
    }

    if (pedestalFiles.empty()) {
        std::cerr << "[BTAPVDatConverter] No pedestal .dat files; set "
                  << "pedestal_input or pedestal_files\n";
        return false;
    }

    ChannelMap channelMap;

    if (!config.contains("map") ||
        !channelMap.LoadCSV(config["map"].get<std::string>())) {
        return false;
    }

    const int samplesPerWaveform = config.value("samples_per_waveform", 30);
    // Keep the old timestamp_debug names as backward-compatible aliases.
    const bool decoderDebug = config.value("decoder_debug", false);
    const int configuredDecoderDebugLimit = config.value("decoder_debug_limit", 20);
    const size_t decoderDebugLimit = static_cast<size_t>(std::max(0, configuredDecoderDebugLimit));
    const bool allowUnmapped = config.value("allow_unmapped", false);

    std::map<HardwareKey, RunningStat> pedestalStats;

    for (const auto& file : pedestalFiles) {
        std::map<uint64_t, DecodedEvent> pedestalEvents;

        if (!DecodeFile(file,
                        samplesPerWaveform,
                        false,
                        0,
                        pedestalEvents)) {
            return false;
        }

        for (const auto& [eventKey, event] : pedestalEvents) {
            (void)eventKey;
            for (const auto& channel : event.channels) {
                const HardwareKey key{
                    channel.hardware.boardID,
                    channel.hardware.chipID,
                    channel.hardware.channelID};
                auto& stat = pedestalStats[key];
                for (const short sample : channel.waveform) {
                    stat.Fill(static_cast<double>(sample));
                }
            }
        }
    }

    size_t missingPedestalChannels = 0;
    for (const auto& mapping : channelMap.Entries()) {
        const HardwareKey key{
            mapping.hardware.boardID,
            mapping.hardware.chipID,
            mapping.hardware.channelID};
        const auto found = pedestalStats.find(key);
        if (found == pedestalStats.end() || found->second.count < 2) {
            ++missingPedestalChannels;
        }
    }

    if (missingPedestalChannels != 0) {
        std::cerr << "[BTAPVDatConverter] Pedestal data missing for "
                  << missingPedestalChannels << " mapped channels\n";
        return false;
    }

    fs::path output(outputPath);

    if (output.has_parent_path()) {
        fs::create_directories(output.parent_path());
    }

    TFile rootFile(outputPath.c_str(), "RECREATE");

    if (rootFile.IsZombie()) {
        std::cerr << "[BTAPVDatConverter] Cannot create "
                  << outputPath << '\n';
        return false;
    }

    ULong64_t eventID = 0;
    ULong64_t timestamp = 0;

    std::vector<int> detectorIDs;
    std::vector<int> planeTypes;
    std::vector<int> stripIDs;

    std::vector<int> boardIDs;
    std::vector<int> chipIDs;
    std::vector<int> channelIDs;

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

    size_t unmappedCount = 0;

    std::map<std::tuple<int, int, int>, size_t> unmappedAddresses;
    std::map<std::tuple<int, int>, size_t> unmappedGroups;

    size_t writtenEvents = 0;

    for (const auto& file : files) {
        std::map<uint64_t, DecodedEvent> fileEvents;

        if (!DecodeFile(file,
                        samplesPerWaveform,
                        decoderDebug,
                        decoderDebugLimit,
                        fileEvents)) {
            return false;
        }

        for (const auto& [id, decodedEvent] : fileEvents) {
            eventID = decodedEvent.eventID;
            timestamp = decodedEvent.timestamp;

            detectorIDs.clear();
            planeTypes.clear();
            stripIDs.clear();

            boardIDs.clear();
            chipIDs.clear();
            channelIDs.clear();

            waveforms.clear();
            channelFlags.clear();

            for (const auto& channel : decodedEvent.channels) {
                const auto* mapping = channelMap.Find(channel.hardware);

                if (!mapping) {
                    ++unmappedCount;

                    ++unmappedAddresses[{
                        channel.hardware.boardID,
                        channel.hardware.chipID,
                        channel.hardware.channelID}];

                    ++unmappedGroups[{
                        channel.hardware.boardID,
                        channel.hardware.chipID}];

                    detectorIDs.push_back(-1);
                    planeTypes.push_back(-1);
                    stripIDs.push_back(-1);

                    boardIDs.push_back(channel.hardware.boardID);
                    chipIDs.push_back(channel.hardware.chipID);
                    channelIDs.push_back(channel.hardware.channelID);

                    waveforms.push_back(channel.waveform);
                    channelFlags.push_back(0x80000000U);

                    continue;
                }

                detectorIDs.push_back(mapping->detectorID);
                planeTypes.push_back(mapping->planeType);
                stripIDs.push_back(mapping->stripID);

                boardIDs.push_back(channel.hardware.boardID);
                chipIDs.push_back(channel.hardware.chipID);
                channelIDs.push_back(channel.hardware.channelID);

                waveforms.push_back(channel.waveform);
                channelFlags.push_back(mapping->status);
            }

            eventTree.Fill();
            ++writtenEvents;
        }
    }

    int mapBoardID = 0;
    int mapChipID = 0;
    int mapChannelID = 0;

    int detectorID = 0;
    int planeType = 0;
    int stripID = 0;
    int polarity = 0;

    double pedestal = 0.0;
    double noiseSigma = 0.0;
    double gain = 1.0;

    unsigned int status = 0;

    TTree channelTree("Channels",
                      "BeamAnalysis channel map and calibration");

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

    for (const auto& mapping : channelMap.Entries()) {
        mapBoardID = mapping.hardware.boardID;
        mapChipID = mapping.hardware.chipID;
        mapChannelID = mapping.hardware.channelID;

        detectorID = mapping.detectorID;
        planeType = mapping.planeType;
        stripID = mapping.stripID;

        const HardwareKey key{
            mapping.hardware.boardID,
            mapping.hardware.chipID,
            mapping.hardware.channelID};
        const auto& stat = pedestalStats.at(key);
        pedestal = stat.mean;
        noiseSigma = stat.Sigma();

        gain = mapping.gain;
        polarity = mapping.polarity;
        status = mapping.status;

        channelTree.Fill();
    }

    int schemaVersion = 1;
    std::string sourceFormat = "BT_APV_DAT";
    std::string converterVersion = "2.0.0";
    std::string eventIDPolicy = "raw_event_word";
    std::string pedestalSource = "dedicated_dat";

    TTree metadataTree("Metadata",
                       "BeamAnalysis canonical format metadata");

    metadataTree.Branch("schema_version", &schemaVersion);
    metadataTree.Branch("source_format", &sourceFormat);
    metadataTree.Branch("converter_version", &converterVersion);
    metadataTree.Branch("event_id_policy", &eventIDPolicy);
    metadataTree.Branch("pedestal_source", &pedestalSource);

    metadataTree.Fill();

    eventTree.Write();
    channelTree.Write();
    metadataTree.Write();

    rootFile.Close();

    std::cout << "[BTAPVDatConverter] Wrote "
              << writtenEvents
              << " events to "
              << outputPath << '\n';

    if (unmappedCount > 0) {
        std::cerr << "[BTAPVDatConverter] Preserved "
                  << unmappedCount
                  << " unmapped channels with detector_id=-1\n";

        for (const auto& [group, count] : unmappedGroups) {
            const auto& [board, chip] = group;

            std::cerr << "  unmapped group board=" << board
                      << " chip=" << chip
                      << " occurrences=" << count << '\n';
        }

        size_t shown = 0;

        for (const auto& [address, count] : unmappedAddresses) {
            const auto& [board, chip, channel] = address;

            std::cerr << "  board=" << board
                      << " chip=" << chip
                      << " channel=" << channel
                      << " occurrences=" << count << '\n';

            if (++shown == 10) {
                break;
            }
        }

        if (!allowUnmapped) {
            return false;
        }
    }

    return true;
}
