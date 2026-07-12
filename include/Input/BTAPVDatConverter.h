#pragma once

#include "Input/ChannelMap.h"
#include "Input/IRawDataConverter.h"
#include <cstdint>
#include <fstream>
#include <map>
#include <string>
#include <vector>

class BTAPVDatConverter : public IRawDataConverter {
   public:
    bool AcquireRawData(const std::filesystem::path& rawDir,
                        const std::string& runID,
                        std::string& error) override;
    bool Convert(const std::string& outputPath) override;

   private:
    struct DecodedChannel {
        HardwareAddress hardware;
        std::vector<short> waveform;
    };

    struct DecodedEvent {
        uint64_t eventID{0};
        uint64_t timestamp{0};
        bool hasTimestamp{false};
        std::vector<DecodedChannel> channels;
    };

    bool DecodeFile(const std::string& path, int samplesPerWaveform,
                    bool decoderDebug,
                    size_t decoderDebugLimit,
                    std::map<uint64_t, DecodedEvent>& events) const;
    static bool ReadWord(std::istream& input, uint16_t& value);

    std::vector<std::string> m_files;
    std::vector<std::string> m_pedestalFiles;
    std::string m_mapPath;
};
