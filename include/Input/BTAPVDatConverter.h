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
    bool Convert(const json& config, const std::string& outputPath) override;

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
};
