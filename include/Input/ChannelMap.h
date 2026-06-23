#pragma once

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

struct HardwareAddress {
    int boardID{0};
    int chipID{0};
    int channelID{0};
};

struct ChannelMapping {
    HardwareAddress hardware;
    int detectorID{0};
    int planeType{0};
    int stripID{0};
    double pedestal{0.0};
    double noiseSigma{-1.0};
    double gain{1.0};
    int polarity{1};
    uint32_t status{0};
};

class ChannelMap {
   public:
    bool LoadCSV(const std::string& path);
    const ChannelMapping* Find(const HardwareAddress& address) const;
    const std::vector<ChannelMapping>& Entries() const { return m_entries; }

   private:
    static uint64_t MakeKey(const HardwareAddress& address);

    std::vector<ChannelMapping> m_entries;
    std::unordered_map<uint64_t, size_t> m_index;
};
