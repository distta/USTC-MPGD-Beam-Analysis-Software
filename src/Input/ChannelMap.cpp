#include "Input/ChannelMap.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <unordered_map>

namespace {
std::string Trim(std::string value) {
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return {};
    const auto last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

std::vector<std::string> SplitCSV(const std::string& line) {
    std::vector<std::string> fields;
    std::stringstream stream(line);
    std::string field;
    while (std::getline(stream, field, ',')) fields.push_back(Trim(field));
    return fields;
}

std::string NormalizeHeader(std::string value) {
    value = Trim(value);
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return value;
}

bool LooksNumeric(const std::string& value) {
    if (value.empty()) return false;
    return std::isdigit(static_cast<unsigned char>(value[0])) || value[0] == '-' || value[0] == '+';
}
}  // namespace

uint64_t ChannelMap::MakeKey(const HardwareAddress& address) {
    return (static_cast<uint64_t>(address.boardID) & 0xffffULL) << 32 |
           (static_cast<uint64_t>(address.chipID) & 0xffffULL) << 16 |
           (static_cast<uint64_t>(address.channelID) & 0xffffULL);
}

bool ChannelMap::LoadCSV(const std::string& path) {
    std::ifstream input(path);
    if (!input) {
        std::cerr << "[ChannelMap] Cannot open " << path << '\n';
        return false;
    }

    m_entries.clear();
    m_index.clear();
    std::unordered_map<std::string, size_t> columns;
    bool hasNamedHeader = false;
    std::string line;
    size_t lineNumber = 0;
    while (std::getline(input, line)) {
        ++lineNumber;
        line = Trim(line);
        if (line.empty() || line[0] == '#') continue;

        const auto fields = SplitCSV(line);
        if (fields.empty()) continue;
        if (!LooksNumeric(fields[0])) {
            columns.clear();
            for (size_t i = 0; i < fields.size(); ++i) {
                columns[NormalizeHeader(fields[i])] = i;
            }
            hasNamedHeader = columns.count("board_id") &&
                             columns.count("chip_id") &&
                             columns.count("channel_id") &&
                             columns.count("detector_id") &&
                             columns.count("plane_type") &&
                             columns.count("strip_id");
            continue;
        }
        if (fields.size() < 6) {
            std::cerr << "[ChannelMap] Expected at least 6 columns at line " << lineNumber << '\n';
            return false;
        }

        try {
            auto column = [&](const std::string& name, size_t fallback) -> size_t {
                if (hasNamedHeader) {
                    const auto found = columns.find(name);
                    return found == columns.end() ? std::numeric_limits<size_t>::max() : found->second;
                }
                return fallback;
            };
            auto hasColumn = [&](const std::string& name, size_t fallback) {
                const size_t index = column(name, fallback);
                return index < fields.size() && !fields[index].empty();
            };
            auto value = [&](const std::string& name, size_t fallback) -> const std::string& {
                return fields[column(name, fallback)];
            };

            ChannelMapping entry;
            entry.hardware.boardID = std::stoi(value("board_id", 0));
            entry.hardware.chipID = std::stoi(value("chip_id", 1));
            entry.hardware.channelID = std::stoi(value("channel_id", 2));
            entry.detectorID = std::stoi(value("detector_id", 3));
            entry.planeType = std::stoi(value("plane_type", 4));
            entry.stripID = std::stoi(value("strip_id", 5));
            if (hasColumn("pedestal", 6)) entry.pedestal = std::stod(value("pedestal", 6));
            if (hasColumn("noise_sigma", 7)) entry.noiseSigma = std::stod(value("noise_sigma", 7));
            if (hasColumn("gain", 8)) entry.gain = std::stod(value("gain", 8));
            if (hasColumn("polarity", 9)) entry.polarity = std::stoi(value("polarity", 9));
            if (hasColumn("status", 10)) entry.status = static_cast<uint32_t>(std::stoul(value("status", 10)));

            const auto key = MakeKey(entry.hardware);
            if (m_index.count(key)) {
                std::cerr << "[ChannelMap] Duplicate hardware address at line " << lineNumber << '\n';
                return false;
            }
            m_index[key] = m_entries.size();
            m_entries.push_back(entry);
        } catch (const std::exception& error) {
            std::cerr << "[ChannelMap] Invalid value at line " << lineNumber << ": " << error.what() << '\n';
            return false;
        }
    }

    std::cout << "[ChannelMap] Loaded " << m_entries.size() << " channels from " << path << '\n';
    return !m_entries.empty();
}

const ChannelMapping* ChannelMap::Find(const HardwareAddress& address) const {
    const auto found = m_index.find(MakeKey(address));
    return found == m_index.end() ? nullptr : &m_entries[found->second];
}
