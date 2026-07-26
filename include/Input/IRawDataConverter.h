#pragma once

#include <filesystem>
#include <nlohmann/json.hpp>
#include <string>

class IRawDataConverter {
   public:
    virtual ~IRawDataConverter() = default;
    virtual void Configure(const nlohmann::json& config) {
        (void)config;
    }
    virtual bool AcquireRawData(const std::filesystem::path& rawDir,
                                const std::string& runID,
                                std::string& error) = 0;
    virtual bool Convert(const std::string& outputPath) = 0;
};
