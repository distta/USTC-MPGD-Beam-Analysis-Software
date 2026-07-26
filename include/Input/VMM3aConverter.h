#pragma once

#include "Input/IRawDataConverter.h"

#include <filesystem>
#include <nlohmann/json.hpp>
#include <string>

class VMM3aConverter : public IRawDataConverter {
   public:
    void Configure(const nlohmann::json& config) override {
        m_config = config;
    }

    bool AcquireRawData(const std::filesystem::path& rawDir,
                        const std::string& runID,
                        std::string& error) override;
    bool Convert(const std::string& outputPath) override;

   private:
    nlohmann::json m_config;
    std::filesystem::path m_inputPath;
    std::filesystem::path m_geometryPath;
};
