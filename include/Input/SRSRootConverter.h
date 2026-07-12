#pragma once

#include "Input/IRawDataConverter.h"

class SRSRootConverter : public IRawDataConverter {
   public:
    bool AcquireRawData(const std::filesystem::path& rawDir,
                        const std::string& runID,
                        std::string& error) override;
    bool Convert(const std::string& outputPath) override;

   private:
    std::filesystem::path m_inputPath;
};
