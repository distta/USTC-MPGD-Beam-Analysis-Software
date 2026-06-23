#pragma once

#include "Input/IRawDataConverter.h"

class SRSRootConverter : public IRawDataConverter {
   public:
    bool Convert(const json& config, const std::string& outputPath) override;
};
