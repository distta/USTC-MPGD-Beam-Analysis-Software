#pragma once

#include <nlohmann/json.hpp>
#include <string>

using json = nlohmann::json;

class IRawDataConverter {
   public:
    virtual ~IRawDataConverter() = default;
    virtual bool Convert(const json& config, const std::string& outputPath) = 0;
};
