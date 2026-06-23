#pragma once

#include "Input/IRawDataConverter.h"
#include <memory>
#include <string>

class ConverterFactory {
   public:
    static std::unique_ptr<IRawDataConverter> Create(const std::string& type);
};
