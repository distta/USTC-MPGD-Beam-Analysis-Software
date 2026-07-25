#include "Input/ConverterFactory.h"
#include "Input/BTAPVDatConverter.h"
#include "Input/APV25SRSConverter.h"

#include <algorithm>
#include <cctype>

std::unique_ptr<IRawDataConverter> ConverterFactory::Create(const std::string& type) {
    std::string normalized = type;
    std::transform(normalized.begin(), normalized.end(), normalized.begin(), ::tolower);
    if (normalized == "btapvdat" || normalized == "bt-apv-dat" || normalized == "bt_apv_dat") {
        return std::make_unique<BTAPVDatConverter>();
    }
    if (normalized == "apv25srs") {
        return std::make_unique<APV25SRSConverter>();
    }
    return nullptr;
}
