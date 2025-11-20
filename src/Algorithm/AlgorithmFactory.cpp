#include "Algorithm/AlgorithmFactory.h"
#include <iostream>
#include <sstream>
#include <stdexcept>

AlgorithmFactory& AlgorithmFactory::Instance() {
    static AlgorithmFactory instance;
    return instance;
}

void AlgorithmFactory::RegisterAlgorithm(const std::string& name, AlgorithmCreateFunc createFunc) {
    std::lock_guard<std::mutex> lock(m_mutex);

    if (m_registry.find(name) != m_registry.end()) {
        std::cerr << "[AlgorithmFactory] Warning: Algorithm '" << name << "' is already registered. Overwriting." << std::endl;
    }

    m_registry[name] = createFunc;
    std::cout << "[AlgorithmFactory] Registered algorithm: " << name << std::endl;
}

std::shared_ptr<IAlgorithm> AlgorithmFactory::CreateAlgorithm(const std::string& name, const json& config) {
    std::lock_guard<std::mutex> lock(m_mutex);

    auto it = m_registry.find(name);
    if (it == m_registry.end()) {
        // 提供有用的错误信息，列出所有可用的算法
        std::ostringstream oss;
        oss << "Algorithm '" << name << "' is not registered. Available algorithms: ";
        bool first = true;
        for (const auto& pair : m_registry) {
            if (!first) oss << ", ";
            oss << pair.first;
            first = false;
        }
        throw std::runtime_error(oss.str());
    }

    // 创建算法实例
    auto algorithm = it->second();

    // 加载配置
    try {
        algorithm->LoadConfig(config);
    } catch (const std::exception& e) {
        std::ostringstream oss;
        oss << "Failed to load config for algorithm '" << name << "': " << e.what();
        throw std::runtime_error(oss.str());
    }

    return algorithm;
}

bool AlgorithmFactory::IsRegistered(const std::string& name) const {
    std::lock_guard<std::mutex> lock(m_mutex);
    return m_registry.find(name) != m_registry.end();
}

std::vector<std::string> AlgorithmFactory::GetRegisteredNames() const {
    std::lock_guard<std::mutex> lock(m_mutex);

    std::vector<std::string> names;
    names.reserve(m_registry.size());

    for (const auto& pair : m_registry) {
        names.push_back(pair.first);
    }

    return names;
}
