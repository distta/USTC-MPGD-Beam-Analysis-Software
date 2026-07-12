#include "Script/Base/ScriptFactory.h"
#include <stdexcept>
#include <iostream>

ScriptFactory& ScriptFactory::Instance() {
    static ScriptFactory instance;
    return instance;
}

void ScriptFactory::RegisterScript(const std::string& name, ScriptCreateFunc createFunc) {
    std::lock_guard<std::mutex> lock(m_mutex);
    
    if (m_registry.find(name) != m_registry.end()) {
        std::cerr << "Warning: Script '" << name << "' already registered. Overwriting." << std::endl;
    }
    
    m_registry[name] = createFunc;
}

std::shared_ptr<IScript> ScriptFactory::CreateScript(const std::string& name, const json& config) {
    std::lock_guard<std::mutex> lock(m_mutex);
    
    auto it = m_registry.find(name);
    if (it == m_registry.end()) {
        throw std::runtime_error("Script type '" + name + "' not registered");
    }
    
    auto script = it->second();
    if (script) {
        script->LoadConfig(config);
    }
    
    return script;
}

bool ScriptFactory::IsRegistered(const std::string& name) const {
    std::lock_guard<std::mutex> lock(m_mutex);
    return m_registry.find(name) != m_registry.end();
}

std::vector<std::string> ScriptFactory::GetRegisteredNames() const {
    std::lock_guard<std::mutex> lock(m_mutex);
    
    std::vector<std::string> names;
    names.reserve(m_registry.size());
    
    for (const auto& pair : m_registry) {
        names.push_back(pair.first);
    }
    
    return names;
}
