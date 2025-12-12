#include "Detector/DetectorFactory.h"
#include "Detector/Cylinder.h"
#include "Detector/Planar.h"

#include <iostream>
#include <stdexcept>

using namespace std;

DetectorFactory& DetectorFactory::GetInstance() {
    static DetectorFactory instance;
    return instance;
}

bool DetectorFactory::Initialize(const json& config) {

    if (!config.contains("detectors")) {
        cerr << "[DetectorFactory] Error: No 'detectors' field in config" << endl;
        return false;
    }

    cout << "[DetectorFactory] Initializing detectors..." << endl;

    try {
        for (const auto& detConfig : config["detectors"]) {
            auto detector = CreateDetector(detConfig);
            if (detector) {
                int id = detector->GetID();
                if (m_detectors.count(id)) {
                    throw runtime_error("Duplicate detector ID: " + to_string(id));
                }
                m_detectors[id] = detector;
                cout << "[DetectorFactory] Created detector: " << detector->GetName()
                     << " (ID=" << id << ", Type=" << detConfig.value("type", "planar") << ")" << endl;
            }
        }

        cout << "[DetectorFactory] Successfully initialized " << m_detectors.size() << " detectors" << endl;
        return true;

    } catch (const exception& e) {
        cerr << "[DetectorFactory] Initialization failed: " << e.what() << endl;
        Clear();
        return false;
    }
}

shared_ptr<Detector> DetectorFactory::CreateDetector(const json& detConfig) {
    if (!detConfig.contains("id")) {
        throw runtime_error("Detector config missing 'id' field");
    }
    if (!detConfig.contains("name")) {
        throw runtime_error("Detector config missing 'name' field");
    }

    int id = detConfig["id"];
    string name = detConfig["name"];
    string type = detConfig.value("type", "planar");

    shared_ptr<Detector> detector;

    if (type == "planar") {
        detector = make_shared<Planar>(id, name, detConfig);
    } else if (type == "cylinder") {
        detector = make_shared<Cylinder>(id, name, detConfig);
    } else {
        throw runtime_error("Unsupported detector type: " + type);
    }

    return detector;
}

shared_ptr<Detector> DetectorFactory::GetDetector(int id) const {
    auto it = m_detectors.find(id);
    if (it != m_detectors.end()) {
        return it->second;
    }
    return nullptr;
}

const map<int, shared_ptr<Detector>>& DetectorFactory::GetAllDetectors() const {
    return m_detectors;
}

vector<shared_ptr<Detector>> DetectorFactory::GetDetectorsByRole(Detector::Role role) const {
    vector<shared_ptr<Detector>> result;
    for (const auto& [id, det] : m_detectors) {
        if ((role == Detector::Role::Tracker && det->isTracker()) ||
            (role == Detector::Role::DUT && det->isDUT()) ||
            (role == Detector::Role::Ignored && !det->isTracker() && !det->isDUT())) {
            result.push_back(det);
        }
    }
    return result;
}

vector<int> DetectorFactory::GetDetectorIDsByRole(Detector::Role role) const {
    vector<int> result;
    for (const auto& [id, det] : m_detectors) {
        if ((role == Detector::Role::Tracker && det->isTracker()) ||
            (role == Detector::Role::DUT && det->isDUT()) ||
            (role == Detector::Role::Ignored && !det->isTracker() && !det->isDUT())) {
            result.push_back(id);
        }
    }
    return result;
}

void DetectorFactory::Clear() {
    m_detectors.clear();
    cout << "[DetectorFactory] Cleared all detectors" << endl;
}
