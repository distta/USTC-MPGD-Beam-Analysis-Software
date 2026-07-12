#include "Detector/DetectorFactory.h"
#include "Detector/Cylinder.h"
#include "Detector/Planar.h"
#include "Detector/PlanarPad.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;

namespace {

string RoleName(const Detector& detector) {
    if (detector.isTracker()) return "Tracker";
    if (detector.isDUT()) return "DUT";
    return "Ignored";
}

string TypeName(const shared_ptr<Detector>& detector) {
    if (dynamic_pointer_cast<PlanarPad>(detector)) return "planar_pad";
    if (dynamic_pointer_cast<Planar>(detector)) return "planar";
    if (dynamic_pointer_cast<Cylinder>(detector)) return "cylinder";
    return "unknown";
}

string FormatVector(const TVector3& value) {
    ostringstream output;
    output << fixed << setprecision(3) << '(' << value.X() << ", "
           << value.Y() << ", " << value.Z() << ')';
    return output.str();
}

void PrintGeometry(const map<int, shared_ptr<Detector>>& detectors) {
    size_t nameWidth = 4;
    size_t trackerCount = 0;
    size_t dutCount = 0;
    size_t ignoredCount = 0;

    for (const auto& [id, detector] : detectors) {
        (void)id;
        nameWidth = max(nameWidth, detector->GetName().size());
        if (detector->isTracker()) {
            ++trackerCount;
        } else if (detector->isDUT()) {
            ++dutCount;
        } else {
            ++ignoredCount;
        }
    }

    const size_t blockWidth = max<size_t>(89, nameWidth + 83);
    const string border(blockWidth, '=');
    const string divider(blockWidth, '-');

    cout << '\n'
         << border << '\n'
         << "  Detector Geometry\n"
         << "  Initialized " << detectors.size() << " detectors: "
         << trackerCount << " Tracker, " << dutCount << " DUT, "
         << ignoredCount << " Ignored\n"
         << divider << '\n'
         << "  " << left << setw(5) << "ID" << setw(9) << "Role"
         << setw(9) << "Type"
         << setw(static_cast<int>(nameWidth + 2)) << "Name"
         << setw(31) << "Position (x, y, z)" << "Rotation (x, y, z)\n"
         << divider << '\n';

    for (const auto& [id, detector] : detectors) {
        cout << "  " << left << setw(5) << id << setw(9) << RoleName(*detector)
             << setw(9) << TypeName(detector)
             << setw(static_cast<int>(nameWidth + 2)) << detector->GetName()
             << setw(31) << FormatVector(detector->GetPos())
             << FormatVector(detector->GetRot()) << '\n';
    }
    cout << right << border << "\n\n" << flush;
}

}  // namespace

DetectorFactory& DetectorFactory::GetInstance() {
    static DetectorFactory instance;
    return instance;
}

bool DetectorFactory::Initialize(const json& config) {

    if (!config.contains("detectors")) {
        cerr << "[DetectorFactory] Error: No 'detectors' field in config" << endl;
        return false;
    }

    try {
        for (const auto& detConfig : config["detectors"]) {
            auto detector = CreateDetector(detConfig);
            if (detector) {
                int id = detector->GetID();
                if (m_detectors.count(id)) {
                    throw runtime_error("Duplicate detector ID: " + to_string(id));
                }
                m_detectors[id] = detector;
            }
        }

        PrintGeometry(m_detectors);
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
    } else if (type == "planar_pad") {
        detector = make_shared<PlanarPad>(id, name, detConfig);
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
}
