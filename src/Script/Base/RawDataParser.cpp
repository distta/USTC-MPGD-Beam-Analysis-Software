#include "Script/Base/RawDataParser.h"
#include "Detector/DetectorFactory.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TDirectory.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <utility>

RawDataParser::RawDataParser(const std::string& rawFile)
    : m_rawFile(rawFile) {
}

RawDataParser::~RawDataParser() {
    if (m_file) {
        m_file->Close();
        delete m_file;
        m_file = nullptr;
    }
}

bool RawDataParser::Initialize() {
    // 打开ROOT文件
    m_file = TFile::Open(m_rawFile.c_str(), "READ");
    if (!m_file || m_file->IsZombie()) {
        std::cerr << "[RawDataParser] Failed to open file: "
                  << m_rawFile << std::endl;
        return false;
    }

    m_tree = static_cast<TTree*>(m_file->Get("Events"));
    if (!m_tree) {
        std::cerr << "[RawDataParser] Events tree not found in "
                  << m_rawFile
                  << "; convert SRS/BT input to the standard ROOT format first"
                  << std::endl;
        return false;
    }

    TTree* metadataTree = static_cast<TTree*>(m_file->Get("Metadata"));
    if (metadataTree && metadataTree->GetBranch("schema_version")) {
        metadataTree->SetBranchAddress("schema_version", &m_schemaVersion);
        metadataTree->GetEntry(0);
        if (m_schemaVersion != 1) {
            std::cerr << "[RawDataParser] Unsupported schema version: "
                      << m_schemaVersion << std::endl;
            return false;
        }
    } else {
        m_schemaVersion = 1;
    }

    const char* required[] = {
        "event_id",
        "timestamp",
        "detector_id",
        "plane_type",
        "id",
        "waveform"
    };

    for (const char* branch : required) {
        if (!m_tree->GetBranch(branch)) {
            std::cerr << "[RawDataParser] Missing branch: "
                      << branch << std::endl;
            return false;
        }
    }

    m_tree->SetBranchAddress("event_id", &m_eventID);
    m_tree->SetBranchAddress("timestamp", &m_timestamp);
    m_tree->SetBranchAddress("detector_id", &m_detectorIDs);
    m_tree->SetBranchAddress("plane_type", &m_planeTypes);
    m_tree->SetBranchAddress("id", &m_stripIDs);
    m_tree->SetBranchAddress("waveform", &m_waveforms);

    m_numOfEvents = m_tree->GetEntries();

    LoadCanonicalChannelData();

    return true;
}

bool RawDataParser::LoadCanonicalChannelData() {
    TTree* channelTree = static_cast<TTree*>(m_file->Get("Channels"));
    if (!channelTree) {
        std::cerr << "[RawDataParser] Channels tree not found; "
                  << "noise sigma is unavailable"
                  << std::endl;
        return false;
    }

    int detectorID = 0;
    int planeType = 0;
    int stripID = 0;

    double pedestal = 0.0;
    double noiseSigma = -1.0;
    double gain = 1.0;

    int polarity = 1;

    channelTree->SetBranchAddress("detector_id", &detectorID);
    channelTree->SetBranchAddress("plane_type", &planeType);
    channelTree->SetBranchAddress("strip_id", &stripID);

    if (channelTree->GetBranch("pedestal")) {
        channelTree->SetBranchAddress("pedestal", &pedestal);
    }

    if (channelTree->GetBranch("noise_sigma")) {
        channelTree->SetBranchAddress("noise_sigma", &noiseSigma);
    }

    if (channelTree->GetBranch("gain")) {
        channelTree->SetBranchAddress("gain", &gain);
    }

    if (channelTree->GetBranch("polarity")) {
        channelTree->SetBranchAddress("polarity", &polarity);
    }

    const auto entries = channelTree->GetEntries();

    for (Long64_t i = 0; i < entries; ++i) {
        channelTree->GetEntry(i);

        if (noiseSigma >= 0) {
            m_pedSigmaMap[detectorID][planeType][stripID] = noiseSigma;
        }

        m_calibrationMap[detectorID][planeType][stripID] = {
            pedestal,
            gain,
            polarity
        };
    }

    return true;
}

double RawDataParser::GetSigma(int detID, int type, int stripID) const {
    auto it1 = m_pedSigmaMap.find(detID);
    if (it1 == m_pedSigmaMap.end()) return -1;

    auto it2 = it1->second.find(type);
    if (it2 == it1->second.end()) return -1;

    auto it3 = it2->second.find(stripID);
    if (it3 == it2->second.end()) return -1;

    return it3->second;
}

bool RawDataParser::WriteDebugRoot(const std::string& outputFile) {
    if (!m_tree) {
        std::cerr << "[RawDataParser] ERROR: TTree not initialized\n";
        return false;
    }

    TFile debugFile(outputFile.c_str(), "RECREATE");
    if (debugFile.IsZombie()) {
        std::cerr << "[RawDataParser] ERROR: cannot create debug file: "
                  << outputFile << '\n';
        return false;
    }

    const Long64_t nEvents = m_tree->GetEntries();

    std::cout << "[RawDataParser] Writing debug ROOT: "
              << outputFile << std::endl;

    for (Long64_t entry = 0; entry < nEvents; ++entry) {
        m_tree->GetEntry(entry);

        if (!m_detectorIDs || !m_planeTypes || !m_stripIDs || !m_waveforms) {
            continue;
        }

        const size_t nHits = m_detectorIDs->size();

        if (m_planeTypes->size() != nHits ||
            m_stripIDs->size() != nHits ||
            m_waveforms->size() != nHits) {
            std::cerr << "[RawDataParser] Inconsistent vectors at entry "
                      << entry << '\n';
            continue;
        }

        std::ostringstream dirName;
        dirName << "event_" << entry;

        TDirectory* eventDir = debugFile.mkdir(dirName.str().c_str());
        if (!eventDir) {
            eventDir = debugFile.GetDirectory(dirName.str().c_str());
        }

        if (!eventDir) {
            continue;
        }

        eventDir->cd();

        std::map<std::pair<int, int>, TH1D*> hMaxMap;
        std::map<std::pair<int, int>, TH2D*> hWaveMap;

        for (size_t i = 0; i < nHits; ++i) {
            const int detID = (*m_detectorIDs)[i];
            const int planeType = (*m_planeTypes)[i];
            const int stripID = (*m_stripIDs)[i];

            if (detID < 0) continue;
            if (planeType < 0) continue;
            if (stripID < 0 || stripID >= 256) continue;

            const auto& waveform = (*m_waveforms)[i];
            if (waveform.empty()) continue;

            const auto key = std::make_pair(detID, planeType);

            if (hMaxMap.find(key) == hMaxMap.end()) {
                std::ostringstream hMaxName;
                hMaxName << "h_maxA_event_" << entry
                         << "_det_" << detID
                         << "_plane_" << planeType;

                std::ostringstream hMaxTitle;
                hMaxTitle << "Event " << entry
                          << ", detID " << detID
                          << ", plane " << planeType
                          << ";stripID;maxA";

                hMaxMap[key] = new TH1D(
                    hMaxName.str().c_str(),
                    hMaxTitle.str().c_str(),
                    256,
                    -0.5,
                    255.5
                );

                const int nSamples = static_cast<int>(waveform.size());

                std::ostringstream hWaveName;
                hWaveName << "h_waveform_event_" << entry
                          << "_det_" << detID
                          << "_plane_" << planeType;

                std::ostringstream hWaveTitle;
                hWaveTitle << "Event " << entry
                           << ", detID " << detID
                           << ", plane " << planeType
                           << ";sample;stripID;ADC";

                hWaveMap[key] = new TH2D(
                    hWaveName.str().c_str(),
                    hWaveTitle.str().c_str(),
                    nSamples,
                    -0.5,
                    nSamples - 0.5,
                    256,
                    -0.5,
                    255.5
                );
            }

            double maxA = -std::numeric_limits<double>::max();

            const auto detIt = m_calibrationMap.find(detID);

            for (size_t s = 0; s < waveform.size(); ++s) {
                double value = static_cast<double>(waveform[s]);

                if (detIt != m_calibrationMap.end()) {
                    const auto planeIt = detIt->second.find(planeType);

                    if (planeIt != detIt->second.end()) {
                        const auto stripIt = planeIt->second.find(stripID);

                        if (stripIt != planeIt->second.end()) {
                            const auto& calibration = stripIt->second;

                            value = (value - calibration.pedestal)
                                    * calibration.gain
                                    * calibration.polarity;
                        }
                    }
                }

                if (value > maxA) {
                    maxA = value;
                }

                TH2D* hWave = hWaveMap[key];

                const int xbin = hWave->GetXaxis()->FindBin(
                    static_cast<double>(s)
                );
                const int ybin = hWave->GetYaxis()->FindBin(stripID);

                hWave->SetBinContent(xbin, ybin, value);
            }

            TH1D* hMax = hMaxMap[key];

            const int bin = hMax->GetXaxis()->FindBin(stripID);
            const double oldValue = hMax->GetBinContent(bin);

            if (maxA > oldValue) {
                hMax->SetBinContent(bin, maxA);
            }
        }

        for (auto& item : hMaxMap) {
            item.second->Write();
            delete item.second;
        }

        for (auto& item : hWaveMap) {
            item.second->Write();
            delete item.second;
        }

        debugFile.cd();

        if ((entry + 1) % 1000 == 0) {
            std::cout << "[RawDataParser] debug events written: "
                      << (entry + 1) << " / " << nEvents << std::endl;
        }
    }

    debugFile.Close();

    std::cout << "[RawDataParser] Debug ROOT written: "
              << outputFile << std::endl;

    return true;
}

std::unordered_map<int, std::vector<RawData>>
RawDataParser::LoadEvent(int eventID) {
    std::unordered_map<int, std::vector<RawData>> result;

    if (!m_tree) {
        std::cerr << "[RawDataParser] ERROR: TTree not initialized\n";
        return result;
    }

    if (eventID < 0 || eventID >= GetTotalEvents()) {
        std::cerr << "[RawDataParser] ERROR: Invalid event index "
                  << eventID << "\n";
        return result;
    }

    // ---- Load TTree Entry ----
    m_tree->GetEntry(eventID);

    auto& factory = DetectorFactory::GetInstance();
    const auto& detectors = factory.GetAllDetectors();

    if (!m_detectorIDs || !m_planeTypes || !m_stripIDs || !m_waveforms) {
        return result;
    }

    const size_t nHits = m_detectorIDs->size();

    if (m_planeTypes->size() != nHits ||
        m_stripIDs->size() != nHits ||
        m_waveforms->size() != nHits) {
        std::cerr << "[RawDataParser] Inconsistent vectors at entry "
                  << eventID << '\n';
        return result;
    }

    result.reserve(16);

    for (size_t i = 0; i < nHits; ++i) {
        const int detID = (*m_detectorIDs)[i];

        if (detectors.find(detID) == detectors.end()) {
            continue;
        }

        const int planeType = (*m_planeTypes)[i];
        const int stripID = (*m_stripIDs)[i];

        std::vector<short> calibrated = (*m_waveforms)[i];
        double gain = 1.0;

        const auto detIt = m_calibrationMap.find(detID);
        if (detIt != m_calibrationMap.end()) {
            const auto planeIt = detIt->second.find(planeType);
            if (planeIt != detIt->second.end()) {
                const auto stripIt = planeIt->second.find(stripID);
                if (stripIt != planeIt->second.end()) {
                    const auto& calibration = stripIt->second;
                    gain = calibration.gain;

                    for (auto& sample : calibrated) {
                        const double value =
                            (sample - calibration.pedestal)
                            * calibration.gain
                            * calibration.polarity;

                        sample = static_cast<short>(std::clamp(
                            std::lround(value),
                            static_cast<long>(
                                std::numeric_limits<short>::min()
                            ),
                            static_cast<long>(
                                std::numeric_limits<short>::max()
                            )
                        ));
                    }
                }
            }
        }

        const double sigma = GetSigma(detID, planeType, stripID);
        if (sigma >= 0) {
            const auto maxIt = std::max_element(
                calibrated.begin(),
                calibrated.end()
            );

            if (maxIt == calibrated.end()) {
                continue;
            }

            const double threshold = 5.0 * sigma * std::abs(gain);
            if (static_cast<double>(*maxIt) <= threshold) {
                continue;
            }
        }

        result[detID].push_back(
            RawData{stripID, planeType, std::move(calibrated)}
        );
    }

    return result;
}
