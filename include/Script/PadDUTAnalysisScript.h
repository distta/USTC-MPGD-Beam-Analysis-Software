#pragma once

#include "Script/Base/IScript.h"

#include <string>
#include <vector>

struct PadDUTAlignmentConfig {
    std::vector<std::string> parameters{"dx", "dy", "rotZ"};
    int minMatches{100};
    double huberDelta{-1.0};
    double maxTranslation{20.0};
    double maxRotation{0.05};
    double tolerance{0.005};
    int maxFunctionCalls{100000};
    int maxPasses{3};
    double convergenceTranslation{0.01};
    double convergenceRotation{0.0001};
};

class PadDUTAnalysisScript : public IScript {
   public:
    std::string GetName() const override { return "PadDUTAnalysisScript"; }
    std::string GetDescription() const override {
        return "Two-dimensional pad DUT residual and efficiency analysis";
    }

    void LoadConfig(const json& config) override;
    void Print() const override;
    bool Execute() override;

   private:
    bool m_runAlignment{false};
    PadDUTAlignmentConfig m_alignment;
    int m_progressInterval{1000};
    int m_maxEvents{-1};

    int m_effXBins{10};
    int m_effYBins{10};
    int m_effMinEntriesPerBin{1};
    double m_effXMin{0.0};
    double m_effXMax{90.0};
    double m_effYMin{0.0};
    double m_effYMax{90.0};
    std::vector<int> m_effExcludedXBins;
    std::vector<int> m_effExcludedYBins;

    double m_margin{0.6};
    double m_marginScanMin{0.0};
    double m_marginScanMax{5.0};
    double m_marginScanStep{0.1};

    bool m_enableFakeEfficiency{true};
    unsigned int m_fakeSeed{12345};
    int m_fakePartnersPerEvent{20};
};
