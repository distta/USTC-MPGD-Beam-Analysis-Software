#include "Algorithm/TrackerAlignment.h"
#include "Detector/Detector.h"
#include "Event/DetectorFrame.h"

#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <tuple>

namespace Tracking {
namespace {

enum class Parameter { X, Y, Rz };
struct Variable { size_t detector; Parameter parameter; };

double Huber(double value, double k) {
    const double a = std::abs(value);
    return a <= k ? 0.5 * value * value : k * (a - 0.5 * k);
}

}  // namespace

Aligner::Aligner(std::vector<std::shared_ptr<Detector>> detectors,
                 const Reconstructor& reconstructor, AlignmentConfig config)
    : m_detectors(std::move(detectors)), m_reconstructor(reconstructor), m_config(config) {
    std::sort(m_detectors.begin(), m_detectors.end(), [](const auto& a, const auto& b) {
        return a->GetPos().Z() < b->GetPos().Z();
    });
}

bool Aligner::Run(const std::vector<Event>& events) {
    m_history.clear();
    if (m_detectors.size() < 3) return false;
    std::vector<Variable> variables;
    for (size_t i = 1; i < m_detectors.size(); ++i) {
        if (i + 1 < m_detectors.size()) {
            variables.push_back({i, Parameter::X});
            variables.push_back({i, Parameter::Y});
        }
        variables.push_back({i, Parameter::Rz});
    }
    if (variables.empty()) return true;

    int stableIterations = 0;
    for (int iteration = 0; iteration < m_config.maxIterations; ++iteration) {
        struct AlignedTrack { const Event* event; Result result; };
        std::vector<AlignedTrack> tracks;
        for (const auto& event : events) {
            for (auto& result : m_reconstructor.Reconstruct(event)) tracks.push_back({&event, std::move(result)});
        }
        if (static_cast<int>(tracks.size()) < m_config.minTracks) {
            std::cerr << "[TrackAlign] stop: only " << tracks.size() << " tracks\n";
            return iteration > 0;
        }

        std::vector<TVector3> basePos, baseRot;
        for (const auto& detector : m_detectors) {
            basePos.push_back(detector->GetAlignPos());
            baseRot.push_back(detector->GetAlignRot());
        }
        auto apply = [&](const double* parameters) {
            auto pos = basePos;
            auto rot = baseRot;
            for (size_t i = 0; i < variables.size(); ++i) {
                const auto& variable = variables[i];
                if (variable.parameter == Parameter::X) pos[variable.detector].SetX(pos[variable.detector].X() + parameters[i]);
                if (variable.parameter == Parameter::Y) pos[variable.detector].SetY(pos[variable.detector].Y() + parameters[i]);
                if (variable.parameter == Parameter::Rz) rot[variable.detector].SetZ(rot[variable.detector].Z() + parameters[i]);
            }
            for (size_t i = 0; i < m_detectors.size(); ++i)
                m_detectors[i]->SetAlignment(pos[i].X(), pos[i].Y(), pos[i].Z(), rot[i].X(), rot[i].Y(), rot[i].Z());
        };

        auto objective = [&](const double* parameters) {
            apply(parameters);
            double loss = 0.0;
            size_t residuals = 0;
            for (const auto& aligned : tracks) {
                for (const auto& [targetID, targetHit] : aligned.result.hitIndices) {
                    std::vector<TVector3> fitHits;
                    for (const auto& [id, hit] : aligned.result.hitIndices) {
                        if (id == targetID) continue;
                        const auto detector = std::find_if(m_detectors.begin(), m_detectors.end(),
                            [id](const auto& d) { return d->GetID() == id; });
                        fitHits.push_back((*detector)->LocalToGlobal(
                            aligned.event->detectorFramesMap.at(id)->LocalHits().at(hit).localPos));
                    }
                    if (fitHits.size() < 2) continue;
                    const Track track = FitWeighted(fitHits, m_reconstructor.GetConfig().resolutionX,
                                                   m_reconstructor.GetConfig().resolutionY);
                    if (!std::isfinite(track.chi2)) continue;
                    const auto target = std::find_if(m_detectors.begin(), m_detectors.end(),
                        [targetID](const auto& d) { return d->GetID() == targetID; });
                    const auto& measured = aligned.event->detectorFramesMap.at(targetID)->LocalHits().at(targetHit).localPos;
                    const auto predicted = (*target)->GlobalToLocal((*target)->CalcHitFromTrack(track));
                    loss += Huber((measured.X() - predicted.X()) / m_reconstructor.GetConfig().resolutionX, m_config.huberK);
                    loss += Huber((measured.Y() - predicted.Y()) / m_reconstructor.GetConfig().resolutionY, m_config.huberK);
                    residuals += 2;
                }
            }
            return residuals ? loss / residuals : 1e30;
        };

        auto minimizer = std::unique_ptr<ROOT::Math::Minimizer>(
            ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
        ROOT::Math::Functor function(objective, variables.size());
        minimizer->SetFunction(function);
        minimizer->SetPrintLevel(0);
        minimizer->SetMaxFunctionCalls(3000);
        minimizer->SetTolerance(1e-4);
        for (size_t i = 0; i < variables.size(); ++i) {
            const bool rotation = variables[i].parameter == Parameter::Rz;
            const double limit = rotation ? m_config.maxRotationStep : m_config.maxShiftStep;
            const std::string axis = variables[i].parameter == Parameter::X ? "dx" :
                                     variables[i].parameter == Parameter::Y ? "dy" : "rz";
            minimizer->SetLimitedVariable(i, axis + std::to_string(m_detectors[variables[i].detector]->GetID()),
                                          0.0, limit / 20.0, -limit, limit);
        }
        const double before = objective(std::vector<double>(variables.size(), 0.0).data());
        const bool ok = minimizer->Minimize();
        const double* solution = minimizer->X();
        const double after = minimizer->MinValue();
        if (!ok || !std::isfinite(after) || after >= before) {
            std::vector<double> zero(variables.size(), 0.0);
            apply(zero.data());
            std::cerr << "[TrackAlign] iteration " << iteration + 1 << " rejected: loss "
                      << before << " -> " << after << '\n';
            return iteration > 0;
        }
        apply(solution);
        double maxShift = 0.0, maxRotation = 0.0;
        for (size_t i = 0; i < variables.size(); ++i) {
            if (variables[i].parameter == Parameter::Rz) maxRotation = std::max(maxRotation, std::abs(solution[i]));
            else maxShift = std::max(maxShift, std::abs(solution[i]));
        }
        AlignmentIteration history;
        history.iteration = iteration + 1;
        history.tracks = tracks.size();
        history.lossBefore = before;
        history.lossAfter = after;
        history.maxShift = maxShift;
        history.maxRotation = maxRotation;
        for (const auto& detector : m_detectors) {
            history.alignPosition[detector->GetID()] = detector->GetAlignPos();
            history.alignRotation[detector->GetID()] = detector->GetAlignRot();
        }
        m_history.push_back(std::move(history));
        if (m_config.debug) {
            std::cout << "[TrackAlign] iter=" << iteration + 1 << " tracks=" << tracks.size()
                      << " loss=" << before << "->" << after << " maxShift=" << maxShift
                      << " mm maxRz=" << maxRotation << " rad\n";
        }
        if (maxShift < m_config.shiftTolerance && maxRotation < m_config.rotationTolerance) ++stableIterations;
        else stableIterations = 0;
        if (stableIterations >= 2) break;
    }
    return true;
}

}  // namespace Tracking
