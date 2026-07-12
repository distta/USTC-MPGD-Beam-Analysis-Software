#include "Algorithm/Track/TrackerAlignment.h"
#include "Detector/Detector.h"
#include "Event/DetectorFrame.h"

#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <unordered_map>

namespace Tracking {
namespace {

enum class Parameter { X, Y, Rz };
struct Variable { size_t detector; Parameter parameter; };

double Huber(double value, double k) {
    const double a = std::abs(value);
    return a <= k ? 0.5 * value * value : k * (a - 0.5 * k);
}

// Detector's general coordinate helpers construct dynamic ROOT matrices on
// every call. Alignment performs millions of transforms, so keep the same ZYX
// rotation as plain scalars for the duration of one objective evaluation.
struct Transform {
    TVector3 position;
    double r00, r01, r02, r10, r11, r12, r20, r21, r22;

    explicit Transform(const Detector& detector) : position(detector.GetPos()) {
        const auto rotation = detector.GetRot();
        const double cx = std::cos(rotation.X()), sx = std::sin(rotation.X());
        const double cy = std::cos(rotation.Y()), sy = std::sin(rotation.Y());
        const double cz = std::cos(rotation.Z()), sz = std::sin(rotation.Z());
        r00 = cz * cy;
        r01 = cz * sy * sx - sz * cx;
        r02 = cz * sy * cx + sz * sx;
        r10 = sz * cy;
        r11 = sz * sy * sx + cz * cx;
        r12 = sz * sy * cx - cz * sx;
        r20 = -sy;
        r21 = cy * sx;
        r22 = cy * cx;
    }

    TVector3 ToGlobal(const TVector3& local) const {
        return position + TVector3(r00 * local.X() + r01 * local.Y() + r02 * local.Z(),
                                   r10 * local.X() + r11 * local.Y() + r12 * local.Z(),
                                   r20 * local.X() + r21 * local.Y() + r22 * local.Z());
    }

    TVector3 ToLocal(const TVector3& global) const {
        const auto shifted = global - position;
        return TVector3(r00 * shifted.X() + r10 * shifted.Y() + r20 * shifted.Z(),
                        r01 * shifted.X() + r11 * shifted.Y() + r21 * shifted.Z(),
                        r02 * shifted.X() + r12 * shifted.Y() + r22 * shifted.Z());
    }
};

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

    // Track finding is substantially more expensive than the alignment fit and the
    // hit association is stable for the deliberately small per-iteration steps.
    // Build it once, spread the sample over the full run, and cache only the local
    // coordinates needed by the objective.
    struct CachedHit { size_t detector; TVector3 local; };
    struct CachedTrack { std::vector<CachedHit> hits; };
    std::unordered_map<int, size_t> detectorIndex;
    for (size_t i = 0; i < m_detectors.size(); ++i) detectorIndex[m_detectors[i]->GetID()] = i;

    std::vector<CachedTrack> tracks;
    if (m_config.maxTracks > 0) tracks.reserve(m_config.maxTracks);
    std::mt19937_64 sampling(0x5eedULL);
    size_t usableTracks = 0;
    const size_t eventBudget = m_config.maxEvents > 0
        ? std::min(events.size(), static_cast<size_t>(m_config.maxEvents)) : events.size();
    for (size_t sample = 0; sample < eventBudget; ++sample) {
        const size_t eventIndex = sample * events.size() / eventBudget;
        const auto& event = events[eventIndex];
        for (const auto& result : m_reconstructor.Reconstruct(event)) {
            CachedTrack track;
            track.hits.reserve(result.hitIndices.size());
            for (const auto& [id, hit] : result.hitIndices) {
                const auto detector = detectorIndex.find(id);
                const auto frame = event.detectorFramesMap.find(id);
                if (detector == detectorIndex.end() || frame == event.detectorFramesMap.end() ||
                    hit < 0 || static_cast<size_t>(hit) >= frame->second->LocalHits().size()) continue;
                track.hits.push_back({detector->second, frame->second->LocalHits()[hit].localPos});
            }
            if (track.hits.size() < 3) continue;
            ++usableTracks;
            if (m_config.maxTracks <= 0 || static_cast<int>(tracks.size()) < m_config.maxTracks) {
                tracks.push_back(std::move(track));
            } else {
                // Reservoir sampling keeps the bounded cache representative of the
                // complete run instead of favoring its first events.
                std::uniform_int_distribution<size_t> replacement(0, usableTracks - 1);
                const size_t slot = replacement(sampling);
                if (slot < tracks.size()) tracks[slot] = std::move(track);
            }
        }
    }
    if (static_cast<int>(tracks.size()) < m_config.minTracks) {
        std::cerr << "[TrackAlign] stop: only " << tracks.size() << " tracks\n";
        return false;
    }
    if (m_config.debug) {
        std::cout << "[TrackAlign] cached " << tracks.size() << " tracks from up to "
                  << eventBudget << " uniformly sampled events\n";
    }

    int stableIterations = 0;
    for (int iteration = 0; iteration < m_config.maxIterations; ++iteration) {
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
            std::vector<TVector3> globalHits(m_detectors.size());
            std::vector<Transform> geometry;
            geometry.reserve(m_detectors.size());
            for (const auto& detector : m_detectors) geometry.emplace_back(*detector);
            const auto& reconstruction = m_reconstructor.GetConfig();
            for (const auto& cached : tracks) {
                double sumZ = 0.0, sumX = 0.0, sumY = 0.0;
                double sumZZ = 0.0, sumZX = 0.0, sumZY = 0.0;
                for (const auto& hit : cached.hits) {
                    const auto global = geometry[hit.detector].ToGlobal(hit.local);
                    globalHits[hit.detector] = global;
                    sumZ += global.Z();
                    sumX += global.X();
                    sumY += global.Y();
                    sumZZ += global.Z() * global.Z();
                    sumZX += global.Z() * global.X();
                    sumZY += global.Z() * global.Y();
                }
                for (const auto& target : cached.hits) {
                    const auto& omitted = globalHits[target.detector];
                    const double n = static_cast<double>(cached.hits.size() - 1);
                    const double z = sumZ - omitted.Z();
                    const double x = sumX - omitted.X();
                    const double y = sumY - omitted.Y();
                    const double szz = sumZZ - omitted.Z() * omitted.Z() - z * z / n;
                    if (szz < 1e-12) continue;
                    Track track{};
                    track.kx = (sumZX - omitted.Z() * omitted.X() - z * x / n) / szz;
                    track.ky = (sumZY - omitted.Z() * omitted.Y() - z * y / n) / szz;
                    track.bx = x / n - track.kx * z / n;
                    track.by = y / n - track.ky * z / n;
                    auto& detector = m_detectors[target.detector];
                    const auto predicted = geometry[target.detector].ToLocal(detector->CalcHitFromTrack(track));
                    loss += Huber((target.local.X() - predicted.X()) / reconstruction.resolutionX, m_config.huberK);
                    loss += Huber((target.local.Y() - predicted.Y()) / reconstruction.resolutionY, m_config.huberK);
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
        minimizer->SetMaxFunctionCalls(std::max(1, m_config.maxFunctionCalls));
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
