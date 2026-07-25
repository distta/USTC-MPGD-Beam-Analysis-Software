#include "Algorithm/Track/TrackReconstruction.h"
#include "Detector/Detector.h"
#include "Detector/Planar.h"
#include "Event/DetectorFrame.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <unordered_set>

namespace Tracking {
namespace {

struct FitSums {
    int n = 0;
    double z = 0, x = 0, y = 0;
    double zz = 0, zx = 0, zy = 0, xx = 0, yy = 0;

    void Add(const TVector3& hit) {
        ++n;
        z += hit.Z(); x += hit.X(); y += hit.Y();
        zz += hit.Z() * hit.Z();
        zx += hit.Z() * hit.X(); zy += hit.Z() * hit.Y();
        xx += hit.X() * hit.X(); yy += hit.Y() * hit.Y();
    }
};

Track Fit(const FitSums& sums, double sigmaX, double sigmaY) {
    Track result{};
    if (sums.n < 2 || sigmaX <= 0 || sigmaY <= 0) {
        result.chi2 = std::numeric_limits<double>::infinity();
        return result;
    }
    const double n = sums.n;
    const double centeredZZ = sums.zz - sums.z * sums.z / n;
    if (centeredZZ < 1e-12) {
        result.chi2 = std::numeric_limits<double>::infinity();
        return result;
    }
    result.kx = (sums.zx - sums.z * sums.x / n) / centeredZZ;
    result.ky = (sums.zy - sums.z * sums.y / n) / centeredZZ;
    result.bx = (sums.x - result.kx * sums.z) / n;
    result.by = (sums.y - result.ky * sums.z) / n;
    // At the least-squares solution SSE = y'y - beta'X'y.  Keeping the
    // sufficient statistics makes extending a candidate O(1), not O(layers).
    const double sseX = std::max(0.0, sums.xx - result.bx * sums.x - result.kx * sums.zx);
    const double sseY = std::max(0.0, sums.yy - result.by * sums.y - result.ky * sums.zy);
    const int ndf = std::max(1, 2 * sums.n - 4);
    result.chi2 = (sseX / (sigmaX * sigmaX) + sseY / (sigmaY * sigmaY)) / ndf;
    return result;
}

struct Candidate {
    std::vector<int> indices;  // indexed by active detector layer; -1 means missing
    FitSums sums;
    Track track{};
    int hitCount = 0;
    double score = -1e30;
};

bool Finite(const Track& track) {
    return std::isfinite(track.kx) && std::isfinite(track.ky) &&
           std::isfinite(track.bx) && std::isfinite(track.by) && std::isfinite(track.chi2);
}

struct IndicesHash {
    size_t operator()(const std::vector<int>& indices) const noexcept {
        size_t hash = 1469598103934665603ULL;
        for (int value : indices) {
            hash ^= static_cast<size_t>(value + 2);
            hash *= 1099511628211ULL;
        }
        return hash;
    }
};

// Detector's generic helpers allocate ROOT matrices on every call.  Cache the
// same ZYX matrix once per layer/event; alignment changes are still picked up
// by every Reconstruct invocation.
struct Geometry {
    const Detector* detector = nullptr;
    TVector3 position;
    double r00, r01, r02, r10, r11, r12, r20, r21, r22;
    bool planar = false;

    explicit Geometry(const Detector& value)
        : detector(&value), position(value.GetPos()),
          planar(value.GetPlanarConfig() != nullptr || value.GetPlanarPadConfig() != nullptr) {
        const auto rotation = value.GetRot();
        const double cx = std::cos(rotation.X()), sx = std::sin(rotation.X());
        const double cy = std::cos(rotation.Y()), sy = std::sin(rotation.Y());
        const double cz = std::cos(rotation.Z()), sz = std::sin(rotation.Z());
        r00 = cz * cy; r01 = cz * sy * sx - sz * cx; r02 = cz * sy * cx + sz * sx;
        r10 = sz * cy; r11 = sz * sy * sx + cz * cx; r12 = sz * sy * cx - cz * sx;
        r20 = -sy;     r21 = cy * sx;                    r22 = cy * cx;
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

    TVector3 PredictLocal(const Track& track) const {
        if (!planar) return ToLocal(detector->CalcHitFromTrack(track));
        // r02/r12/r22 is the plane normal.  The unnormalised track direction
        // avoids two Unit() calls performed by Planar::CalcHitFromTrack.
        const double denominator = track.kx * r02 + track.ky * r12 + r22;
        if (std::abs(denominator) < 1e-12) return {};
        const double numerator = (position.X() - track.bx) * r02 +
                                 (position.Y() - track.by) * r12 + position.Z() * r22;
        const double z = numerator / denominator;
        if (z < 0) return {}; // preserve Planar's forward-intersection check
        return ToLocal({track.kx * z + track.bx, track.ky * z + track.by, z});
    }
};

}  // namespace

Track FitWeighted(const std::vector<TVector3>& hits, double sigmaX, double sigmaY) {
    FitSums sums;
    for (const auto& hit : hits) sums.Add(hit);
    return Fit(sums, sigmaX, sigmaY);
}

Reconstructor::Reconstructor(std::vector<std::shared_ptr<Detector>> detectors, Config config)
    : m_detectors(std::move(detectors)), m_config(config) {
    std::sort(m_detectors.begin(), m_detectors.end(), [](const auto& a, const auto& b) {
        return a->GetPos().Z() < b->GetPos().Z();
    });
}

std::vector<Result> Reconstructor::Reconstruct(const Event& event, ReconstructionStats* stats) const {
    if (stats) *stats = {};
    if (m_detectors.size() < 3 || m_config.maxCandidates <= 0 || m_config.maxTracks <= 0) return {};
    const int minHits = m_detectors.size() <= 3 ? m_detectors.size() : m_detectors.size() - 1;

    struct Layer {
        std::shared_ptr<Detector> detector;
        const std::vector<LocalHit>* localHits;
        Geometry geometry;
        std::vector<TVector3> globalHits;

        Layer(std::shared_ptr<Detector> detector_, const std::vector<LocalHit>& hits)
            : detector(std::move(detector_)), localHits(&hits), geometry(*detector) {
            globalHits.reserve(hits.size());
            for (const auto& hit : hits) globalHits.push_back(geometry.ToGlobal(hit.localPos));
        }
    };

    std::vector<Layer> active;
    active.reserve(m_detectors.size());
    for (const auto& detector : m_detectors) {
        const auto frame = event.detectorFramesMap.find(detector->GetID());
        if (frame != event.detectorFramesMap.end() && !frame->second->LocalHits().empty())
            active.emplace_back(detector, frame->second->LocalHits());
    }
    if (static_cast<int>(active.size()) < minHits) return {};

    const auto& firstHits = *active.front().localHits;
    const auto& lastHits = *active.back().localHits;
    std::vector<Candidate> candidates;
    const size_t seedCount = std::min<size_t>(m_config.maxCandidates, firstHits.size() * lastHits.size());
    candidates.reserve(seedCount);
    for (size_t i = 0; i < firstHits.size() && candidates.size() < seedCount; ++i) {
        for (size_t j = 0; j < lastHits.size() && candidates.size() < seedCount; ++j) {
            Candidate candidate;
            candidate.indices.assign(active.size(), -1);
            candidate.indices.front() = static_cast<int>(i);
            candidate.indices.back() = static_cast<int>(j);
            candidate.sums.Add(active.front().globalHits[i]);
            candidate.sums.Add(active.back().globalHits[j]);
            candidate.track = Fit(candidate.sums, m_config.resolutionX, m_config.resolutionY);
            candidate.hitCount = 2;
            if (Finite(candidate.track)) candidates.push_back(std::move(candidate));
        }
    }
    if (stats) stats->seedCandidates = candidates.size();

    const auto better = [](const Candidate& a, const Candidate& b) {
        if (a.hitCount != b.hitCount) return a.hitCount > b.hitCount;
        return a.track.chi2 < b.track.chi2;
    };
    for (size_t layerIndex = 1; layerIndex + 1 < active.size() && !candidates.empty(); ++layerIndex) {
        const auto& layer = active[layerIndex];
        std::vector<Candidate> next;
        next.reserve(std::min<size_t>(m_config.maxCandidates,
            candidates.size() * static_cast<size_t>(std::max(1, m_config.maxBranchesPerLayer))));
        for (const auto& candidate : candidates) {
            const auto predicted = layer.geometry.PredictLocal(candidate.track);
            std::vector<std::pair<double, int>> matches;
            const int branchLimit = std::max(0, m_config.maxBranchesPerLayer);
            matches.reserve(std::min<size_t>(branchLimit, layer.localHits->size()));
            const double gate2 = m_config.gateSigma * m_config.gateSigma;
            for (size_t hit = 0; hit < layer.localHits->size(); ++hit) {
                const double dx = ((*layer.localHits)[hit].localPos.X() - predicted.X()) / m_config.resolutionX;
                const double dy = ((*layer.localHits)[hit].localPos.Y() - predicted.Y()) / m_config.resolutionY;
                const double distance2 = dx * dx + dy * dy;
                if (distance2 > gate2 || branchLimit == 0) continue;
                const auto position = std::lower_bound(matches.begin(), matches.end(), distance2,
                    [](const auto& match, double value) { return match.first < value; });
                matches.insert(position, {distance2, static_cast<int>(hit)});
                if (static_cast<int>(matches.size()) > branchLimit) matches.pop_back();
            }
            for (const auto& match : matches) {
                Candidate child = candidate;
                child.indices[layerIndex] = match.second;
                child.sums.Add(layer.globalHits[match.second]);
                child.track = Fit(child.sums, m_config.resolutionX, m_config.resolutionY);
                ++child.hitCount;
                if (Finite(child.track)) next.push_back(std::move(child));
            }
            if (minHits < static_cast<int>(m_detectors.size())) next.push_back(candidate);
        }
        if (next.size() > static_cast<size_t>(m_config.maxCandidates)) {
            std::nth_element(next.begin(), next.begin() + m_config.maxCandidates, next.end(), better);
            next.resize(m_config.maxCandidates);
        }
        candidates = std::move(next);
    }

    std::unordered_set<std::vector<int>, IndicesHash> seen;
    seen.reserve(candidates.size());
    std::vector<Candidate> valid;
    valid.reserve(candidates.size());
    for (auto& candidate : candidates) {
        if (candidate.hitCount < minHits || candidate.track.chi2 > m_config.maxChi2Ndf ||
            !seen.insert(candidate.indices).second) continue;
        candidate.score = 100.0 * candidate.hitCount - std::min(99.0, candidate.track.chi2);
        valid.push_back(std::move(candidate));
    }
    std::sort(valid.begin(), valid.end(), [](const Candidate& a, const Candidate& b) {
        return a.score > b.score;
    });
    if (stats) stats->finalCandidates = valid.size();
    if (valid.size() > 80) valid.resize(80);
    if (valid.empty()) return {};

    // Precompute conflicts into two machine words.  This replaces map lookups
    // throughout the bounded branch-and-bound search.
    std::vector<std::array<uint64_t, 2>> conflicts(valid.size(), {0, 0});
    for (size_t i = 0; i < valid.size(); ++i) {
        for (size_t j = i + 1; j < valid.size(); ++j) {
            bool conflict = false;
            for (size_t layer = 0; layer < active.size(); ++layer) {
                if (valid[i].indices[layer] >= 0 && valid[i].indices[layer] == valid[j].indices[layer]) {
                    conflict = true;
                    break;
                }
            }
            if (conflict) {
                conflicts[i][j / 64] |= uint64_t{1} << (j % 64);
                conflicts[j][i / 64] |= uint64_t{1} << (i % 64);
            }
        }
    }

    std::vector<int> best, current;
    std::array<uint64_t, 2> greedyBlocked{0, 0};
    double bestScore = 0;
    for (size_t i = 0; i < valid.size() && best.size() < static_cast<size_t>(m_config.maxTracks); ++i) {
        if (greedyBlocked[i / 64] & (uint64_t{1} << (i % 64))) continue;
        best.push_back(static_cast<int>(i));
        bestScore += valid[i].score;
        greedyBlocked[0] |= conflicts[i][0];
        greedyBlocked[1] |= conflicts[i][1];
    }

    int nodes = 0;
    const int nodeLimit = std::max(1, m_config.conflictSearchNodes);
    auto search = [&](auto&& self, size_t pos, double score, std::array<uint64_t, 2> blocked) -> void {
        if (++nodes > nodeLimit) return;
        // current is always a valid solution.  Retain it even when the node
        // budget expires before this branch reaches a leaf.
        if (score > bestScore) { bestScore = score; best = current; }
        while (pos < valid.size() && (blocked[pos / 64] & (uint64_t{1} << (pos % 64)))) ++pos;
        if (pos == valid.size() || current.size() >= static_cast<size_t>(m_config.maxTracks)) {
            if (score > bestScore) { bestScore = score; best = current; }
            return;
        }
        double upper = score;
        size_t slots = static_cast<size_t>(m_config.maxTracks) - current.size();
        for (size_t i = pos; i < valid.size() && slots; ++i) {
            if (!(blocked[i / 64] & (uint64_t{1} << (i % 64)))) {
                upper += valid[i].score;
                --slots;
            }
        }
        if (upper <= bestScore) return;

        current.push_back(static_cast<int>(pos));
        auto includeBlocked = blocked;
        includeBlocked[0] |= conflicts[pos][0];
        includeBlocked[1] |= conflicts[pos][1];
        self(self, pos + 1, score + valid[pos].score, includeBlocked);
        current.pop_back();
        blocked[pos / 64] |= uint64_t{1} << (pos % 64);
        self(self, pos + 1, score, blocked);
    };
    search(search, 0, 0.0, {0, 0});
    if (stats) stats->conflictSearchNodes = nodes;

    std::vector<Result> results;
    results.reserve(best.size());
    for (int index : best) {
        std::map<int, int> hitIndices;
        for (size_t layer = 0; layer < active.size(); ++layer) {
            if (valid[index].indices[layer] >= 0)
                hitIndices[active[layer].detector->GetID()] = valid[index].indices[layer];
        }
        results.push_back({valid[index].track, std::move(hitIndices), valid[index].track.chi2});
    }
    return results;
}

}  // namespace Tracking
