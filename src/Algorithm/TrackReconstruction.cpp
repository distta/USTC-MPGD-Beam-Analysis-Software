#include "Algorithm/TrackReconstruction.h"
#include "Detector/Detector.h"
#include "Event/DetectorFrame.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <set>
#include <sstream>

namespace Tracking {
namespace {

struct Candidate {
    std::map<int, int> indices;
    std::vector<TVector3> hits;
    Track track{};
    double score = -1e30;
};

bool Finite(const Track& t) {
    return std::isfinite(t.kx) && std::isfinite(t.ky) && std::isfinite(t.bx) &&
           std::isfinite(t.by) && std::isfinite(t.chi2);
}

std::string Key(const Candidate& c) {
    std::ostringstream out;
    for (const auto& [id, hit] : c.indices) out << id << ':' << hit << ';';
    return out.str();
}

bool Conflict(const Candidate& a, const Candidate& b) {
    for (const auto& [id, hit] : a.indices) {
        const auto it = b.indices.find(id);
        if (it != b.indices.end() && it->second == hit) return true;
    }
    return false;
}

}  // namespace

Track FitWeighted(const std::vector<TVector3>& hits, double sigmaX, double sigmaY) {
    Track result{};
    if (hits.size() < 2 || sigmaX <= 0 || sigmaY <= 0) {
        result.chi2 = std::numeric_limits<double>::infinity();
        return result;
    }

    double zMean = 0.0;
    for (const auto& hit : hits) zMean += hit.Z();
    zMean /= hits.size();

    double szz = 0.0, szx = 0.0, szy = 0.0, xMean = 0.0, yMean = 0.0;
    for (const auto& hit : hits) {
        xMean += hit.X();
        yMean += hit.Y();
    }
    xMean /= hits.size();
    yMean /= hits.size();
    for (const auto& hit : hits) {
        const double dz = hit.Z() - zMean;
        szz += dz * dz;
        szx += dz * (hit.X() - xMean);
        szy += dz * (hit.Y() - yMean);
    }
    if (szz < 1e-12) {
        result.chi2 = std::numeric_limits<double>::infinity();
        return result;
    }

    result.kx = szx / szz;
    result.ky = szy / szz;
    result.bx = xMean - result.kx * zMean;
    result.by = yMean - result.ky * zMean;
    double chi2 = 0.0;
    for (const auto& hit : hits) {
        const double rx = (hit.X() - result.kx * hit.Z() - result.bx) / sigmaX;
        const double ry = (hit.Y() - result.ky * hit.Z() - result.by) / sigmaY;
        chi2 += rx * rx + ry * ry;
    }
    const int ndf = std::max(1, static_cast<int>(2 * hits.size()) - 4);
    result.chi2 = chi2 / ndf;
    return result;
}

Reconstructor::Reconstructor(std::vector<std::shared_ptr<Detector>> detectors, Config config)
    : m_detectors(std::move(detectors)), m_config(config) {
    std::sort(m_detectors.begin(), m_detectors.end(), [](const auto& a, const auto& b) {
        return a->GetPos().Z() < b->GetPos().Z();
    });
}

std::vector<Result> Reconstructor::Reconstruct(const Event& event, ReconstructionStats* stats) const {
    if (stats) *stats = {};
    if (m_detectors.size() < 3) return {};
    const int minHits = m_detectors.size() <= 3 ? m_detectors.size() : m_detectors.size() - 1;
    std::vector<std::shared_ptr<Detector>> active;
    for (const auto& detector : m_detectors) {
        const auto frame = event.detectorFramesMap.find(detector->GetID());
        if (frame != event.detectorFramesMap.end() && !frame->second->LocalHits().empty()) active.push_back(detector);
    }
    if (static_cast<int>(active.size()) < minHits) return {};
    const auto& firstHits = event.detectorFramesMap.at(active.front()->GetID())->LocalHits();
    const auto& lastHits = event.detectorFramesMap.at(active.back()->GetID())->LocalHits();
    if (firstHits.empty() || lastHits.empty()) return {};

    std::vector<Candidate> candidates;
    candidates.reserve(std::min<int>(m_config.maxCandidates, firstHits.size() * lastHits.size()));
    for (size_t i = 0; i < firstHits.size(); ++i) {
        for (size_t j = 0; j < lastHits.size(); ++j) {
            Candidate c;
            c.indices[active.front()->GetID()] = i;
            c.indices[active.back()->GetID()] = j;
            c.hits = {active.front()->LocalToGlobal(firstHits[i].localPos),
                      active.back()->LocalToGlobal(lastHits[j].localPos)};
            c.track = FitWeighted(c.hits, m_config.resolutionX, m_config.resolutionY);
            if (Finite(c.track)) candidates.push_back(std::move(c));
            if (static_cast<int>(candidates.size()) >= m_config.maxCandidates) break;
        }
        if (static_cast<int>(candidates.size()) >= m_config.maxCandidates) break;
    }
    if (stats) stats->seedCandidates = candidates.size();

    for (size_t layer = 1; layer + 1 < active.size() && !candidates.empty(); ++layer) {
        const auto& detector = active[layer];
        const auto& localHits = event.detectorFramesMap.at(detector->GetID())->LocalHits();
        std::vector<Candidate> next;
        for (const auto& candidate : candidates) {
            const auto predictedLocal = detector->GlobalToLocal(detector->CalcHitFromTrack(candidate.track));
            std::vector<std::pair<double, int>> matches;
            for (size_t hit = 0; hit < localHits.size(); ++hit) {
                const double dx = (localHits[hit].localPos.X() - predictedLocal.X()) / m_config.resolutionX;
                const double dy = (localHits[hit].localPos.Y() - predictedLocal.Y()) / m_config.resolutionY;
                const double distance = std::hypot(dx, dy);
                if (distance <= m_config.gateSigma) matches.emplace_back(distance, hit);
            }
            std::sort(matches.begin(), matches.end());
            const int branches = std::min<int>(m_config.maxBranchesPerLayer, matches.size());
            for (int k = 0; k < branches; ++k) {
                Candidate child = candidate;
                child.indices[detector->GetID()] = matches[k].second;
                child.hits.push_back(detector->LocalToGlobal(localHits[matches[k].second].localPos));
                child.track = FitWeighted(child.hits, m_config.resolutionX, m_config.resolutionY);
                if (Finite(child.track)) next.push_back(std::move(child));
            }
            if (minHits < static_cast<int>(m_detectors.size())) next.push_back(candidate);
        }
        std::sort(next.begin(), next.end(), [](const Candidate& a, const Candidate& b) {
            if (a.indices.size() != b.indices.size()) return a.indices.size() > b.indices.size();
            return a.track.chi2 < b.track.chi2;
        });
        if (static_cast<int>(next.size()) > m_config.maxCandidates) next.resize(m_config.maxCandidates);
        candidates = std::move(next);
    }

    std::set<std::string> seen;
    std::vector<Candidate> valid;
    for (auto& candidate : candidates) {
        if (static_cast<int>(candidate.indices.size()) < minHits ||
            candidate.track.chi2 > m_config.maxChi2Ndf ||
            !seen.insert(Key(candidate)).second) continue;
        candidate.score = 100.0 * candidate.indices.size() - std::min(99.0, candidate.track.chi2);
        valid.push_back(std::move(candidate));
    }
    std::sort(valid.begin(), valid.end(), [](const Candidate& a, const Candidate& b) { return a.score > b.score; });
    if (stats) stats->finalCandidates = valid.size();
    if (valid.size() > 80) valid.resize(80);  // bounded exact conflict search

    std::vector<int> best, current;
    double bestScore = 0.0;
    int nodes = 0;
    std::function<void(size_t, double)> search = [&](size_t pos, double score) {
        if (++nodes > m_config.conflictSearchNodes || current.size() >= static_cast<size_t>(m_config.maxTracks) || pos == valid.size()) {
            if (score > bestScore) { bestScore = score; best = current; }
            return;
        }
        double upper = score;
        for (size_t i = pos; i < valid.size(); ++i) upper += std::max(0.0, valid[i].score);
        if (upper <= bestScore) return;
        bool conflict = false;
        for (int selected : current) conflict |= Conflict(valid[pos], valid[selected]);
        if (!conflict) {
            current.push_back(pos);
            search(pos + 1, score + valid[pos].score);
            current.pop_back();
        }
        search(pos + 1, score);
    };
    search(0, 0.0);
    if (stats) stats->conflictSearchNodes = nodes;

    std::vector<Result> results;
    for (int index : best) results.push_back({valid[index].track, valid[index].indices, valid[index].track.chi2});
    return results;
}

}  // namespace Tracking
