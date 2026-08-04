// SPDX-License-Identifier: MIT

#include "delaunay32/delaunay.hpp"
#include "internal.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace delaunay32 {
namespace {

bool opposite_nonzero_signs(std::int64_t a, std::int64_t b) {
    return (a > 0 && b < 0) || (a < 0 && b > 0);
}

}  // namespace

std::vector<Triangle> Triangulator::triangulate_constrained_int(
    const std::vector<Point>& points,
    const std::vector<Constraint>& constraints) {
    if (constraints.empty()) {
        return triangulate_int(points);
    }
    require_point_count(points.size());
    for (const Constraint constraint : constraints) {
        if (constraint.i0 >= points.size() ||
            constraint.i1 >= points.size()) {
            throw std::invalid_argument(
                "constraint endpoint is outside the point array");
        }
        if (constraint.i0 == constraint.i1) {
            throw std::invalid_argument(
                "constraint endpoints are coincident");
        }
    }
    load_int_points(points);
    std::vector<std::uint32_t> representatives;
    build_loaded_topology(&representatives);
    build_constraint_indices(representatives, points.size());
    recover_constraints(constraints);
    return finish_triangle_export();
}

void Triangulator::build_constraint_indices(
    const std::vector<std::uint32_t>& representatives,
    std::size_t original_point_count) {
    original_to_site_.assign(original_point_count, kDeletedEdge);
    site_edge_.assign(points_.size(), kDeletedEdge);
    std::uint32_t marker_count = 0;
    for (const EdgeRange range : edge_ranges_) {
        marker_count = std::max(marker_count, range.last);
    }
    // Keep constraint-only storage and initialization off the ordinary
    // triangulation path.
    edge_constrained_.assign(marker_count, 0);

    for (std::size_t site = 0; site < points_.size(); ++site) {
        original_to_site_[points_[site].original] =
            static_cast<std::uint32_t>(site);
    }
    for (std::size_t original = 0;
         original < original_point_count;
         ++original) {
        if (original_to_site_[original] != kDeletedEdge) {
            continue;
        }
        const std::uint32_t representative = representatives[original];
        original_to_site_[original] = original_to_site_[representative];
    }
}

std::uint32_t Triangulator::find_edge(
    std::uint32_t origin,
    std::uint32_t destination) const {
    const std::uint32_t start = site_edge_[origin];
    if (start == kDeletedEdge || org(start) != origin) {
        throw std::logic_error("invalid site-to-edge mapping");
    }
    std::uint32_t edge = start;
    do {
        if (dest(edge) == destination) {
            return edge;
        }
        edge = onext(edge);
    } while (edge != start);
    return kDeletedEdge;
}

std::uint32_t Triangulator::first_collinear_edge(
    std::uint32_t origin,
    std::uint32_t destination) const {
    const Site& a = points_[origin];
    const Site& b = points_[destination];
    const std::uint32_t start = site_edge_[origin];
    std::uint32_t best = kDeletedEdge;
    std::uint64_t best_distance =
        std::numeric_limits<std::uint64_t>::max();
    std::uint32_t edge = start;
    do {
        const std::uint32_t candidate = dest(edge);
        if (candidate != destination &&
            orient(origin, destination, candidate) == 0) {
            const Site& p = points_[candidate];
            const bool between =
                p.x >= std::min(a.x, b.x) &&
                p.x <= std::max(a.x, b.x) &&
                p.y >= std::min(a.y, b.y) &&
                p.y <= std::max(a.y, b.y) &&
                (p.x != a.x || p.y != a.y);
            if (between) {
                const std::int64_t dx =
                    static_cast<std::int64_t>(p.x) - a.x;
                const std::int64_t dy =
                    static_cast<std::int64_t>(p.y) - a.y;
                const std::uint64_t ux = static_cast<std::uint64_t>(
                    dx < 0 ? -dx : dx);
                const std::uint64_t uy = static_cast<std::uint64_t>(
                    dy < 0 ? -dy : dy);
                const std::uint64_t distance = ux * ux + uy * uy;
                if (distance < best_distance) {
                    best = edge;
                    best_distance = distance;
                }
            }
        }
        edge = onext(edge);
    } while (edge != start);
    return best;
}

bool Triangulator::properly_intersects(
    std::uint32_t edge,
    std::uint32_t a,
    std::uint32_t b) const {
    const std::uint32_t u = org(edge);
    const std::uint32_t v = dest(edge);
    return opposite_nonzero_signs(
               orient(a, b, u), orient(a, b, v)) &&
           opposite_nonzero_signs(
               orient(u, v, a), orient(u, v, b));
}

std::vector<std::uint32_t> Triangulator::crossed_edges(
    std::uint32_t a,
    std::uint32_t b,
    std::uint32_t& reached) const {
    reached = b;
    std::uint32_t crossed = kDeletedEdge;
    const std::uint32_t start = site_edge_[a];
    std::uint32_t outgoing = start;
    do {
        std::uint32_t opposite = 0;
        if (left_face_opposite(outgoing, opposite)) {
            const std::uint32_t candidate = lnext(outgoing);
            if (properly_intersects(candidate, a, b)) {
                crossed = candidate;
                break;
            }
        }
        outgoing = onext(outgoing);
    } while (outgoing != start);
    if (crossed == kDeletedEdge) {
        throw std::logic_error(
            "could not locate the first edge crossed by a constraint");
    }

    std::vector<std::uint32_t> result;
    result.reserve(16);
    std::size_t step_count = 0;
    const std::size_t maximum_steps = points_.size() * 3;
    while (true) {
        result.push_back(crossed & ~1U);
        const std::uint32_t entered = sym(crossed);
        const std::uint32_t first = lnext(entered);
        const std::uint32_t second = lnext(first);
        const std::uint32_t third = dest(first);
        if (third == b) {
            return result;
        }
        if (orient(a, b, third) == 0) {
            // An existing site on the requested segment becomes the endpoint
            // of this recovery pass. The caller continues from it, producing
            // a constraint chain without inserting a Steiner point.
            reached = third;
            return result;
        }
        if (properly_intersects(first, a, b)) {
            crossed = first;
        } else if (properly_intersects(second, a, b)) {
            crossed = second;
        } else {
            throw std::logic_error(
                "constraint walk did not leave the current triangle");
        }
        ++step_count;
        if (step_count > maximum_steps) {
            throw std::logic_error("constraint walk did not terminate");
        }
    }
}

void Triangulator::recover_constraint(
    std::uint32_t a,
    std::uint32_t b,
    std::vector<std::uint32_t>& legalization_queue) {
    while (a != b) {
        const std::uint32_t existing = find_edge(a, b);
        if (existing != kDeletedEdge) {
            mark_constrained(existing);
            return;
        }

        const std::uint32_t collinear = first_collinear_edge(a, b);
        if (collinear != kDeletedEdge) {
            mark_constrained(collinear);
            const std::uint32_t reached = dest(collinear);
            site_edge_[reached] = sym(collinear);
            a = reached;
            continue;
        }

        std::uint32_t reached = b;
        const std::vector<std::uint32_t> crossed =
            crossed_edges(a, b, reached);
        // A crossed diagonal can initially be nonflippable. Cycling those
        // diagonals lets neighboring flips make their quadrilaterals convex.
        std::deque<std::uint32_t> pending(
            crossed.begin(), crossed.end());
        std::size_t blocked = 0;
        std::size_t flips = 0;
        const std::size_t crossing_count = crossed.size();
        const std::size_t maximum =
            std::numeric_limits<std::size_t>::max();
        const bool flip_bound_overflows =
            crossing_count != 0 &&
            crossing_count > (maximum - 64) / 8 / crossing_count;
        const std::size_t maximum_flips = flip_bound_overflows
                                              ? maximum
                                              : 64 + crossing_count *
                                                         crossing_count * 8;

        while (!pending.empty()) {
            const std::uint32_t edge = pending.front();
            pending.pop_front();
            if (!properly_intersects(edge, a, reached)) {
                blocked = 0;
                continue;
            }
            if (is_constrained(edge)) {
                throw std::invalid_argument(
                    "constraints intersect away from an existing site");
            }
            if (!can_flip(edge)) {
                pending.push_back(edge);
                ++blocked;
                if (blocked >= pending.size()) {
                    throw std::logic_error(
                        "constraint recovery has no flippable crossed edge");
                }
                continue;
            }

            flip_edge(edge);
            seed_constraint_legalization(edge, legalization_queue);
            blocked = 0;
            ++flips;
            if (flips > maximum_flips) {
                throw std::logic_error(
                    "constraint recovery exceeded its flip bound");
            }
            if (properly_intersects(edge, a, reached)) {
                pending.push_back(edge);
            }
        }

        const std::uint32_t recovered = find_edge(a, reached);
        if (recovered == kDeletedEdge) {
            throw std::logic_error(
                "constraint recovery did not create the requested edge");
        }
        mark_constrained(recovered);
        site_edge_[reached] = sym(recovered);
        a = reached;
    }
}

void Triangulator::queue_constraint_legalization(
    std::uint32_t edge,
    std::vector<std::uint32_t>& legalization_queue) {
    const std::uint32_t pair = edge & ~1U;
    if ((edge_constrained_[pair] & kLegalizationQueuedBit) != 0) {
        return;
    }
    edge_constrained_[pair] |= kLegalizationQueuedBit;
    edge_constrained_[pair + 1U] |= kLegalizationQueuedBit;
    legalization_queue.push_back(pair);
}

void Triangulator::seed_constraint_legalization(
    std::uint32_t edge,
    std::vector<std::uint32_t>& legalization_queue) {
    // A diagonal swap can change local Delaunay legality only for the new
    // diagonal and the four sides of its quadrilateral. The input topology is
    // already legal, so the union of these neighborhoods is a complete Lawson
    // worklist after constraint recovery.
    queue_constraint_legalization(edge, legalization_queue);
    const std::uint32_t left_next = lnext(edge);
    const std::uint32_t left_previous = lnext(left_next);
    const std::uint32_t right_next = lnext(sym(edge));
    const std::uint32_t right_previous = lnext(right_next);
    queue_constraint_legalization(left_next, legalization_queue);
    queue_constraint_legalization(left_previous, legalization_queue);
    queue_constraint_legalization(right_next, legalization_queue);
    queue_constraint_legalization(right_previous, legalization_queue);
}

void Triangulator::legalize_unconstrained_edges(
    std::vector<std::uint32_t>& legalization_queue) {
    // Lawson flips restore local Delaunay legality without ever removing a
    // recovered segment. A strict in-circle test leaves cocircular choices
    // deterministic.
    while (!legalization_queue.empty()) {
        const std::uint32_t edge = legalization_queue.back();
        legalization_queue.pop_back();
        const std::uint32_t pair = edge & ~1U;
        edge_constrained_[pair] &=
            static_cast<std::uint8_t>(~kLegalizationQueuedBit);
        edge_constrained_[pair + 1U] &=
            static_cast<std::uint8_t>(~kLegalizationQueuedBit);
        if (!can_flip(edge)) {
            continue;
        }
        std::uint32_t left = 0;
        std::uint32_t right = 0;
        if (!left_face_opposite(edge, left) ||
            !left_face_opposite(sym(edge), right) ||
            !active_in_circle(org(edge), dest(edge), left, right)) {
            continue;
        }

        flip_edge(edge);
        const std::uint32_t left_next = lnext(edge);
        const std::uint32_t left_previous = lnext(left_next);
        const std::uint32_t right_next = lnext(sym(edge));
        const std::uint32_t right_previous = lnext(right_next);
        queue_constraint_legalization(left_next, legalization_queue);
        queue_constraint_legalization(left_previous, legalization_queue);
        queue_constraint_legalization(right_next, legalization_queue);
        queue_constraint_legalization(right_previous, legalization_queue);
    }
}

void Triangulator::recover_constraints(
    const std::vector<Constraint>& constraints) {
    std::vector<std::pair<std::uint32_t, std::uint32_t>> normalized;
    normalized.reserve(constraints.size());
    for (const Constraint constraint : constraints) {
        if (constraint.i0 >= original_to_site_.size() ||
            constraint.i1 >= original_to_site_.size()) {
            throw std::invalid_argument(
                "constraint endpoint is outside the point array");
        }
        std::uint32_t a = original_to_site_[constraint.i0];
        std::uint32_t b = original_to_site_[constraint.i1];
        if (a == b) {
            throw std::invalid_argument(
                "constraint endpoints are coincident");
        }
        if (points_[b].original < points_[a].original) {
            std::swap(a, b);
        }
        normalized.emplace_back(a, b);
    }
    std::sort(
        normalized.begin(), normalized.end(),
        [&](const auto& lhs, const auto& rhs) {
            const auto lhs_key = std::pair<std::uint32_t, std::uint32_t>{
                points_[lhs.first].original,
                points_[lhs.second].original,
            };
            const auto rhs_key = std::pair<std::uint32_t, std::uint32_t>{
                points_[rhs.first].original,
                points_[rhs.second].original,
            };
            return lhs_key < rhs_key;
        });
    normalized.erase(
        std::unique(normalized.begin(), normalized.end()),
        normalized.end());
    if (normalized.empty()) {
        return;
    }

    // Only requested endpoints need an outgoing dart initially. Collinear
    // sites encountered while walking a segment acquire one lazily.
    std::vector<std::uint8_t> endpoint_required(points_.size(), 0);
    std::size_t missing_endpoints = 0;
    const auto require_site_edge = [&](std::uint32_t site) {
        if (endpoint_required[site] == 0) {
            endpoint_required[site] = 1;
            ++missing_endpoints;
        }
    };
    for (const auto constraint : normalized) {
        require_site_edge(constraint.first);
        require_site_edge(constraint.second);
    }

    const std::size_t endpoint_workers = std::min(
        active_thread_count_, edge_ranges_.size());
    if (endpoint_workers > 1) {
        // Workers keep the shared maps read-only. An endpoint can occur in
        // several edge ranges, so each worker publishes local candidates and
        // the caller resolves duplicates after the scan.
        std::vector<std::vector<std::pair<std::uint32_t, std::uint32_t>>>
            endpoint_hits(endpoint_workers);
        const auto scan = [&](std::size_t worker) {
            std::vector<std::pair<std::uint32_t, std::uint32_t>>& hits =
                endpoint_hits[worker];
            const std::size_t first_range =
                edge_ranges_.size() * worker / endpoint_workers;
            const std::size_t last_range =
                edge_ranges_.size() * (worker + 1) / endpoint_workers;
            for (std::size_t index = first_range;
                 index < last_range;
                 ++index) {
                const EdgeRange range = edge_ranges_[index];
                for (std::uint32_t dart = range.first;
                     dart < range.last;
                     ++dart) {
                    const std::uint32_t origin = edge_origin_[dart];
                    if (origin != kDeletedEdge &&
                        endpoint_required[origin] != 0) {
                        hits.emplace_back(origin, dart);
                    }
                }
            }
        };
        worker_team_->run(endpoint_workers, scan);
        for (const auto& hits : endpoint_hits) {
            for (const auto hit : hits) {
                if (site_edge_[hit.first] == kDeletedEdge) {
                    site_edge_[hit.first] = hit.second;
                    --missing_endpoints;
                }
            }
        }
    } else {
        for (const EdgeRange range : edge_ranges_) {
            for (std::uint32_t dart = range.first;
                 dart < range.last;
                 ++dart) {
                const std::uint32_t origin = edge_origin_[dart];
                if (origin == kDeletedEdge ||
                    endpoint_required[origin] == 0 ||
                    site_edge_[origin] != kDeletedEdge) {
                    continue;
                }
                site_edge_[origin] = dart;
                if (--missing_endpoints == 0) {
                    break;
                }
            }
            if (missing_endpoints == 0) {
                break;
            }
        }
    }
    if (missing_endpoints != 0) {
        throw std::logic_error(
            "constructed topology contains an isolated constraint endpoint");
    }

    std::vector<std::uint32_t> legalization_queue;
    for (const auto constraint : normalized) {
        recover_constraint(
            constraint.first,
            constraint.second,
            legalization_queue);
    }
    legalize_unconstrained_edges(legalization_queue);
}

}  // namespace delaunay32
