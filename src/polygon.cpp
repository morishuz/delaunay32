// SPDX-License-Identifier: MIT

#include "delaunay32/delaunay.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <iterator>
#include <stdexcept>
#include <unordered_set>
#include <vector>

namespace delaunay32 {

std::vector<std::vector<std::uint32_t>>
Triangulator::prepare_polygon_rings(
    const std::vector<std::uint32_t>& outer_ring,
    const std::vector<std::vector<std::uint32_t>>& holes) const {
    std::vector<std::vector<std::uint32_t>> rings;
    rings.reserve(holes.size() + 1);

    const auto map_ring = [&](const std::vector<std::uint32_t>& input) {
        if (input.size() < 3) {
            throw std::invalid_argument(
                "polygon rings need at least three indices");
        }
        std::vector<std::uint32_t> ring;
        ring.reserve(input.size());
        for (const std::uint32_t index : input) {
            if (index >= original_to_site_.size()) {
                throw std::invalid_argument(
                    "polygon ring index is outside the point array");
            }
            ring.push_back(original_to_site_[index]);
        }
        if (ring.size() > 1 && ring.front() == ring.back()) {
            ring.pop_back();
        }
        if (ring.size() < 3) {
            throw std::invalid_argument(
                "polygon rings need three distinct coordinates");
        }
        std::unordered_set<std::uint32_t> unique;
        unique.reserve(ring.size());
        for (const std::uint32_t site : ring) {
            if (!unique.insert(site).second) {
                throw std::invalid_argument(
                    "polygon ring repeats a coordinate");
            }
        }
        return ring;
    };

    rings.push_back(map_ring(outer_ring));
    for (const std::vector<std::uint32_t>& hole : holes) {
        rings.push_back(map_ring(hole));
    }

    const auto point_on_segment = [&](std::uint32_t point,
                                      std::uint32_t a,
                                      std::uint32_t b) {
        if (orient(a, b, point) != 0) {
            return false;
        }
        const Site& p = points_[point];
        const Site& u = points_[a];
        const Site& v = points_[b];
        return p.x >= std::min(u.x, v.x) &&
               p.x <= std::max(u.x, v.x) &&
               p.y >= std::min(u.y, v.y) &&
               p.y <= std::max(u.y, v.y);
    };
    const auto segments_intersect = [&](std::uint32_t a,
                                        std::uint32_t b,
                                        std::uint32_t c,
                                        std::uint32_t d) {
        const std::int64_t ab_c = orient(a, b, c);
        const std::int64_t ab_d = orient(a, b, d);
        const std::int64_t cd_a = orient(c, d, a);
        const std::int64_t cd_b = orient(c, d, b);
        const auto opposite = [](std::int64_t lhs, std::int64_t rhs) {
            return (lhs > 0 && rhs < 0) ||
                   (lhs < 0 && rhs > 0);
        };
        return (opposite(ab_c, ab_d) && opposite(cd_a, cd_b)) ||
               (ab_c == 0 && point_on_segment(c, a, b)) ||
               (ab_d == 0 && point_on_segment(d, a, b)) ||
               (cd_a == 0 && point_on_segment(a, c, d)) ||
               (cd_b == 0 && point_on_segment(b, c, d));
    };

    for (std::vector<std::uint32_t>& ring : rings) {
        const std::size_t size = ring.size();
        for (std::size_t i = 0; i < size; ++i) {
            const std::uint32_t a = ring[i];
            const std::uint32_t b = ring[(i + 1) % size];
            for (std::size_t j = i + 1; j < size; ++j) {
                const std::uint32_t c = ring[j];
                const std::uint32_t d = ring[(j + 1) % size];
                const bool adjacent =
                    j == i + 1 || (i == 0 && j + 1 == size);
                if (!adjacent) {
                    if (segments_intersect(a, b, c, d)) {
                        throw std::invalid_argument(
                            "polygon ring is not simple");
                    }
                    continue;
                }

                const std::uint32_t shared =
                    a == c || a == d ? a : b;
                const std::uint32_t first_other =
                    a == shared ? b : a;
                const std::uint32_t second_other =
                    c == shared ? d : c;
                if (point_on_segment(
                        first_other, shared, second_other) ||
                    point_on_segment(
                        second_other, shared, first_other)) {
                    throw std::invalid_argument(
                        "polygon ring has overlapping adjacent edges");
                }
            }
        }

        const auto extreme = std::min_element(
            ring.begin(), ring.end(), [&](std::uint32_t lhs,
                                          std::uint32_t rhs) {
                const Site& a = points_[lhs];
                const Site& b = points_[rhs];
                return a.x != b.x ? a.x < b.x : a.y < b.y;
            });
        const std::size_t index = static_cast<std::size_t>(
            std::distance(ring.begin(), extreme));
        const std::int64_t winding = orient(
            ring[(index + size - 1) % size],
            ring[index],
            ring[(index + 1) % size]);
        if (winding == 0) {
            throw std::invalid_argument(
                "polygon ring has no well-defined winding");
        }
        const bool want_counterclockwise = &ring == &rings.front();
        if ((winding > 0) != want_counterclockwise) {
            std::reverse(ring.begin(), ring.end());
        }
    }

    for (std::size_t first = 0; first < rings.size(); ++first) {
        const std::vector<std::uint32_t>& a = rings[first];
        for (std::size_t second = first + 1;
             second < rings.size();
             ++second) {
            const std::vector<std::uint32_t>& b = rings[second];
            for (std::size_t i = 0; i < a.size(); ++i) {
                for (std::size_t j = 0; j < b.size(); ++j) {
                    if (segments_intersect(
                            a[i],
                            a[(i + 1) % a.size()],
                            b[j],
                            b[(j + 1) % b.size()])) {
                        throw std::invalid_argument(
                            "polygon rings intersect or touch");
                    }
                }
            }
        }
    }

    const auto point_in_ring = [&](std::uint32_t point,
                                   const std::vector<std::uint32_t>& ring) {
        bool inside = false;
        const Site& p = points_[point];
        for (std::size_t i = 0; i < ring.size(); ++i) {
            const std::uint32_t a = ring[i];
            const std::uint32_t b = ring[(i + 1) % ring.size()];
            const Site& u = points_[a];
            const Site& v = points_[b];
            if ((u.y <= p.y && p.y < v.y && orient(a, b, point) > 0) ||
                (v.y <= p.y && p.y < u.y && orient(a, b, point) < 0)) {
                inside = !inside;
            }
        }
        return inside;
    };

    for (std::size_t hole = 1; hole < rings.size(); ++hole) {
        if (!point_in_ring(rings[hole].front(), rings.front())) {
            throw std::invalid_argument(
                "polygon hole is not strictly inside the outer ring");
        }
    }
    for (std::size_t first = 1; first < rings.size(); ++first) {
        for (std::size_t second = first + 1;
             second < rings.size();
             ++second) {
            if (point_in_ring(rings[first].front(), rings[second]) ||
                point_in_ring(rings[second].front(), rings[first])) {
                throw std::invalid_argument(
                    "polygon holes overlap or are nested");
            }
        }
    }
    return rings;
}

std::uint32_t Triangulator::first_boundary_edge(
    std::uint32_t origin,
    std::uint32_t destination) const {
    std::uint32_t edge = find_edge(origin, destination);
    if (edge == kDeletedEdge) {
        edge = first_collinear_edge(origin, destination);
    }
    if (edge == kDeletedEdge || !is_constrained(edge)) {
        throw std::logic_error(
            "polygon boundary recovery did not create an edge chain");
    }
    return edge;
}

void Triangulator::mark_polygon_excluded_faces(
    const std::vector<std::vector<std::uint32_t>>& rings) {
    std::vector<std::uint8_t> excluded(edge_constrained_.size(), 0);

    const auto exclude_component = [&](std::uint32_t initial_seed) {
        std::deque<std::uint32_t> pending = {initial_seed};
        while (!pending.empty()) {
            const std::uint32_t seed = pending.front();
            pending.pop_front();
            if (seed >= excluded.size() || excluded[seed] != 0) {
                continue;
            }

            std::uint32_t edge = seed;
            std::size_t face_size = 0;
            do {
                if (edge >= excluded.size() || !is_live_edge(edge)) {
                    throw std::logic_error(
                        "polygon face flood reached invalid topology");
                }
                excluded[edge] = 1;
                edge = lnext(edge);
                ++face_size;
                if (face_size > excluded.size()) {
                    throw std::logic_error(
                        "polygon face flood did not close a face");
                }
            } while (edge != seed);

            edge = seed;
            do {
                const std::uint32_t neighbor = sym(edge);
                if (!is_constrained(edge) &&
                    excluded[neighbor] == 0) {
                    pending.push_back(neighbor);
                }
                edge = lnext(edge);
            } while (edge != seed);
        }
    };

    exclude_component(outer_seed_);
    for (std::size_t hole = 1; hole < rings.size(); ++hole) {
        const std::vector<std::uint32_t>& ring = rings[hole];
        // Holes are normalized clockwise, so their interior is the right face
        // of every directed boundary edge.
        exclude_component(sym(first_boundary_edge(ring[0], ring[1])));
    }

    for (const EdgeRange range : edge_ranges_) {
        for (std::uint32_t edge = range.first;
             edge < range.last;
             ++edge) {
            if (excluded[edge] != 0) {
                edge_origin_[edge] |= kVisitedBit;
            }
        }
    }
}

std::vector<Triangle> Triangulator::triangulate_polygon_int(
    const std::vector<Point>& points,
    const std::vector<std::uint32_t>& outer_ring,
    const std::vector<std::vector<std::uint32_t>>& holes) {
    require_point_count(points.size());
    load_int_points(points);
    std::vector<std::uint32_t> representatives;
    build_loaded_topology(&representatives);
    build_constraint_indices(representatives, points.size());
    const std::vector<std::vector<std::uint32_t>> rings =
        prepare_polygon_rings(outer_ring, holes);

    std::vector<Constraint> boundaries;
    std::size_t boundary_count = 0;
    for (const std::vector<std::uint32_t>& ring : rings) {
        boundary_count += ring.size();
    }
    boundaries.reserve(boundary_count);
    for (const std::vector<std::uint32_t>& ring : rings) {
        for (std::size_t i = 0; i < ring.size(); ++i) {
            boundaries.push_back({
                points_[ring[i]].original,
                points_[ring[(i + 1) % ring.size()]].original,
            });
        }
    }
    recover_constraints(boundaries);
    mark_polygon_excluded_faces(rings);

    if (active_thread_count_ > 1) {
        export_triangles_parallel(
            active_thread_count_, *worker_team_);
    } else {
        export_triangles();
    }
    return std::move(triangles_out_);
}

}  // namespace delaunay32
