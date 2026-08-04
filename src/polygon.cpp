// SPDX-License-Identifier: MIT

#include "delaunay32/delaunay.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <iterator>
#include <stdexcept>
#include <vector>

namespace delaunay32 {
namespace {

using SiteIndex = std::uint32_t;
using Ring = std::vector<SiteIndex>;
using Rings = std::vector<Ring>;

void require_ring_indices(
    const std::vector<std::uint32_t>& ring,
    std::size_t point_count) {
    if (ring.size() < 3) {
        throw std::invalid_argument(
            "polygon rings need at least three indices");
    }
    for (const std::uint32_t index : ring) {
        if (index >= point_count) {
            throw std::invalid_argument(
                "polygon ring index is outside the point array");
        }
    }
}

Ring map_ring_indices(
    const std::vector<std::uint32_t>& input,
    const std::vector<std::uint32_t>& original_to_site) {
    Ring ring;
    ring.reserve(input.size());
    for (const std::uint32_t index : input) {
        ring.push_back(original_to_site[index]);
    }
    if (ring.size() > 1 && ring.front() == ring.back()) {
        ring.pop_back();
    }
    if (ring.size() < 3) {
        throw std::invalid_argument(
            "polygon rings need three distinct coordinates");
    }

    Ring sorted = ring;
    std::sort(sorted.begin(), sorted.end());
    if (std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end()) {
        throw std::invalid_argument(
            "polygon ring repeats a coordinate");
    }
    return ring;
}

template <typename Sites, typename Orientation>
bool point_on_segment(
    const Sites& sites,
    const Orientation& orientation,
    SiteIndex point,
    SiteIndex a,
    SiteIndex b) {
    if (orientation(a, b, point) != 0) {
        return false;
    }
    const auto& p = sites[point];
    const auto& u = sites[a];
    const auto& v = sites[b];
    return p.x >= std::min(u.x, v.x) &&
           p.x <= std::max(u.x, v.x) &&
           p.y >= std::min(u.y, v.y) &&
           p.y <= std::max(u.y, v.y);
}

template <typename Sites, typename Orientation>
bool segments_intersect_or_touch(
    const Sites& sites,
    const Orientation& orientation,
    SiteIndex a,
    SiteIndex b,
    SiteIndex c,
    SiteIndex d) {
    const std::int64_t ab_c = orientation(a, b, c);
    const std::int64_t ab_d = orientation(a, b, d);
    const std::int64_t cd_a = orientation(c, d, a);
    const std::int64_t cd_b = orientation(c, d, b);
    const auto opposite = [](std::int64_t lhs, std::int64_t rhs) {
        return (lhs > 0 && rhs < 0) ||
               (lhs < 0 && rhs > 0);
    };
    return (opposite(ab_c, ab_d) && opposite(cd_a, cd_b)) ||
           (ab_c == 0 &&
            point_on_segment(sites, orientation, c, a, b)) ||
           (ab_d == 0 &&
            point_on_segment(sites, orientation, d, a, b)) ||
           (cd_a == 0 &&
            point_on_segment(sites, orientation, a, c, d)) ||
           (cd_b == 0 &&
            point_on_segment(sites, orientation, b, c, d));
}

template <typename Sites, typename Orientation>
void validate_simple_ring(
    const Sites& sites,
    const Orientation& orientation,
    const Ring& ring) {
    const std::size_t size = ring.size();
    for (std::size_t i = 0; i < size; ++i) {
        const SiteIndex a = ring[i];
        const SiteIndex b = ring[(i + 1) % size];
        for (std::size_t j = i + 1; j < size; ++j) {
            const SiteIndex c = ring[j];
            const SiteIndex d = ring[(j + 1) % size];
            const bool adjacent =
                j == i + 1 || (i == 0 && j + 1 == size);
            if (!adjacent) {
                if (segments_intersect_or_touch(
                        sites, orientation, a, b, c, d)) {
                    throw std::invalid_argument(
                        "polygon ring is not simple");
                }
                continue;
            }

            const SiteIndex shared = a == c || a == d ? a : b;
            const SiteIndex first_other = a == shared ? b : a;
            const SiteIndex second_other = c == shared ? d : c;
            if (point_on_segment(
                    sites,
                    orientation,
                    first_other,
                    shared,
                    second_other) ||
                point_on_segment(
                    sites,
                    orientation,
                    second_other,
                    shared,
                    first_other)) {
                throw std::invalid_argument(
                    "polygon ring has overlapping adjacent edges");
            }
        }
    }
}

template <typename Sites, typename Orientation>
void normalize_ring_winding(
    const Sites& sites,
    const Orientation& orientation,
    Ring& ring,
    bool want_counterclockwise) {
    const auto extreme = std::min_element(
        ring.begin(), ring.end(), [&](SiteIndex lhs, SiteIndex rhs) {
            const auto& a = sites[lhs];
            const auto& b = sites[rhs];
            return a.x != b.x ? a.x < b.x : a.y < b.y;
        });
    const std::size_t index = static_cast<std::size_t>(
        std::distance(ring.begin(), extreme));
    const std::size_t size = ring.size();
    const std::int64_t winding = orientation(
        ring[(index + size - 1) % size],
        ring[index],
        ring[(index + 1) % size]);
    if (winding == 0) {
        throw std::invalid_argument(
            "polygon ring has no well-defined winding");
    }
    if ((winding > 0) != want_counterclockwise) {
        std::reverse(ring.begin(), ring.end());
    }
}

template <typename Sites, typename Orientation>
bool rings_intersect_or_touch(
    const Sites& sites,
    const Orientation& orientation,
    const Ring& first,
    const Ring& second) {
    for (std::size_t i = 0; i < first.size(); ++i) {
        for (std::size_t j = 0; j < second.size(); ++j) {
            if (segments_intersect_or_touch(
                    sites,
                    orientation,
                    first[i],
                    first[(i + 1) % first.size()],
                    second[j],
                    second[(j + 1) % second.size()])) {
                return true;
            }
        }
    }
    return false;
}

template <typename Sites, typename Orientation>
bool point_in_ring(
    const Sites& sites,
    const Orientation& orientation,
    SiteIndex point,
    const Ring& ring) {
    bool inside = false;
    const auto& p = sites[point];
    for (std::size_t i = 0; i < ring.size(); ++i) {
        const SiteIndex a = ring[i];
        const SiteIndex b = ring[(i + 1) % ring.size()];
        const auto& u = sites[a];
        const auto& v = sites[b];
        if ((u.y <= p.y && p.y < v.y && orientation(a, b, point) > 0) ||
            (v.y <= p.y && p.y < u.y && orientation(a, b, point) < 0)) {
            inside = !inside;
        }
    }
    return inside;
}

}  // namespace

std::vector<std::vector<std::uint32_t>>
Triangulator::prepare_polygon_rings(
    const std::vector<std::uint32_t>& outer_ring,
    const std::vector<std::vector<std::uint32_t>>& holes) const {
    Rings rings;
    rings.reserve(holes.size() + 1);
    rings.push_back(map_ring_indices(outer_ring, original_to_site_));
    for (const std::vector<std::uint32_t>& hole : holes) {
        rings.push_back(map_ring_indices(hole, original_to_site_));
    }

    const auto orientation = [&](SiteIndex a, SiteIndex b, SiteIndex c) {
        return orient(a, b, c);
    };

    for (std::size_t i = 0; i < rings.size(); ++i) {
        validate_simple_ring(points_, orientation, rings[i]);
        normalize_ring_winding(
            points_, orientation, rings[i], i == 0);
    }

    for (std::size_t first = 0; first < rings.size(); ++first) {
        for (std::size_t second = first + 1;
             second < rings.size();
             ++second) {
            if (rings_intersect_or_touch(
                    points_,
                    orientation,
                    rings[first],
                    rings[second])) {
                throw std::invalid_argument(
                    "polygon rings intersect or touch");
            }
        }
    }

    for (std::size_t hole = 1; hole < rings.size(); ++hole) {
        if (!point_in_ring(
                points_,
                orientation,
                rings[hole].front(),
                rings.front())) {
            throw std::invalid_argument(
                "polygon hole is not strictly inside the outer ring");
        }
    }
    for (std::size_t first = 1; first < rings.size(); ++first) {
        for (std::size_t second = first + 1;
             second < rings.size();
             ++second) {
            if (point_in_ring(
                    points_,
                    orientation,
                    rings[first].front(),
                    rings[second]) ||
                point_in_ring(
                    points_,
                    orientation,
                    rings[second].front(),
                    rings[first])) {
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
    require_ring_indices(outer_ring, points.size());
    for (const std::vector<std::uint32_t>& hole : holes) {
        require_ring_indices(hole, points.size());
    }
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
    return finish_triangle_export();
}

}  // namespace delaunay32
