// SPDX-License-Identifier: MIT

#include "triangulator_impl.hpp"

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

struct EdgeBounds {
    std::int32_t min_x;
    std::int32_t max_x;
    std::int32_t min_y;
    std::int32_t max_y;

    bool overlaps(const EdgeBounds& other) const {
        return min_x <= other.max_x && other.min_x <= max_x &&
               min_y <= other.max_y && other.min_y <= max_y;
    }
};

template <typename Sites>
EdgeBounds edge_bounds(const Sites& sites, SiteIndex a, SiteIndex b) {
    return {
        std::min(sites[a].x, sites[b].x),
        std::max(sites[a].x, sites[b].x),
        std::min(sites[a].y, sites[b].y),
        std::max(sites[a].y, sites[b].y),
    };
}

constexpr std::size_t kIndexedRingSize = 64;

// A balanced bounding-box tree only rejects disjoint edge pairs. Candidate
// pairs still use the exact predicates below, including endpoint contacts.
class RingEdgeIndex {
public:
    template <typename Sites>
    RingEdgeIndex(const Sites& sites, const Ring& ring) {
        edges_.reserve(ring.size());
        for (std::size_t i = 0; i < ring.size(); ++i) {
            edges_.push_back({
                edge_bounds(sites, ring[i], ring[(i + 1) % ring.size()]), i});
        }
        nodes_.reserve(ring.size());
        build(0, edges_.size());
    }

    template <typename Predicate>
    bool any_overlap(
        const EdgeBounds& bounds,
        const Predicate& predicate) const {
        return query(0, bounds, predicate);
    }

private:
    static constexpr std::size_t kLeafSize = 8;
    struct Edge {
        EdgeBounds bounds;
        std::size_t index;
    };
    struct Node {
        EdgeBounds bounds;
        std::size_t first;
        std::size_t last;
        std::size_t right = 0;
    };
    std::vector<Edge> edges_;
    std::vector<Node> nodes_;

    std::size_t build(std::size_t first, std::size_t last) {
        EdgeBounds bounds = edges_[first].bounds;
        for (std::size_t i = first + 1; i < last; ++i) {
            const EdgeBounds& edge = edges_[i].bounds;
            bounds.min_x = std::min(bounds.min_x, edge.min_x);
            bounds.max_x = std::max(bounds.max_x, edge.max_x);
            bounds.min_y = std::min(bounds.min_y, edge.min_y);
            bounds.max_y = std::max(bounds.max_y, edge.max_y);
        }
        const std::size_t node = nodes_.size();
        nodes_.push_back({bounds, first, last});
        if (last - first > kLeafSize) {
            const bool split_x =
                static_cast<std::int64_t>(bounds.max_x) - bounds.min_x >=
                static_cast<std::int64_t>(bounds.max_y) - bounds.min_y;
            const std::size_t middle = first + (last - first) / 2;
            std::nth_element(
                edges_.begin() + static_cast<std::ptrdiff_t>(first),
                edges_.begin() + static_cast<std::ptrdiff_t>(middle),
                edges_.begin() + static_cast<std::ptrdiff_t>(last),
                [split_x](const Edge& a, const Edge& b) {
                    const auto center = [split_x](const EdgeBounds& box) {
                        return split_x
                                   ? static_cast<std::int64_t>(box.min_x) +
                                         box.max_x
                                   : static_cast<std::int64_t>(box.min_y) +
                                         box.max_y;
                    };
                    const std::int64_t ac = center(a.bounds);
                    const std::int64_t bc = center(b.bounds);
                    return ac != bc ? ac < bc : a.index < b.index;
                });
            // The left child immediately follows its parent.
            build(first, middle);
            const std::size_t right = build(middle, last);
            nodes_[node].right = right;
        }
        return node;
    }

    template <typename Predicate>
    bool query(
        std::size_t index,
        const EdgeBounds& bounds,
        const Predicate& predicate) const {
        const Node& node = nodes_[index];
        if (!node.bounds.overlaps(bounds)) {
            return false;
        }
        if (node.right != 0) {
            return query(index + 1, bounds, predicate) ||
                   query(node.right, bounds, predicate);
        }
        for (std::size_t i = node.first; i < node.last; ++i) {
            if (edges_[i].bounds.overlaps(bounds) &&
                predicate(edges_[i].index)) {
                return true;
            }
        }
        return false;
    }
};

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
    if (!edge_bounds(sites, a, b).overlaps(edge_bounds(sites, c, d))) {
        return false;
    }
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
    const auto intersects = [&](std::size_t i, std::size_t j) {
        const SiteIndex a = ring[i];
        const SiteIndex b = ring[(i + 1) % size];
        const SiteIndex c = ring[j];
        const SiteIndex d = ring[(j + 1) % size];
        const bool adjacent =
            j == i + 1 || (i == 0 && j + 1 == size);
        if (!adjacent) {
            return segments_intersect_or_touch(
                sites, orientation, a, b, c, d);
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
        return false;
    };

    if (size >= kIndexedRingSize) {
        const RingEdgeIndex index(sites, ring);
        for (std::size_t i = 0; i < size; ++i) {
            if (index.any_overlap(
                    edge_bounds(sites, ring[i], ring[(i + 1) % size]),
                    [&](std::size_t j) { return j > i && intersects(i, j); })) {
                throw std::invalid_argument("polygon ring is not simple");
            }
        }
    } else {
        for (std::size_t i = 0; i < size; ++i) {
            for (std::size_t j = i + 1; j < size; ++j) {
                if (intersects(i, j)) {
                    throw std::invalid_argument("polygon ring is not simple");
                }
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
    if (std::max(first.size(), second.size()) >= kIndexedRingSize) {
        const Ring& indexed = first.size() >= second.size() ? first : second;
        const Ring& queried = first.size() >= second.size() ? second : first;
        const RingEdgeIndex index(sites, indexed);
        for (std::size_t i = 0; i < queried.size(); ++i) {
            const SiteIndex a = queried[i];
            const SiteIndex b = queried[(i + 1) % queried.size()];
            if (index.any_overlap(edge_bounds(sites, a, b), [&](std::size_t j) {
                    return segments_intersect_or_touch(
                        sites, orientation, a, b,
                        indexed[j], indexed[(j + 1) % indexed.size()]);
                })) {
                return true;
            }
        }
        return false;
    }
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
Triangulator::Impl::prepare_polygon_rings(
    const std::vector<std::uint32_t>& outer_ring,
    const std::vector<std::vector<std::uint32_t>>& holes) const {
    Rings rings;
    rings.reserve(holes.size() + 1);
    rings.push_back(map_ring_indices(outer_ring, constraints_.original_to_site));
    for (const std::vector<std::uint32_t>& hole : holes) {
        rings.push_back(map_ring_indices(hole, constraints_.original_to_site));
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

std::vector<std::vector<std::vector<std::uint32_t>>>
Triangulator::Impl::prepare_polygon_domains() const {
    std::vector<Rings> domains;
    domains.reserve(constraints_.polygons.size());
    for (const PolygonDomain& polygon : constraints_.polygons) {
        domains.push_back(prepare_polygon_rings(
            polygon.outer_ring, polygon.holes));
    }

    const auto orientation = [&](SiteIndex a, SiteIndex b, SiteIndex c) {
        return orient(a, b, c);
    };
    for (std::size_t first = 0; first < domains.size(); ++first) {
        const Ring& first_outer = domains[first].front();
        for (std::size_t second = first + 1;
             second < domains.size();
             ++second) {
            const Ring& second_outer = domains[second].front();
            if (rings_intersect_or_touch(
                    points_,
                    orientation,
                    first_outer,
                    second_outer) ||
                point_in_ring(
                    points_,
                    orientation,
                    first_outer.front(),
                    second_outer) ||
                point_in_ring(
                    points_,
                    orientation,
                    second_outer.front(),
                    first_outer)) {
                throw std::invalid_argument(
                    "polygon outer domains overlap, nest, or touch");
            }
        }
    }
    return domains;
}

std::uint32_t Triangulator::Impl::first_boundary_edge(
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

void Triangulator::Impl::mark_polygon_excluded_faces(
    const std::vector<Rings>& domains) {
    constexpr std::uint8_t kBoundaryBit = 1;
    constexpr std::uint8_t kExcludedBit = 2;
    std::vector<std::uint8_t> marks(constraints_.edge_flags.size(), 0);

    // Standalone constraints protect edges during recovery and legalization,
    // but only polygon boundaries stop the exclusion flood. Mark every dart
    // in each recovered boundary chain, including unlisted collinear sites.
    for (const Rings& rings : domains) {
        for (const Ring& ring : rings) {
            for (std::size_t i = 0; i < ring.size(); ++i) {
                std::uint32_t origin = ring[i];
                const std::uint32_t destination =
                    ring[(i + 1) % ring.size()];
                while (origin != destination) {
                    const std::uint32_t edge =
                        first_boundary_edge(origin, destination);
                    marks[edge] |= kBoundaryBit;
                    marks[sym(edge)] |= kBoundaryBit;
                    origin = dest(edge);
                }
            }
        }
    }

    const auto exclude_component = [&](std::uint32_t initial_seed) {
        std::deque<std::uint32_t> pending = {initial_seed};
        while (!pending.empty()) {
            const std::uint32_t seed = pending.front();
            pending.pop_front();
            if (seed >= marks.size() ||
                (marks[seed] & kExcludedBit) != 0) {
                continue;
            }

            std::uint32_t edge = seed;
            std::size_t face_size = 0;
            do {
                if (edge >= marks.size() || !is_live_edge(edge)) {
                    throw std::logic_error(
                        "polygon face flood reached invalid topology");
                }
                marks[edge] |= kExcludedBit;
                edge = lnext(edge);
                ++face_size;
                if (face_size > marks.size()) {
                    throw std::logic_error(
                        "polygon face flood did not close a face");
                }
            } while (edge != seed);

            edge = seed;
            do {
                const std::uint32_t neighbor = sym(edge);
                if ((marks[edge] & kBoundaryBit) == 0 &&
                    (marks[neighbor] & kExcludedBit) == 0) {
                    pending.push_back(neighbor);
                }
                edge = lnext(edge);
            } while (edge != seed);
        }
    };

    exclude_component(outer_seed_);
    for (const Rings& rings : domains) {
        for (std::size_t hole = 1; hole < rings.size(); ++hole) {
            const Ring& ring = rings[hole];
            // Holes are normalized clockwise, so their interior is the right
            // face of every directed boundary edge.
            exclude_component(sym(first_boundary_edge(ring[0], ring[1])));
        }
    }

    for (const EdgeRange range : arena_.ranges) {
        for (std::uint32_t edge = range.first;
             edge < range.last;
             ++edge) {
            if ((marks[edge] & kExcludedBit) != 0) {
                arena_.origin[edge] |= kVisitedBit;
            }
        }
    }
}

}  // namespace delaunay32
