// SPDX-License-Identifier: MIT

#include "delaunay32/extras/sampling.hpp"
#include "internal.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

namespace delaunay32::extras {
namespace {

std::uint64_t point_key(const Point& point) {
    return (static_cast<std::uint64_t>(
                static_cast<std::uint32_t>(point.x))
            << 32U) |
           static_cast<std::uint32_t>(point.y);
}

void validate(const IntBounds& bounds) {
    if (bounds.min_x > bounds.max_x || bounds.min_y > bounds.max_y) {
        throw std::invalid_argument(
            "integer sampling bounds must be ordered");
    }
}

void validate(const FloatBounds& bounds) {
    if (!std::isfinite(bounds.min_x) ||
        !std::isfinite(bounds.max_x) ||
        !std::isfinite(bounds.min_y) ||
        !std::isfinite(bounds.max_y) ||
        bounds.min_x > bounds.max_x || bounds.min_y > bounds.max_y) {
        throw std::invalid_argument(
            "float sampling bounds must be finite and ordered");
    }
}

std::vector<Point> integer_corners(const IntBounds& bounds) {
    const Point candidates[] = {
        {bounds.min_x, bounds.min_y},
        {bounds.max_x, bounds.min_y},
        {bounds.max_x, bounds.max_y},
        {bounds.min_x, bounds.max_y},
    };
    std::vector<Point> corners;
    for (const Point& candidate : candidates) {
        if (std::none_of(
                corners.begin(),
                corners.end(),
                [&](const Point& point) {
                    return point.x == candidate.x &&
                           point.y == candidate.y;
                })) {
            corners.push_back(candidate);
        }
    }
    return corners;
}

std::vector<FloatPoint> float_corners(const FloatBounds& bounds) {
    const FloatPoint candidates[] = {
        {bounds.min_x, bounds.min_y},
        {bounds.max_x, bounds.min_y},
        {bounds.max_x, bounds.max_y},
        {bounds.min_x, bounds.max_y},
    };
    std::vector<FloatPoint> corners;
    for (const FloatPoint& candidate : candidates) {
        if (std::none_of(
                corners.begin(),
                corners.end(),
                [&](const FloatPoint& point) {
                    return point.x == candidate.x &&
                           point.y == candidate.y;
                })) {
            corners.push_back(candidate);
        }
    }
    return corners;
}

bool grid_can_hold(
    std::uint64_t width,
    std::uint64_t height,
    std::size_t point_count) {
    const std::uint64_t count = static_cast<std::uint64_t>(point_count);
    const std::uint64_t rows =
        count / width + (count % width == 0U ? 0U : 1U);
    return rows <= height;
}

std::uint64_t point_offset(
    const Point& point,
    const IntBounds& bounds,
    std::uint64_t width) {
    const std::uint64_t x = static_cast<std::uint64_t>(
        static_cast<std::int64_t>(point.x) - bounds.min_x);
    const std::uint64_t y = static_cast<std::uint64_t>(
        static_cast<std::int64_t>(point.y) - bounds.min_y);
    return y * width + x;
}

Point point_from_offset(
    std::uint64_t offset,
    const IntBounds& bounds,
    std::uint64_t width) {
    const std::uint64_t x = offset % width;
    const std::uint64_t y = offset / width;
    return {
        static_cast<std::int32_t>(
            static_cast<std::int64_t>(bounds.min_x) +
            static_cast<std::int64_t>(x)),
        static_cast<std::int32_t>(
            static_cast<std::int64_t>(bounds.min_y) +
            static_cast<std::int64_t>(y)),
    };
}

std::vector<Point> sample_integer_grid_without_replacement(
    const UniformIntOptions& options,
    const std::vector<Point>& corners,
    std::uint64_t width,
    std::uint64_t capacity) {
    std::vector<std::uint64_t> corner_offsets;
    corner_offsets.reserve(corners.size());
    for (const Point& corner : corners) {
        corner_offsets.push_back(
            point_offset(corner, options.bounds, width));
    }
    std::sort(corner_offsets.begin(), corner_offsets.end());

    const std::size_t random_count = options.point_count - corners.size();
    const std::uint64_t available =
        capacity - static_cast<std::uint64_t>(corners.size());
    std::mt19937_64 random(options.seed);
    std::unordered_set<std::uint64_t> selected;
    selected.reserve(random_count);

    // Floyd's algorithm draws a uniform subset in O(point_count) time even
    // when the requested sample occupies most of a small integer grid.
    const std::uint64_t start =
        available - static_cast<std::uint64_t>(random_count);
    for (std::size_t i = 0; i < random_count; ++i) {
        const std::uint64_t upper = start + static_cast<std::uint64_t>(i);
        std::uniform_int_distribution<std::uint64_t> offset(0, upper);
        const std::uint64_t candidate = offset(random);
        if (!selected.insert(candidate).second) {
            selected.insert(upper);
        }
    }

    std::vector<Point> points = corners;
    points.reserve(options.point_count);
    for (const std::uint64_t compressed_offset : selected) {
        std::uint64_t offset = compressed_offset;
        for (const std::uint64_t corner_offset : corner_offsets) {
            if (offset < corner_offset) {
                break;
            }
            ++offset;
        }
        points.push_back(point_from_offset(offset, options.bounds, width));
    }
    return points;
}

long double squared_distance_to_segment(
    const Point& point,
    const Point& a,
    const Point& b) {
    const long double dx = static_cast<long double>(b.x) - a.x;
    const long double dy = static_cast<long double>(b.y) - a.y;
    const long double length_squared = dx * dx + dy * dy;
    if (length_squared == 0.0L) {
        const long double point_dx =
            static_cast<long double>(point.x) - a.x;
        const long double point_dy =
            static_cast<long double>(point.y) - a.y;
        return point_dx * point_dx + point_dy * point_dy;
    }
    const long double projection = std::clamp(
        ((static_cast<long double>(point.x) - a.x) * dx +
         (static_cast<long double>(point.y) - a.y) * dy) /
            length_squared,
        0.0L,
        1.0L);
    const long double nearest_x = a.x + projection * dx;
    const long double nearest_y = a.y + projection * dy;
    const long double point_dx = point.x - nearest_x;
    const long double point_dy = point.y - nearest_y;
    return point_dx * point_dx + point_dy * point_dy;
}

long double minimum_boundary_distance_squared(
    const Point& point,
    const PolygonDomain& domain,
    const std::vector<Point>& points) {
    long double minimum = std::numeric_limits<long double>::max();
    const auto inspect_ring = [&](const std::vector<std::uint32_t>& ring) {
        for (std::size_t i = 0; i < ring.size(); ++i) {
            minimum = std::min(
                minimum,
                squared_distance_to_segment(
                    point,
                    points[ring[i]],
                    points[ring[(i + 1) % ring.size()]]));
        }
    };
    inspect_ring(domain.outer_ring);
    for (const std::vector<std::uint32_t>& hole : domain.holes) {
        inspect_ring(hole);
    }
    return minimum;
}

struct Candidate {
    Point point;
    std::size_t domain = 0;
};

}  // namespace

std::vector<Point> generate_uniform_int_points(
    const UniformIntOptions& options) {
    validate(options.bounds);
    const std::uint64_t width =
        static_cast<std::uint64_t>(
            static_cast<std::int64_t>(options.bounds.max_x) -
            options.bounds.min_x) +
        1U;
    const std::uint64_t height =
        static_cast<std::uint64_t>(
            static_cast<std::int64_t>(options.bounds.max_y) -
            options.bounds.min_y) +
        1U;
    if (!grid_can_hold(width, height, options.point_count)) {
        throw std::invalid_argument(
            "point_count exceeds the number of integer coordinates in bounds");
    }

    std::vector<Point> corners;
    if (options.include_corners) {
        corners = integer_corners(options.bounds);
        if (options.point_count < corners.size()) {
            throw std::invalid_argument(
                "point_count is smaller than the number of distinct corners");
        }
    }

    if (width <= std::numeric_limits<std::uint64_t>::max() / height) {
        return sample_integer_grid_without_replacement(
            options, corners, width, width * height);
    }

    // The only signed-32-bit rectangle whose cell count does not fit uint64_t
    // is the complete 2^32 by 2^32 grid. Any representable vector is sparse
    // in that domain, so ordinary rejection sampling remains efficient.
    std::vector<Point> points = corners;
    points.reserve(options.point_count);

    std::unordered_set<std::uint64_t> occupied;
    occupied.reserve(options.point_count);
    for (const Point& point : points) {
        occupied.insert(point_key(point));
    }

    std::mt19937_64 random(options.seed);
    std::uniform_int_distribution<std::int32_t> x_coordinate(
        options.bounds.min_x, options.bounds.max_x);
    std::uniform_int_distribution<std::int32_t> y_coordinate(
        options.bounds.min_y, options.bounds.max_y);
    while (points.size() < options.point_count) {
        const Point point = {x_coordinate(random), y_coordinate(random)};
        if (occupied.insert(point_key(point)).second) {
            points.push_back(point);
        }
    }
    return points;
}

std::vector<FloatPoint> generate_uniform_float_points(
    const UniformFloatOptions& options) {
    validate(options.bounds);
    std::vector<FloatPoint> points;
    if (options.include_corners) {
        points = float_corners(options.bounds);
        if (options.point_count < points.size()) {
            throw std::invalid_argument(
                "point_count is smaller than the number of distinct corners");
        }
    }
    points.reserve(options.point_count);

    std::mt19937_64 random(options.seed);
    std::uniform_real_distribution<float> x_coordinate(
        options.bounds.min_x, options.bounds.max_x);
    std::uniform_real_distribution<float> y_coordinate(
        options.bounds.min_y, options.bounds.max_y);
    while (points.size() < options.point_count) {
        points.push_back({x_coordinate(random), y_coordinate(random)});
    }
    return points;
}

std::vector<Point> sample_polygon_interiors(
    const std::vector<Point>& points,
    const std::vector<PolygonDomain>& domains,
    const BestCandidateOptions& options) {
    if (options.point_count == 0) {
        return {};
    }
    if (points.empty() || domains.empty()) {
        throw std::invalid_argument(
            "polygon sampling requires points and at least one domain");
    }
    if (options.candidates_per_point == 0 ||
        options.attempts_per_candidate == 0) {
        throw std::invalid_argument(
            "candidate and attempt counts must be positive");
    }
    for (std::size_t i = 0; i < domains.size(); ++i) {
        detail::validate_domain(
            domains[i], points.size(), "domain " + std::to_string(i));
    }

    std::int32_t minimum_x = std::numeric_limits<std::int32_t>::max();
    std::int32_t minimum_y = std::numeric_limits<std::int32_t>::max();
    std::int32_t maximum_x = std::numeric_limits<std::int32_t>::min();
    std::int32_t maximum_y = std::numeric_limits<std::int32_t>::min();
    const auto inspect_ring = [&](const std::vector<std::uint32_t>& ring) {
        for (const std::uint32_t index : ring) {
            minimum_x = std::min(minimum_x, points[index].x);
            minimum_y = std::min(minimum_y, points[index].y);
            maximum_x = std::max(maximum_x, points[index].x);
            maximum_y = std::max(maximum_y, points[index].y);
        }
    };
    for (const PolygonDomain& domain : domains) {
        inspect_ring(domain.outer_ring);
        for (const std::vector<std::uint32_t>& hole : domain.holes) {
            inspect_ring(hole);
        }
    }

    std::unordered_set<std::uint64_t> occupied;
    if (options.point_count >
        std::numeric_limits<std::size_t>::max() - points.size()) {
        throw std::invalid_argument("requested sample count is too large");
    }
    occupied.reserve(points.size() + options.point_count);
    for (const Point& point : points) {
        occupied.insert(point_key(point));
    }

    std::mt19937_64 random(options.seed);
    std::uniform_int_distribution<std::int32_t> x_coordinate(
        minimum_x, maximum_x);
    std::uniform_int_distribution<std::int32_t> y_coordinate(
        minimum_y, maximum_y);

    std::vector<Point> generated;
    generated.reserve(options.point_count);
    std::vector<std::vector<Point>> accepted_by_domain(domains.size());
    while (generated.size() < options.point_count) {
        Candidate best;
        long double best_distance_squared = -1.0L;
        for (std::size_t sample = 0;
             sample < options.candidates_per_point;
             ++sample) {
            Candidate candidate;
            bool found = false;
            for (std::size_t attempt = 0;
                 attempt < options.attempts_per_candidate;
                 ++attempt) {
                candidate.point = {
                    x_coordinate(random),
                    y_coordinate(random),
                };
                if (occupied.find(point_key(candidate.point)) !=
                    occupied.end()) {
                    continue;
                }
                for (std::size_t domain = 0;
                     domain < domains.size();
                     ++domain) {
                    if (detail::point_is_strictly_inside_domain_unchecked(
                            candidate.point, domains[domain], points)) {
                        candidate.domain = domain;
                        found = true;
                        break;
                    }
                }
                if (found) {
                    break;
                }
            }
            if (!found) {
                throw std::runtime_error(
                    "could not sample a unique point inside the domains");
            }

            long double minimum_distance_squared =
                minimum_boundary_distance_squared(
                    candidate.point,
                    domains[candidate.domain],
                    points);
            for (const Point& accepted :
                 accepted_by_domain[candidate.domain]) {
                const long double dx =
                    static_cast<long double>(candidate.point.x) - accepted.x;
                const long double dy =
                    static_cast<long double>(candidate.point.y) - accepted.y;
                minimum_distance_squared = std::min(
                    minimum_distance_squared,
                    dx * dx + dy * dy);
            }
            if (minimum_distance_squared > best_distance_squared) {
                best = candidate;
                best_distance_squared = minimum_distance_squared;
            }
        }

        occupied.insert(point_key(best.point));
        accepted_by_domain[best.domain].push_back(best.point);
        generated.push_back(best.point);
    }
    return generated;
}

}  // namespace delaunay32::extras
