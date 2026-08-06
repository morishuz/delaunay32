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
#include <utility>
#include <vector>

namespace delaunay32::extras {
namespace {

void validate_bounds(const SamplingBounds& bounds) {
    if (!std::isfinite(bounds.min_x) ||
        !std::isfinite(bounds.max_x) ||
        !std::isfinite(bounds.min_y) ||
        !std::isfinite(bounds.max_y) ||
        bounds.min_x > bounds.max_x || bounds.min_y > bounds.max_y) {
        throw std::invalid_argument(
            "sampling bounds must be finite and ordered");
    }
}

void validate_polygon_region(
    const std::vector<FloatPoint>& points,
    const std::vector<PolygonDomain>& domains) {
    if (points.empty() || domains.empty()) {
        throw std::invalid_argument(
            "polygon sampling requires points and at least one domain");
    }
    for (const FloatPoint& point : points) {
        if (!std::isfinite(point.x) || !std::isfinite(point.y)) {
            throw std::invalid_argument(
                "polygon sampling coordinates must be finite");
        }
    }
    for (std::size_t i = 0; i < domains.size(); ++i) {
        detail::validate_domain(
            domains[i], points.size(), "domain " + std::to_string(i));
    }
}

SamplingBounds polygon_bounds(
    const std::vector<FloatPoint>& points,
    const std::vector<PolygonDomain>& domains) {
    SamplingBounds bounds{
        std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity(),
    };
    const auto inspect_ring = [&](const std::vector<std::uint32_t>& ring) {
        for (const std::uint32_t index : ring) {
            bounds.min_x = std::min(bounds.min_x, points[index].x);
            bounds.max_x = std::max(bounds.max_x, points[index].x);
            bounds.min_y = std::min(bounds.min_y, points[index].y);
            bounds.max_y = std::max(bounds.max_y, points[index].y);
        }
    };
    for (const PolygonDomain& domain : domains) {
        inspect_ring(domain.outer_ring);
        for (const std::vector<std::uint32_t>& hole : domain.holes) {
            inspect_ring(hole);
        }
    }
    return bounds;
}

std::vector<FloatPoint> distinct_corners(const SamplingBounds& bounds) {
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

struct PolygonCandidate {
    FloatPoint point;
    std::size_t domain = 0;
};

PolygonCandidate polygon_candidate(
    std::mt19937_64& random,
    const SamplingBounds& bounds,
    const std::vector<FloatPoint>& points,
    const std::vector<PolygonDomain>& domains,
    std::size_t attempts) {
    std::uniform_real_distribution<double> x_coordinate(
        bounds.min_x, bounds.max_x);
    std::uniform_real_distribution<double> y_coordinate(
        bounds.min_y, bounds.max_y);
    for (std::size_t attempt = 0; attempt < attempts; ++attempt) {
        const FloatPoint candidate = {
            x_coordinate(random),
            y_coordinate(random),
        };
        for (std::size_t domain = 0; domain < domains.size(); ++domain) {
            if (detail::point_is_strictly_inside_domain_unchecked(
                    candidate, domains[domain], points)) {
                return {candidate, domain};
            }
        }
    }
    throw std::runtime_error(
        "could not sample a point inside the configured polygon domains");
}

long double squared_distance(
    const FloatPoint& first,
    const FloatPoint& second) {
    const long double dx =
        static_cast<long double>(first.x) - second.x;
    const long double dy =
        static_cast<long double>(first.y) - second.y;
    return dx * dx + dy * dy;
}

long double squared_distance_to_segment(
    const FloatPoint& point,
    const FloatPoint& a,
    const FloatPoint& b) {
    const long double dx = static_cast<long double>(b.x) - a.x;
    const long double dy = static_cast<long double>(b.y) - a.y;
    const long double length_squared = dx * dx + dy * dy;
    if (length_squared == 0.0L) {
        return squared_distance(point, a);
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

long double polygon_boundary_distance_squared(
    const FloatPoint& point,
    const PolygonDomain& domain,
    const std::vector<FloatPoint>& points) {
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

long double bounds_boundary_distance_squared(
    const FloatPoint& point,
    const SamplingBounds& bounds) {
    const long double distance = std::min({
        static_cast<long double>(point.x - bounds.min_x),
        static_cast<long double>(bounds.max_x - point.x),
        static_cast<long double>(point.y - bounds.min_y),
        static_cast<long double>(bounds.max_y - point.y),
    });
    return distance * distance;
}

void include_sample_distances(
    long double& minimum_distance_squared,
    const FloatPoint& candidate,
    const std::vector<FloatPoint>& accepted) {
    for (const FloatPoint& point : accepted) {
        minimum_distance_squared = std::min(
            minimum_distance_squared,
            squared_distance(candidate, point));
    }
}

}  // namespace

void PointSampler::set_bounds(SamplingBounds bounds) {
    validate_bounds(bounds);
    bounds_ = bounds;
    polygon_points_.clear();
    domains_.clear();
    region_ = Region::Bounds;
}

void PointSampler::set_polygon_interiors(
    std::vector<Point> points,
    std::vector<PolygonDomain> domains) {
    std::vector<FloatPoint> converted;
    converted.reserve(points.size());
    for (const Point& point : points) {
        converted.push_back({
            static_cast<double>(point.x),
            static_cast<double>(point.y),
        });
    }
    set_polygon_interiors(std::move(converted), std::move(domains));
}

void PointSampler::set_polygon_interiors(
    std::vector<FloatPoint> points,
    std::vector<PolygonDomain> domains) {
    validate_polygon_region(points, domains);
    bounds_ = polygon_bounds(points, domains);
    polygon_points_ = std::move(points);
    domains_ = std::move(domains);
    region_ = Region::PolygonInteriors;
}

std::vector<FloatPoint> PointSampler::generate_uniform(
    const UniformSamplingOptions& options) const {
    if (region_ == Region::None) {
        throw std::logic_error(
            "point sampler requires configured bounds or polygon interiors");
    }
    if (options.attempts_per_point == 0) {
        throw std::invalid_argument(
            "uniform sampling attempts_per_point must be positive");
    }
    if (region_ == Region::PolygonInteriors &&
        options.include_bounds_corners) {
        throw std::invalid_argument(
            "include_bounds_corners is only valid for bounds sampling");
    }

    std::vector<FloatPoint> generated;
    if (region_ == Region::Bounds && options.include_bounds_corners) {
        generated = distinct_corners(bounds_);
        if (options.point_count < generated.size()) {
            throw std::invalid_argument(
                "point_count is smaller than the number of distinct bounds "
                "corners");
        }
    }
    generated.reserve(options.point_count);

    std::mt19937_64 random(options.seed);
    if (region_ == Region::Bounds) {
        std::uniform_real_distribution<double> x_coordinate(
            bounds_.min_x, bounds_.max_x);
        std::uniform_real_distribution<double> y_coordinate(
            bounds_.min_y, bounds_.max_y);
        while (generated.size() < options.point_count) {
            generated.push_back({x_coordinate(random), y_coordinate(random)});
        }
        return generated;
    }

    while (generated.size() < options.point_count) {
        generated.push_back(
            polygon_candidate(
                random,
                bounds_,
                polygon_points_,
                domains_,
                options.attempts_per_point)
                .point);
    }
    return generated;
}

std::vector<FloatPoint> PointSampler::generate_blue_noise(
    const BlueNoiseSamplingOptions& options) const {
    if (region_ == Region::None) {
        throw std::logic_error(
            "point sampler requires configured bounds or polygon interiors");
    }
    if (options.candidates_per_point == 0 ||
        options.attempts_per_candidate == 0) {
        throw std::invalid_argument(
            "blue-noise candidate and attempt counts must be positive");
    }

    std::mt19937_64 random(options.seed);
    std::uniform_real_distribution<double> x_coordinate(
        bounds_.min_x, bounds_.max_x);
    std::uniform_real_distribution<double> y_coordinate(
        bounds_.min_y, bounds_.max_y);
    std::vector<FloatPoint> generated;
    generated.reserve(options.point_count);
    while (generated.size() < options.point_count) {
        FloatPoint best;
        long double best_distance_squared = -1.0L;
        for (std::size_t sample = 0;
             sample < options.candidates_per_point;
             ++sample) {
            FloatPoint candidate;
            long double minimum_distance_squared = 0.0L;
            if (region_ == Region::Bounds) {
                candidate = {x_coordinate(random), y_coordinate(random)};
                minimum_distance_squared =
                    bounds_boundary_distance_squared(candidate, bounds_);
            } else {
                const PolygonCandidate polygon = polygon_candidate(
                    random,
                    bounds_,
                    polygon_points_,
                    domains_,
                    options.attempts_per_candidate);
                candidate = polygon.point;
                minimum_distance_squared =
                    polygon_boundary_distance_squared(
                        candidate,
                        domains_[polygon.domain],
                        polygon_points_);
            }
            include_sample_distances(
                minimum_distance_squared, candidate, generated);
            if (minimum_distance_squared > best_distance_squared) {
                best = candidate;
                best_distance_squared = minimum_distance_squared;
            }
        }
        generated.push_back(best);
    }
    return generated;
}

}  // namespace delaunay32::extras
