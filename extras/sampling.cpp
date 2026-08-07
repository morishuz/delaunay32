// SPDX-License-Identifier: MIT

#include "delaunay32/extras/sampling.hpp"
#include "internal.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace delaunay32::extras {
namespace {

using SamplingScalar = detail::SamplingScalar;

SamplingBounds empty_bounds() {
    return {
        std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity(),
    };
}

void expand_bounds(SamplingBounds& bounds, const FloatPoint& point) {
    bounds.min_x = std::min(bounds.min_x, point.x);
    bounds.max_x = std::max(bounds.max_x, point.x);
    bounds.min_y = std::min(bounds.min_y, point.y);
    bounds.max_y = std::max(bounds.max_y, point.y);
}

void expand_bounds(
    SamplingBounds& bounds,
    const SamplingBounds& extension) {
    bounds.min_x = std::min(bounds.min_x, extension.min_x);
    bounds.max_x = std::max(bounds.max_x, extension.max_x);
    bounds.min_y = std::min(bounds.min_y, extension.min_y);
    bounds.max_y = std::max(bounds.max_y, extension.max_y);
}

template <typename Function>
void for_each_ring(const PolygonDomain& domain, const Function& function) {
    function(domain.outer_ring);
    for (const std::vector<std::uint32_t>& hole : domain.holes) {
        function(hole);
    }
}

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

void validate_options(
    const UniformSamplingOptions& options,
    bool polygon_region) {
    if (options.attempts_per_point == 0) {
        throw std::invalid_argument(
            "uniform sampling attempts_per_point must be positive");
    }
    if (polygon_region && options.include_bounds_corners) {
        throw std::invalid_argument(
            "include_bounds_corners is only valid for bounds sampling");
    }
}

void validate_options(const BlueNoiseSamplingOptions& options) {
    if (options.candidates_per_point == 0 ||
        options.attempts_per_candidate == 0) {
        throw std::invalid_argument(
            "blue-noise candidate and attempt counts must be positive");
    }
}

void validate_options(const JitteredGridSamplingOptions& options) {
    if (!std::isfinite(options.jitter) ||
        options.jitter < 0.0 || options.jitter > 1.0) {
        throw std::invalid_argument(
            "jittered-grid jitter must be between zero and one");
    }
    if (options.attempts_per_point == 0) {
        throw std::invalid_argument(
            "jittered-grid attempts_per_point must be positive");
    }
}

std::vector<FloatPoint> convert_points(const std::vector<Point>& points) {
    std::vector<FloatPoint> converted;
    converted.reserve(points.size());
    for (const Point& point : points) {
        converted.push_back({
            static_cast<double>(point.x),
            static_cast<double>(point.y),
        });
    }
    return converted;
}

SamplingBounds polygon_bounds(
    const std::vector<FloatPoint>& points,
    const PolygonDomain& domain) {
    SamplingBounds bounds = empty_bounds();
    const auto inspect_ring = [&](const std::vector<std::uint32_t>& ring) {
        for (const std::uint32_t index : ring) {
            expand_bounds(bounds, points[index]);
        }
    };
    for_each_ring(domain, inspect_ring);
    return bounds;
}

std::vector<SamplingBounds> calculate_domain_bounds(
    const std::vector<FloatPoint>& points,
    const std::vector<PolygonDomain>& domains) {
    std::vector<SamplingBounds> bounds;
    bounds.reserve(domains.size());
    for (const PolygonDomain& domain : domains) {
        bounds.push_back(polygon_bounds(points, domain));
    }
    return bounds;
}

SamplingBounds combined_bounds(
    const std::vector<SamplingBounds>& domain_bounds) {
    SamplingBounds bounds = empty_bounds();
    for (const SamplingBounds& domain : domain_bounds) {
        expand_bounds(bounds, domain);
    }
    return bounds;
}

bool bounds_contain(
    const SamplingBounds& bounds,
    const FloatPoint& point) {
    return point.x >= bounds.min_x && point.x <= bounds.max_x &&
           point.y >= bounds.min_y && point.y <= bounds.max_y;
}

std::size_t find_containing_domain(
    const FloatPoint& candidate,
    const std::vector<FloatPoint>& points,
    const std::vector<PolygonDomain>& domains,
    const std::vector<SamplingBounds>& domain_bounds) {
    for (std::size_t domain = 0; domain < domains.size(); ++domain) {
        if (bounds_contain(domain_bounds[domain], candidate) &&
            detail::point_is_strictly_inside_domain_unchecked(
                candidate, domains[domain], points)) {
            return domain;
        }
    }
    return domains.size();
}

std::vector<FloatPoint> distinct_corners(const SamplingBounds& bounds) {
    const std::array<FloatPoint, 4> candidates = {{
        {bounds.min_x, bounds.min_y},
        {bounds.max_x, bounds.min_y},
        {bounds.max_x, bounds.max_y},
        {bounds.min_x, bounds.max_y},
    }};
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

struct RegionCandidate {
    FloatPoint point;
    std::size_t domain = 0;
};

RegionCandidate random_polygon_candidate(
    std::mt19937_64& random,
    const SamplingBounds& bounds,
    const std::vector<FloatPoint>& points,
    const std::vector<PolygonDomain>& domains,
    const std::vector<SamplingBounds>& domain_bounds,
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
        const std::size_t domain = find_containing_domain(
            candidate, points, domains, domain_bounds);
        if (domain != domains.size()) {
            return {candidate, domain};
        }
    }
    throw std::runtime_error(
        "could not sample a point inside the configured polygon domains");
}

SamplingScalar squared_distance(
    const FloatPoint& first,
    const FloatPoint& second) {
    const SamplingScalar dx =
        static_cast<SamplingScalar>(first.x) - second.x;
    const SamplingScalar dy =
        static_cast<SamplingScalar>(first.y) - second.y;
    return dx * dx + dy * dy;
}

SamplingScalar squared_distance_to_segment(
    const FloatPoint& point,
    const FloatPoint& a,
    const FloatPoint& b) {
    const SamplingScalar dx = static_cast<SamplingScalar>(b.x) - a.x;
    const SamplingScalar dy = static_cast<SamplingScalar>(b.y) - a.y;
    const SamplingScalar length_squared = dx * dx + dy * dy;
    if (length_squared == 0.0) {
        return squared_distance(point, a);
    }
    const SamplingScalar projection = std::clamp(
        ((static_cast<SamplingScalar>(point.x) - a.x) * dx +
         (static_cast<SamplingScalar>(point.y) - a.y) * dy) /
            length_squared,
        0.0,
        1.0);
    const SamplingScalar nearest_x = a.x + projection * dx;
    const SamplingScalar nearest_y = a.y + projection * dy;
    const SamplingScalar point_dx = point.x - nearest_x;
    const SamplingScalar point_dy = point.y - nearest_y;
    return point_dx * point_dx + point_dy * point_dy;
}

SamplingScalar polygon_boundary_distance_squared(
    const FloatPoint& point,
    const PolygonDomain& domain,
    const std::vector<FloatPoint>& points) {
    SamplingScalar minimum = std::numeric_limits<SamplingScalar>::max();
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
    for_each_ring(domain, inspect_ring);
    return minimum;
}

SamplingScalar bounds_boundary_distance_squared(
    const FloatPoint& point,
    const SamplingBounds& bounds) {
    const SamplingScalar distance = std::min({
        static_cast<SamplingScalar>(point.x - bounds.min_x),
        static_cast<SamplingScalar>(bounds.max_x - point.x),
        static_cast<SamplingScalar>(point.y - bounds.min_y),
        static_cast<SamplingScalar>(bounds.max_y - point.y),
    });
    return distance * distance;
}

class NearestNeighborGrid {
public:
    NearestNeighborGrid(
        const SamplingBounds& bounds,
        std::size_t point_count)
        : bounds_(bounds) {
        const double width = bounds.max_x - bounds.min_x;
        const double height = bounds.max_y - bounds.min_y;
        if (point_count == 0 || !std::isfinite(width) ||
            !std::isfinite(height) || width <= 0.0 || height <= 0.0) {
            return;
        }

        const SamplingScalar horizontal = std::ceil(std::sqrt(
            static_cast<SamplingScalar>(point_count) * width / height));
        const bool single_row =
            !std::isfinite(horizontal) ||
            horizontal >= static_cast<SamplingScalar>(point_count);
        x_cells_ = single_row
                       ? point_count
                       : std::max<std::size_t>(
                             1, static_cast<std::size_t>(horizontal));
        y_cells_ = std::max<std::size_t>(
            1, (point_count + x_cells_ - 1) / x_cells_);
        cell_width_ = width / static_cast<double>(x_cells_);
        cell_height_ = height / static_cast<double>(y_cells_);
        inverse_cell_width_ = 1.0 / cell_width_;
        inverse_cell_height_ = 1.0 / cell_height_;
        heads_.assign(x_cells_ * y_cells_, kNoPoint);
        next_.reserve(point_count);
    }

    void insert(const FloatPoint& point) {
        if (heads_.empty()) {
            return;
        }
        const auto [x, y] = cell(point);
        const std::size_t point_index = next_.size();
        const std::size_t index = flat_index(x, y);
        next_.push_back(heads_[index]);
        heads_[index] = point_index;
    }

    SamplingScalar nearest_distance_squared(
        const FloatPoint& candidate,
        const std::vector<FloatPoint>& accepted,
        SamplingScalar upper_bound,
        std::size_t ignored_point =
            std::numeric_limits<std::size_t>::max()) const {
        if (accepted.size() <= (ignored_point < accepted.size() ? 1U : 0U) ||
            upper_bound <= 0.0) {
            return upper_bound;
        }
        if (heads_.empty()) {
            for (std::size_t point = 0; point < accepted.size(); ++point) {
                if (point != ignored_point) {
                    upper_bound = std::min(
                        upper_bound,
                        squared_distance(candidate, accepted[point]));
                }
            }
            return upper_bound;
        }

        const auto [center_x, center_y] = cell(candidate);
        const std::size_t maximum_ring = std::max({
            center_x,
            x_cells_ - 1 - center_x,
            center_y,
            y_cells_ - 1 - center_y,
        });
        for (std::size_t ring = 0; ring <= maximum_ring; ++ring) {
            const CellRing cells = cell_ring(center_x, center_y, ring);
            inspect_ring(
                cells,
                candidate,
                accepted,
                upper_bound,
                ignored_point);
            if (upper_bound == 0.0 ||
                outside_distance_squared(cells, candidate) >= upper_bound) {
                break;
            }
        }
        return upper_bound;
    }

private:
    struct CellRing {
        std::size_t min_x = 0;
        std::size_t max_x = 0;
        std::size_t min_y = 0;
        std::size_t max_y = 0;
        bool has_bottom = false;
        bool has_top = false;
        bool has_left = false;
        bool has_right = false;
    };

    static constexpr std::size_t kNoPoint =
        std::numeric_limits<std::size_t>::max();

    std::pair<std::size_t, std::size_t> cell(
        const FloatPoint& point) const {
        const auto coordinate = [](double value, std::size_t count) {
            if (value <= 0.0) {
                return std::size_t{0};
            }
            return std::min(
                static_cast<std::size_t>(value), count - 1);
        };
        return {
            coordinate(
                (point.x - bounds_.min_x) * inverse_cell_width_,
                x_cells_),
            coordinate(
                (point.y - bounds_.min_y) * inverse_cell_height_,
                y_cells_),
        };
    }

    std::size_t flat_index(std::size_t x, std::size_t y) const {
        return y * x_cells_ + x;
    }

    CellRing cell_ring(
        std::size_t center_x,
        std::size_t center_y,
        std::size_t ring) const {
        return {
            ring <= center_x ? center_x - ring : 0,
            std::min(
                x_cells_ - 1,
                center_x + std::min(ring, x_cells_ - 1)),
            ring <= center_y ? center_y - ring : 0,
            std::min(
                y_cells_ - 1,
                center_y + std::min(ring, y_cells_ - 1)),
            ring <= center_y,
            ring <= y_cells_ - 1 - center_y,
            ring <= center_x,
            ring <= x_cells_ - 1 - center_x,
        };
    }

    void inspect_cell(
        std::size_t x,
        std::size_t y,
        const FloatPoint& candidate,
        const std::vector<FloatPoint>& accepted,
        SamplingScalar& minimum,
        std::size_t ignored_point) const {
        std::size_t point = heads_[flat_index(x, y)];
        while (point != kNoPoint) {
            if (point != ignored_point) {
                minimum = std::min(
                    minimum, squared_distance(candidate, accepted[point]));
            }
            point = next_[point];
        }
    }

    void inspect_ring(
        const CellRing& cells,
        const FloatPoint& candidate,
        const std::vector<FloatPoint>& accepted,
        SamplingScalar& minimum,
        std::size_t ignored_point) const {
        if (cells.has_bottom) {
            for (std::size_t x = cells.min_x; x <= cells.max_x; ++x) {
                inspect_cell(
                    x,
                    cells.min_y,
                    candidate,
                    accepted,
                    minimum,
                    ignored_point);
            }
        }
        if (cells.has_top &&
            (!cells.has_bottom || cells.max_y != cells.min_y)) {
            for (std::size_t x = cells.min_x; x <= cells.max_x; ++x) {
                inspect_cell(
                    x,
                    cells.max_y,
                    candidate,
                    accepted,
                    minimum,
                    ignored_point);
            }
        }

        const std::size_t first_y =
            cells.min_y + (cells.has_bottom ? 1 : 0);
        const std::size_t end_y =
            cells.max_y + 1 - (cells.has_top ? 1 : 0);
        if (first_y < end_y) {
            if (cells.has_left) {
                for (std::size_t y = first_y; y < end_y; ++y) {
                    inspect_cell(
                        cells.min_x,
                        y,
                        candidate,
                        accepted,
                        minimum,
                        ignored_point);
                }
            }
            if (cells.has_right &&
                (!cells.has_left || cells.max_x != cells.min_x)) {
                for (std::size_t y = first_y; y < end_y; ++y) {
                    inspect_cell(
                        cells.max_x,
                        y,
                        candidate,
                        accepted,
                        minimum,
                        ignored_point);
                }
            }
        }
    }

    SamplingScalar outside_distance_squared(
        const CellRing& cells,
        const FloatPoint& candidate) const {
        SamplingScalar distance = std::numeric_limits<SamplingScalar>::max();
        if (cells.min_x != 0) {
            distance = std::min(
                distance,
                static_cast<SamplingScalar>(candidate.x) -
                    (bounds_.min_x +
                     static_cast<double>(cells.min_x) * cell_width_));
        }
        if (cells.max_x + 1 != x_cells_) {
            distance = std::min(
                distance,
                static_cast<SamplingScalar>(
                    (bounds_.min_x +
                     static_cast<double>(cells.max_x + 1) * cell_width_) -
                    candidate.x));
        }
        if (cells.min_y != 0) {
            distance = std::min(
                distance,
                static_cast<SamplingScalar>(candidate.y) -
                    (bounds_.min_y +
                     static_cast<double>(cells.min_y) * cell_height_));
        }
        if (cells.max_y + 1 != y_cells_) {
            distance = std::min(
                distance,
                static_cast<SamplingScalar>(
                    (bounds_.min_y +
                     static_cast<double>(cells.max_y + 1) * cell_height_) -
                    candidate.y));
        }
        return distance * distance;
    }

    SamplingBounds bounds_;
    std::size_t x_cells_ = 0;
    std::size_t y_cells_ = 0;
    double cell_width_ = 0.0;
    double cell_height_ = 0.0;
    double inverse_cell_width_ = 0.0;
    double inverse_cell_height_ = 0.0;
    std::vector<std::size_t> heads_;
    std::vector<std::size_t> next_;
};

constexpr std::size_t kGapCandidates = 8;
constexpr std::size_t kMaximumDensityAttempts = 6;
constexpr std::size_t kAdjustmentDivisor = 100;
constexpr std::uint64_t kGapSeedMix = 0xd1b54a32d192ed03ULL;
constexpr std::array<double, 4> kRemovalSeparationFactors = {
    2.0,
    1.5,
    1.0,
    0.0,
};
constexpr SamplingScalar kDensityIncrease = 1.002;
constexpr SamplingScalar kDensityDecrease = 0.998;
constexpr SamplingScalar kPi =
    3.141592653589793238462643383279502884;
constexpr SamplingScalar kSqrtThree =
    1.732050807568877293527446341505872367;

SamplingScalar ring_area_twice(
    const std::vector<std::uint32_t>& ring,
    const std::vector<FloatPoint>& points) {
    SamplingScalar area = 0.0;
    for (std::size_t i = 0; i < ring.size(); ++i) {
        const FloatPoint& first = points[ring[i]];
        const FloatPoint& second = points[ring[(i + 1) % ring.size()]];
        area += static_cast<SamplingScalar>(first.x) * second.y -
                static_cast<SamplingScalar>(first.y) * second.x;
    }
    return std::abs(area);
}

SamplingScalar polygon_area(
    const std::vector<FloatPoint>& points,
    const std::vector<PolygonDomain>& domains) {
    SamplingScalar area_twice = 0.0;
    for (const PolygonDomain& domain : domains) {
        SamplingScalar domain_area = ring_area_twice(domain.outer_ring, points);
        for (const std::vector<std::uint32_t>& hole : domain.holes) {
            domain_area -= ring_area_twice(hole, points);
        }
        area_twice += std::max(0.0, domain_area);
    }
    return area_twice * 0.5;
}

class SamplingRegionView {
public:
    SamplingRegionView(
        const SamplingBounds& bounds,
        const std::vector<FloatPoint>& polygon_points,
        const std::vector<PolygonDomain>& domains,
        const std::vector<SamplingBounds>& domain_bounds)
        : bounds_(bounds),
          polygon_points_(polygon_points),
          domains_(domains),
          domain_bounds_(domain_bounds) {}

    const SamplingBounds& bounds() const { return bounds_; }
    bool is_polygon() const { return !domains_.empty(); }

    bool contains(const FloatPoint& point) const {
        return is_polygon()
                   ? find_containing_domain(
                         point, polygon_points_, domains_, domain_bounds_) !=
                         domains_.size()
                   : bounds_contain(bounds_, point);
    }

    SamplingScalar area() const {
        if (is_polygon()) {
            return polygon_area(polygon_points_, domains_);
        }
        return static_cast<SamplingScalar>(bounds_.max_x - bounds_.min_x) *
               (bounds_.max_y - bounds_.min_y);
    }

    RegionCandidate random_candidate(
        std::mt19937_64& random,
        std::size_t attempts) const {
        if (is_polygon()) {
            return random_polygon_candidate(
                random,
                bounds_,
                polygon_points_,
                domains_,
                domain_bounds_,
                attempts);
        }
        std::uniform_real_distribution<double> x_coordinate(
            bounds_.min_x, bounds_.max_x);
        std::uniform_real_distribution<double> y_coordinate(
            bounds_.min_y, bounds_.max_y);
        return {{x_coordinate(random), y_coordinate(random)}, 0};
    }

    SamplingScalar boundary_distance_squared(
        const RegionCandidate& candidate) const {
        return is_polygon()
                   ? polygon_boundary_distance_squared(
                         candidate.point,
                         domains_[candidate.domain],
                         polygon_points_)
                   : bounds_boundary_distance_squared(candidate.point, bounds_);
    }

private:
    const SamplingBounds& bounds_;
    const std::vector<FloatPoint>& polygon_points_;
    const std::vector<PolygonDomain>& domains_;
    const std::vector<SamplingBounds>& domain_bounds_;
};

FloatPoint choose_best_candidate(
    std::mt19937_64& random,
    std::size_t candidate_count,
    std::size_t attempts_per_candidate,
    const SamplingRegionView& region,
    const NearestNeighborGrid& grid,
    const std::vector<FloatPoint>& accepted) {
    FloatPoint best;
    SamplingScalar best_distance_squared = -1.0;
    for (std::size_t sample = 0; sample < candidate_count; ++sample) {
        const RegionCandidate candidate =
            region.random_candidate(random, attempts_per_candidate);
        const SamplingScalar minimum_distance_squared =
            grid.nearest_distance_squared(
                candidate.point,
                accepted,
                region.boundary_distance_squared(candidate));
        if (minimum_distance_squared > best_distance_squared) {
            best = candidate.point;
            best_distance_squared = minimum_distance_squared;
        }
    }
    return best;
}

std::vector<FloatPoint> trim_crowded_points(
    std::vector<FloatPoint> points,
    std::size_t point_count,
    const SamplingBounds& bounds,
    double nominal_spacing) {
    if (points.size() <= point_count) {
        return points;
    }

    struct Score {
        SamplingScalar nearest_distance_squared = 0.0;
        std::size_t point = 0;
    };

    NearestNeighborGrid grid(bounds, points.size());
    for (const FloatPoint& point : points) {
        grid.insert(point);
    }
    std::vector<Score> scores;
    scores.reserve(points.size());
    for (std::size_t point = 0; point < points.size(); ++point) {
        scores.push_back({
            grid.nearest_distance_squared(
                points[point],
                points,
                std::numeric_limits<SamplingScalar>::max(),
                point),
            point,
        });
    }
    std::sort(
        scores.begin(),
        scores.end(),
        [](const Score& first, const Score& second) {
            if (first.nearest_distance_squared !=
                second.nearest_distance_squared) {
                return first.nearest_distance_squared <
                       second.nearest_distance_squared;
            }
            return first.point < second.point;
        });

    const std::size_t removal_count = points.size() - point_count;
    std::vector<bool> removed(points.size(), false);
    std::vector<std::size_t> removals;
    removals.reserve(removal_count);
    for (const double factor : kRemovalSeparationFactors) {
        const SamplingScalar separation =
            static_cast<SamplingScalar>(factor) * nominal_spacing;
        const SamplingScalar separation_squared = separation * separation;
        for (const Score& score : scores) {
            if (removed[score.point]) {
                continue;
            }
            bool separated = true;
            for (const std::size_t earlier : removals) {
                if (squared_distance(
                        points[score.point], points[earlier]) <
                    separation_squared) {
                    separated = false;
                    break;
                }
            }
            if (!separated) {
                continue;
            }
            removed[score.point] = true;
            removals.push_back(score.point);
            if (removals.size() == removal_count) {
                break;
            }
        }
        if (removals.size() == removal_count) {
            break;
        }
    }

    std::vector<FloatPoint> retained;
    retained.reserve(point_count);
    for (std::size_t point = 0; point < points.size(); ++point) {
        if (!removed[point]) {
            retained.push_back(points[point]);
        }
    }
    return retained;
}

std::vector<FloatPoint> fill_sample_gaps(
    std::vector<FloatPoint> points,
    std::size_t point_count,
    const SamplingRegionView& region,
    std::size_t attempts_per_candidate,
    std::uint64_t seed) {
    if (points.size() >= point_count) {
        return points;
    }

    std::mt19937_64 random(seed ^ kGapSeedMix);
    NearestNeighborGrid grid(region.bounds(), point_count);
    for (const FloatPoint& point : points) {
        grid.insert(point);
    }

    while (points.size() < point_count) {
        const FloatPoint best = choose_best_candidate(
            random,
            kGapCandidates,
            attempts_per_candidate,
            region,
            grid,
            points);
        points.push_back(best);
        grid.insert(best);
    }
    return points;
}

std::size_t checked_product(
    std::size_t first,
    std::size_t second,
    const char* message) {
    if (first != 0 && second > std::numeric_limits<std::size_t>::max() / first) {
        throw std::runtime_error(message);
    }
    return first * second;
}

struct LatticeAttempt {
    std::vector<FloatPoint> points;
    double spacing = 0.0;
    std::size_t candidates_examined = 0;
};

[[noreturn]] void throw_candidate_limit_error() {
    throw std::runtime_error(
        "could not sample enough jittered-grid polygon points "
        "within the configured attempt limit");
}

LatticeAttempt generate_lattice_attempt(
    const SamplingRegionView& region,
    SamplingScalar region_area,
    std::size_t target_count,
    SamplingScalar target_density,
    double jitter,
    std::uint64_t seed,
    std::size_t remaining_candidates) {
    const SamplingScalar spacing_value = std::sqrt(
        2.0 * region_area / (kSqrtThree * target_density));
    if (!std::isfinite(spacing_value) || spacing_value <= 0.0 ||
        spacing_value > std::numeric_limits<double>::max()) {
        throw_candidate_limit_error();
    }

    LatticeAttempt attempt;
    attempt.spacing = static_cast<double>(spacing_value);
    const double row_spacing =
        static_cast<double>(0.5 * kSqrtThree) * attempt.spacing;
    const double jitter_radius = 0.5 * jitter * attempt.spacing;

    std::mt19937_64 random(seed);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    std::uniform_real_distribution<double> signed_unit(-1.0, 1.0);
    const double angle = unit(random) * static_cast<double>(kPi / 3.0);
    const double cosine = std::cos(angle);
    const double sine = std::sin(angle);
    const double phase_x = unit(random) * attempt.spacing;
    const double phase_y = unit(random) * row_spacing;

    const SamplingBounds& bounds = region.bounds();
    const double width = bounds.max_x - bounds.min_x;
    const double height = bounds.max_y - bounds.min_y;
    const double half_width = 0.5 * width;
    const double half_height = 0.5 * height;
    const double horizontal_extent =
        std::abs(cosine) * half_width + std::abs(sine) * half_height +
        jitter_radius;
    const double vertical_extent =
        std::abs(sine) * half_width + std::abs(cosine) * half_height +
        jitter_radius;
    const double first_y =
        std::floor((-vertical_extent - phase_y) / row_spacing) * row_spacing +
        phase_y;
    const SamplingScalar row_count_value = std::floor(
        (static_cast<SamplingScalar>(vertical_extent) - first_y) /
        row_spacing) + 1.0;
    if (!std::isfinite(row_count_value) || row_count_value <= 0.0 ||
        row_count_value > static_cast<SamplingScalar>(remaining_candidates)) {
        throw_candidate_limit_error();
    }
    const std::size_t row_count = static_cast<std::size_t>(row_count_value);
    const std::size_t reserve_limit =
        target_count > std::numeric_limits<std::size_t>::max() / 2
            ? std::numeric_limits<std::size_t>::max()
            : target_count * 2;
    attempt.points.reserve(std::min(remaining_candidates, reserve_limit));

    for (std::size_t row = 0; row < row_count; ++row) {
        const double lattice_y =
            first_y + static_cast<double>(row) * row_spacing;
        const double row_offset =
            (row & 1U) == 0U ? 0.0 : 0.5 * attempt.spacing;
        const double first_x =
            std::floor(
                (-horizontal_extent - phase_x - row_offset) /
                attempt.spacing) *
                attempt.spacing +
            phase_x + row_offset;
        const SamplingScalar column_count_value = std::floor(
            (static_cast<SamplingScalar>(horizontal_extent) - first_x) /
            attempt.spacing) + 1.0;
        if (!std::isfinite(column_count_value) ||
            column_count_value <= 0.0 ||
            column_count_value >
                static_cast<SamplingScalar>(remaining_candidates -
                                         attempt.candidates_examined)) {
            throw_candidate_limit_error();
        }
        const std::size_t column_count =
            static_cast<std::size_t>(column_count_value);
        for (std::size_t column = 0; column < column_count; ++column) {
            ++attempt.candidates_examined;
            const double lattice_x =
                first_x + static_cast<double>(column) * attempt.spacing;
            const FloatPoint lattice_point = {
                bounds.min_x + half_width + cosine * lattice_x -
                    sine * lattice_y,
                bounds.min_y + half_height + sine * lattice_x +
                    cosine * lattice_y,
            };
            if (lattice_point.x < bounds.min_x - jitter_radius ||
                lattice_point.x > bounds.max_x + jitter_radius ||
                lattice_point.y < bounds.min_y - jitter_radius ||
                lattice_point.y > bounds.max_y + jitter_radius) {
                continue;
            }

            double jitter_x = 0.0;
            double jitter_y = 0.0;
            if (jitter_radius != 0.0) {
                do {
                    jitter_x = signed_unit(random);
                    jitter_y = signed_unit(random);
                } while (jitter_x * jitter_x + jitter_y * jitter_y > 1.0);
                jitter_x *= jitter_radius;
                jitter_y *= jitter_radius;
            }
            const FloatPoint candidate = {
                lattice_point.x + cosine * jitter_x - sine * jitter_y,
                lattice_point.y + sine * jitter_x + cosine * jitter_y,
            };
            if (region.contains(candidate)) {
                attempt.points.push_back(candidate);
            }
        }
    }
    return attempt;
}

std::optional<std::vector<FloatPoint>> reconcile_sample_count(
    LatticeAttempt& attempt,
    std::size_t point_count,
    const SamplingRegionView& region,
    std::size_t attempts_per_candidate,
    std::uint64_t seed) {
    const std::size_t acceptable_adjustment =
        std::max<std::size_t>(1, point_count / kAdjustmentDivisor);
    if (attempt.points.size() == point_count) {
        return std::move(attempt.points);
    }
    if (attempt.points.size() < point_count &&
        point_count - attempt.points.size() <= acceptable_adjustment) {
        return fill_sample_gaps(
            std::move(attempt.points),
            point_count,
            region,
            attempts_per_candidate,
            seed);
    }
    if (attempt.points.size() > point_count &&
        attempt.points.size() - point_count <= acceptable_adjustment) {
        return trim_crowded_points(
            std::move(attempt.points),
            point_count,
            region.bounds(),
            attempt.spacing);
    }
    return std::nullopt;
}

SamplingScalar adjusted_density(
    SamplingScalar current_density,
    std::size_t observed_count,
    std::size_t target_count) {
    const SamplingScalar observed =
        std::max(1.0, static_cast<SamplingScalar>(observed_count));
    const SamplingScalar correction =
        static_cast<SamplingScalar>(target_count) / observed;
    if (observed_count < target_count) {
        return std::max(
            current_density * kDensityIncrease,
            current_density * correction * kDensityIncrease);
    }
    return std::min(
        current_density * kDensityDecrease,
        current_density * correction * kDensityDecrease);
}

std::vector<FloatPoint> generate_jittered_samples(
    const SamplingRegionView& region,
    SamplingScalar region_area,
    const JitteredGridSamplingOptions& options,
    std::size_t maximum_candidates) {
    SamplingScalar target_density =
        static_cast<SamplingScalar>(options.point_count);
    std::size_t candidates_examined = 0;
    std::optional<LatticeAttempt> best_surplus;

    for (std::size_t density_attempt = 0;; ++density_attempt) {
        LatticeAttempt attempt = generate_lattice_attempt(
            region,
            region_area,
            options.point_count,
            target_density,
            options.jitter,
            options.seed,
            maximum_candidates - candidates_examined);
        candidates_examined += attempt.candidates_examined;

        if (auto reconciled = reconcile_sample_count(
                attempt,
                options.point_count,
                region,
                options.attempts_per_point,
                options.seed)) {
            return std::move(*reconciled);
        }

        const std::size_t generated_count = attempt.points.size();
        if (generated_count > options.point_count) {
            const std::size_t surplus =
                generated_count - options.point_count;
            if (!best_surplus ||
                surplus <
                    best_surplus->points.size() - options.point_count) {
                best_surplus = std::move(attempt);
            }
        }
        if (density_attempt + 1 >= kMaximumDensityAttempts &&
            best_surplus) {
            return trim_crowded_points(
                std::move(best_surplus->points),
                options.point_count,
                region.bounds(),
                best_surplus->spacing);
        }

        target_density = adjusted_density(
            target_density, generated_count, options.point_count);
        if (!std::isfinite(target_density) || target_density <= 0.0 ||
            target_density > static_cast<SamplingScalar>(maximum_candidates)) {
            throw_candidate_limit_error();
        }
    }
}

}  // namespace

void PointSampler::require_configured_region() const {
    if (region_ == RegionKind::None) {
        throw std::logic_error(
            "point sampler requires configured bounds or polygon interiors");
    }
}

void PointSampler::set_bounds(SamplingBounds bounds) {
    validate_bounds(bounds);
    bounds_ = bounds;
    polygon_points_.clear();
    domains_.clear();
    domain_bounds_.clear();
    region_ = RegionKind::Bounds;
}

void PointSampler::set_polygon_interiors(
    std::vector<Point> points,
    std::vector<PolygonDomain> domains) {
    set_polygon_interiors(convert_points(points), std::move(domains));
}

void PointSampler::set_polygon_interiors(
    std::vector<FloatPoint> points,
    std::vector<PolygonDomain> domains) {
    validate_polygon_region(points, domains);
    std::vector<SamplingBounds> domain_bounds =
        calculate_domain_bounds(points, domains);
    const SamplingBounds bounds = combined_bounds(domain_bounds);

    bounds_ = bounds;
    polygon_points_ = std::move(points);
    domains_ = std::move(domains);
    domain_bounds_ = std::move(domain_bounds);
    region_ = RegionKind::PolygonInteriors;
}

std::vector<FloatPoint> PointSampler::generate_uniform(
    const UniformSamplingOptions& options) const {
    require_configured_region();
    validate_options(
        options, region_ == RegionKind::PolygonInteriors);

    std::vector<FloatPoint> generated;
    if (region_ == RegionKind::Bounds && options.include_bounds_corners) {
        generated = distinct_corners(bounds_);
        if (options.point_count < generated.size()) {
            throw std::invalid_argument(
                "point_count is smaller than the number of distinct bounds "
                "corners");
        }
    }
    generated.reserve(options.point_count);

    std::mt19937_64 random(options.seed);
    if (region_ == RegionKind::Bounds) {
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
            random_polygon_candidate(
                random,
                bounds_,
                polygon_points_,
                domains_,
                domain_bounds_,
                options.attempts_per_point)
                .point);
    }
    return generated;
}

std::vector<FloatPoint> PointSampler::generate_blue_noise(
    const BlueNoiseSamplingOptions& options) const {
    require_configured_region();
    validate_options(options);

    const SamplingRegionView region(
        bounds_, polygon_points_, domains_, domain_bounds_);
    std::mt19937_64 random(options.seed);
    std::vector<FloatPoint> generated;
    generated.reserve(options.point_count);
    NearestNeighborGrid sample_grid(bounds_, options.point_count);
    while (generated.size() < options.point_count) {
        const FloatPoint best = choose_best_candidate(
            random,
            options.candidates_per_point,
            options.attempts_per_candidate,
            region,
            sample_grid,
            generated);
        generated.push_back(best);
        sample_grid.insert(best);
    }
    return generated;
}

std::vector<FloatPoint> PointSampler::generate_jittered_grid(
    const JitteredGridSamplingOptions& options) const {
    require_configured_region();
    validate_options(options);
    if (options.point_count == 0) {
        return {};
    }

    const double width = bounds_.max_x - bounds_.min_x;
    const double height = bounds_.max_y - bounds_.min_y;
    if (!std::isfinite(width) || !std::isfinite(height) ||
        width <= 0.0 || height <= 0.0) {
        throw std::runtime_error(
            "could not sample a jittered grid in a zero-area region");
    }

    const std::size_t maximum_candidates = checked_product(
        options.point_count,
        options.attempts_per_point,
        "jittered-grid candidate limit overflows size_t");
    const SamplingRegionView region(
        bounds_, polygon_points_, domains_, domain_bounds_);
    const SamplingScalar region_area = region.area();
    if (!std::isfinite(region_area) || region_area <= 0.0) {
        throw std::runtime_error(
            "could not sample a point inside the configured polygon domains");
    }

    return generate_jittered_samples(
        region, region_area, options, maximum_candidates);
}

}  // namespace delaunay32::extras
