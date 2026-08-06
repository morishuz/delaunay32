// SPDX-License-Identifier: MIT

#pragma once

#include "delaunay32/delaunay.hpp"
#include "delaunay32/quantization.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace delaunay32::extras {

struct SamplingBounds {
    double min_x = 0.0;
    double max_x = 1.0;
    double min_y = 0.0;
    double max_y = 1.0;
};

struct UniformSamplingOptions {
    std::size_t point_count = 1000;
    std::uint64_t seed = 1;
    bool include_bounds_corners = false;
    std::size_t attempts_per_point = 10000;
};

struct BlueNoiseSamplingOptions {
    std::size_t point_count = 1000;
    std::size_t candidates_per_point = 16;
    std::size_t attempts_per_candidate = 10000;
    std::uint64_t seed = 1;
};

// Stateful floating-point sampler for either rectangular bounds or indexed
// polygon interiors. Setting one region replaces the previous region.
class PointSampler {
public:
    void set_bounds(SamplingBounds bounds);

    // The sampler owns its configured polygon data. Integer coordinates are
    // exactly representable by FloatPoint and are converted internally.
    void set_polygon_interiors(
        std::vector<Point> points,
        std::vector<PolygonDomain> domains);
    void set_polygon_interiors(
        std::vector<FloatPoint> points,
        std::vector<PolygonDomain> domains);

    std::vector<FloatPoint> generate_uniform(
        const UniformSamplingOptions& options = {}) const;
    std::vector<FloatPoint> generate_blue_noise(
        const BlueNoiseSamplingOptions& options = {}) const;

private:
    enum class Region {
        None,
        Bounds,
        PolygonInteriors,
    };

    Region region_ = Region::None;
    SamplingBounds bounds_;
    std::vector<FloatPoint> polygon_points_;
    std::vector<PolygonDomain> domains_;
};

}  // namespace delaunay32::extras
