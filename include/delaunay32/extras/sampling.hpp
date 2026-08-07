// SPDX-License-Identifier: MIT

#pragma once

#include "delaunay32/delaunay.hpp"
#include "delaunay32/quantization.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace delaunay32::extras {

// Axis-aligned sampling rectangle, in {minimum x, maximum x, minimum y,
// maximum y} order.
struct SamplingBounds {
    double min_x = 0.0;
    double max_x = 1.0;
    double min_y = 0.0;
    double max_y = 1.0;
};

struct UniformSamplingOptions {
    std::size_t point_count = 1000;
    std::uint64_t seed = 1;
    // Bounds sampling can reserve the first sites for distinct corners.
    bool include_bounds_corners = false;
    // Maximum rejection attempts for each polygon-interior point.
    std::size_t attempts_per_point = 10000;
};

struct BlueNoiseSamplingOptions {
    std::size_t point_count = 1000;
    // More candidates improve spacing at a proportional runtime cost.
    std::size_t candidates_per_point = 16;
    // Maximum rejection attempts for each polygon-interior candidate.
    std::size_t attempts_per_candidate = 10000;
    std::uint64_t seed = 1;
};

struct JitteredGridSamplingOptions {
    std::size_t point_count = 1000;
    // Zero preserves triangular-lattice sites. One permits a displacement of
    // half the lattice spacing; the displacement direction is isotropic.
    double jitter = 0.75;
    // Bounds the work needed to cover sparse polygon regions.
    std::size_t attempts_per_point = 10000;
    std::uint64_t seed = 1;
};

// Stateful floating-point sampler for either rectangular bounds or indexed
// polygon interiors. Setting one region replaces the previous region. Each
// generation call owns its random engine, so equal regions and options are
// deterministic.
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
    std::vector<FloatPoint> generate_jittered_grid(
        const JitteredGridSamplingOptions& options = {}) const;

private:
    enum class RegionKind {
        None,
        Bounds,
        PolygonInteriors,
    };

    void require_configured_region() const;

    RegionKind region_ = RegionKind::None;
    SamplingBounds bounds_;
    std::vector<FloatPoint> polygon_points_;
    std::vector<PolygonDomain> domains_;
    std::vector<SamplingBounds> domain_bounds_;
};

}  // namespace delaunay32::extras
