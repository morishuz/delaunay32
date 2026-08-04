// SPDX-License-Identifier: MIT

#pragma once

#include "delaunay32/delaunay.hpp"
#include "delaunay32/extras/geometry.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace delaunay32::extras {

struct IntBounds {
    std::int32_t min_x = 0;
    std::int32_t max_x = 999;
    std::int32_t min_y = 0;
    std::int32_t max_y = 999;
};

struct FloatBounds {
    float min_x = 0.0F;
    float max_x = 1.0F;
    float min_y = 0.0F;
    float max_y = 1.0F;
};

struct UniformIntOptions {
    std::size_t point_count = 1000;
    IntBounds bounds;
    std::uint64_t seed = 1;
    bool include_corners = true;
};

struct UniformFloatOptions {
    std::size_t point_count = 1000;
    FloatBounds bounds;
    std::uint64_t seed = 1;
    bool include_corners = true;
};

// Generates unique integer points. point_count includes any requested corners.
std::vector<Point> generate_uniform_int_points(
    const UniformIntOptions& options = {});

// Generates floating-point samples. point_count includes requested corners.
std::vector<FloatPoint> generate_uniform_float_points(
    const UniformFloatOptions& options = {});

struct BestCandidateOptions {
    std::size_t point_count = 1000;
    std::size_t candidates_per_point = 16;
    std::size_t attempts_per_candidate = 10000;
    std::uint64_t seed = 1;
};

// Returns new unique points strictly inside the supplied polygon domains.
// Candidates maximize their distance from domain boundaries and earlier
// samples in the same domain, producing blue-noise-style spacing. The returned
// vector contains only generated samples; callers can append it to `points`.
std::vector<Point> sample_polygon_interiors(
    const std::vector<Point>& points,
    const std::vector<PolygonDomain>& domains,
    const BestCandidateOptions& options = {});

}  // namespace delaunay32::extras
