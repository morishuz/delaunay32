// SPDX-License-Identifier: MIT

#pragma once

#include "delaunay32/delaunay.hpp"

#include <cstddef>
#include <vector>

namespace delaunay32 {

// A floating-point source site for explicit conversion to Point. The returned
// integer vector preserves the source vector's length and index order.
struct FloatPoint {
    double x = 0.0;
    double y = 0.0;
};

enum class QuantizationMode {
    // Use the finest uniform grid supported by the predicate backend.
    Automatic,
    // Use QuantizationOptions::grid_step and the input minima as the origin.
    GridStep,
    // Use QuantizationOptions::origin_x, origin_y, and scale verbatim.
    FixedScale,
};

enum class QuantizationCollisionPolicy {
    Allow,
    Reject,
};

struct QuantizationOptions {
    QuantizationMode mode = QuantizationMode::Automatic;
    double grid_step = 0.0;
    double origin_x = 0.0;
    double origin_y = 0.0;
    double scale = 0.0;
    double max_coordinate_error = 0.0;
    QuantizationCollisionPolicy collision_policy =
        QuantizationCollisionPolicy::Allow;
};

// Quantized coordinates are round((value - origin) * scale).
struct QuantizationReport {
    double origin_x = 0.0;
    double origin_y = 0.0;
    double scale = 0.0;
    double grid_step = 0.0;
    double max_coordinate_error = 0.0;
    std::size_t unique_points = 0;
    std::size_t collapsed_points = 0;
};

struct QuantizationResult {
    std::vector<Point> points;
    QuantizationReport report;
};

QuantizationResult quantize(
    const std::vector<FloatPoint>& points,
    const QuantizationOptions& options = {});

}  // namespace delaunay32
