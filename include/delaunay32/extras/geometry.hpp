// SPDX-License-Identifier: MIT

#pragma once

#include "delaunay32/delaunay.hpp"
#include "delaunay32/quantization.hpp"

#include <cstdint>
#include <optional>
#include <vector>

namespace delaunay32::extras {

// Data model used by the Delaunay32 geometry JSON helpers. `polygon` describes
// one domain; `polygons` describes independent domains over the same points.
struct Geometry {
    std::vector<Point> points;
    std::vector<Constraint> constraints;
    std::optional<PolygonDomain> polygon;
    std::vector<PolygonDomain> polygons;
};

// Returns false for points on an outer or hole boundary. Ring indices must be
// valid for `points`; rings with fewer than three indices throw
// std::invalid_argument.
bool point_is_strictly_inside_domain(
    const Point& point,
    const PolygonDomain& domain,
    const std::vector<Point>& points);

// Floating-point equivalent of the exact integer query above. Coordinates
// must be finite.
bool point_is_strictly_inside_domain(
    const FloatPoint& point,
    const PolygonDomain& domain,
    const std::vector<FloatPoint>& points);

}  // namespace delaunay32::extras
