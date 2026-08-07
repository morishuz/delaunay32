// SPDX-License-Identifier: MIT

#pragma once

#include "delaunay32/extras/geometry.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace delaunay32::extras::detail {

// Sampling inputs are already binary64. Keeping their geometry in binary64
// also avoids software-emulated IEEE-128 arithmetic on WebAssembly targets.
using SamplingScalar = double;

void validate_domain(
    const PolygonDomain& domain,
    std::size_t point_count,
    const std::string& label);

bool point_is_strictly_inside_domain_unchecked(
    const Point& point,
    const PolygonDomain& domain,
    const std::vector<Point>& points);

bool point_is_strictly_inside_domain_unchecked(
    const FloatPoint& point,
    const PolygonDomain& domain,
    const std::vector<FloatPoint>& points);

}  // namespace delaunay32::extras::detail
