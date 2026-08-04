// SPDX-License-Identifier: MIT

#pragma once

#include "delaunay32/extras/geometry.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace delaunay32::extras::detail {

void validate_domain(
    const PolygonDomain& domain,
    std::size_t point_count,
    const std::string& label);

bool point_is_strictly_inside_domain_unchecked(
    const Point& point,
    const PolygonDomain& domain,
    const std::vector<Point>& points);

}  // namespace delaunay32::extras::detail
