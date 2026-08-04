// SPDX-License-Identifier: MIT

#pragma once

#include "delaunay32/delaunay.hpp"
#include "delaunay32/extras/geometry.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace delaunay32_example {

// Presentation-specific renderers used only by the comparison and logo demos.
void write_constrained_comparison_svg(
    const std::string& output_path,
    const std::vector<delaunay32::Point>& points,
    const std::vector<delaunay32::Constraint>& constraints,
    const std::vector<delaunay32::Triangle>& ordinary_triangles,
    const std::vector<delaunay32::Triangle>& constrained_triangles);

void write_logo_polygon_svg(
    const std::string& output_path,
    const std::vector<delaunay32::Point>& points,
    std::size_t interior_point_count,
    const std::vector<delaunay32::extras::PolygonDomain>& domains,
    const std::vector<delaunay32::Triangle>& triangles);

}  // namespace delaunay32_example
