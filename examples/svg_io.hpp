// SPDX-License-Identifier: MIT

#pragma once

#include "delaunay32/delaunay.hpp"
#include "geometry_io.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace delaunay32_example {

void write_svg(
    const std::string& output_path,
    const std::vector<delaunay32::Point>& points,
    const std::vector<delaunay32::Triangle>& triangles);
void write_svg(
    const std::string& output_path,
    const std::vector<delaunay32::FloatPoint>& points,
    const std::vector<delaunay32::Triangle>& triangles,
    const delaunay32::QuantizationReport& report);
void write_constrained_comparison_svg(
    const std::string& output_path,
    const std::vector<delaunay32::Point>& points,
    const std::vector<delaunay32::Constraint>& constraints,
    const std::vector<delaunay32::Triangle>& ordinary_triangles,
    const std::vector<delaunay32::Triangle>& constrained_triangles);
void write_polygon_svg(
    const std::string& output_path,
    const std::vector<delaunay32::Point>& points,
    const std::vector<std::uint32_t>& outer_ring,
    const std::vector<std::vector<std::uint32_t>>& holes,
    const std::vector<delaunay32::Triangle>& triangles);
void write_logo_polygon_svg(
    const std::string& output_path,
    const std::vector<delaunay32::Point>& points,
    std::size_t interior_point_count,
    const std::vector<PolygonDomain>& domains,
    const std::vector<delaunay32::Triangle>& triangles);

}  // namespace delaunay32_example
