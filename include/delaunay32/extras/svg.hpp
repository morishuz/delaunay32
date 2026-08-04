// SPDX-License-Identifier: MIT

#pragma once

#include "delaunay32/delaunay.hpp"
#include "delaunay32/extras/geometry.hpp"

#include <string>
#include <vector>

namespace delaunay32::extras {

void write_mesh_svg(
    const std::string& output_path,
    const std::vector<Point>& points,
    const std::vector<Triangle>& triangles);

void write_mesh_svg(
    const std::string& output_path,
    const std::vector<FloatPoint>& points,
    const std::vector<Triangle>& triangles,
    const QuantizationReport& report);

void write_polygon_mesh_svg(
    const std::string& output_path,
    const std::vector<Point>& points,
    const PolygonDomain& domain,
    const std::vector<Triangle>& triangles);

}  // namespace delaunay32::extras
