// SPDX-License-Identifier: MIT

#pragma once

#include "delaunay32/delaunay.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace delaunay32_example {

struct PolygonDomain {
    std::vector<std::uint32_t> outer_ring;
    std::vector<std::vector<std::uint32_t>> holes;
};

// The example JSON schema keeps coordinates separate from topology so the
// same fixture can describe ordinary points, constrained edges, and polygon
// rings without changing point indices. `polygon` describes one domain;
// `polygons` describes multiple independent domains over the same points.
struct GeometryInput {
    std::vector<delaunay32::Point> points;
    std::vector<delaunay32::Constraint> constraints;
    std::vector<std::uint32_t> outer_ring;
    std::vector<std::vector<std::uint32_t>> holes;
    std::vector<PolygonDomain> polygons;
};

struct Options {
    bool input_mode = false;
    std::string input_path;
    std::string output_path = "mesh.svg";
    std::size_t point_count = 1000;
    std::uint64_t seed = 1;
};

Options parse_options(int argc, char** argv);
void print_usage(const char* executable);
GeometryInput load_geometry(const Options& options);
GeometryInput read_geometry_json(const std::string& input_path);
bool point_is_strictly_inside_domain(
    const delaunay32::Point& point,
    const PolygonDomain& domain,
    const std::vector<delaunay32::Point>& points);

}  // namespace delaunay32_example
