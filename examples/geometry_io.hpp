// SPDX-License-Identifier: MIT

#pragma once

#include "delaunay32/delaunay.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace delaunay32_example {

// The example JSON schema keeps coordinates separate from topology so the
// same fixture can describe ordinary points, constrained edges, and polygon
// rings without changing point indices.
struct GeometryInput {
    std::vector<delaunay32::Point> points;
    std::vector<delaunay32::Constraint> constraints;
    std::vector<std::uint32_t> outer_ring;
    std::vector<std::vector<std::uint32_t>> holes;
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

}  // namespace delaunay32_example
