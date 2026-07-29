// SPDX-License-Identifier: MIT

#pragma once

#include "delaunay32/delaunay.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace delaunay32_example {

struct Options {
    bool csv_mode = false;
    std::string input_path;
    std::string output_path = "mesh.svg";
    std::size_t point_count = 1000;
    std::uint64_t seed = 1;
};

Options parse_options(int argc, char** argv);
void print_usage(const char* executable);
std::vector<delaunay32::Point> load_points(const Options& options);
void write_svg(
    const std::string& output_path,
    const std::vector<delaunay32::Point>& points,
    const std::vector<delaunay32::Triangle>& triangles);

}  // namespace delaunay32_example
