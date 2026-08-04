// SPDX-License-Identifier: MIT

#pragma once

#include "delaunay32/extras/geometry.hpp"

#include <cstddef>
#include <cstdint>
#include <string>

namespace delaunay32_example {

struct Options {
    bool input_mode = false;
    std::string input_path;
    std::string output_path = "mesh.svg";
    std::size_t point_count = 1000;
    std::uint64_t seed = 1;
};

Options parse_options(int argc, char** argv);
void print_usage(const char* executable);
delaunay32::extras::Geometry load_geometry(const Options& options);

}  // namespace delaunay32_example
