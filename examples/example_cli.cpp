// SPDX-License-Identifier: MIT

#include "example_cli.hpp"

#include "delaunay32/extras/json.hpp"
#include "delaunay32/extras/sampling.hpp"

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace delaunay32_example {
namespace {

std::size_t parse_point_count(const char* value) {
    std::size_t consumed = 0;
    const std::string text = value;
    const unsigned long long parsed = std::stoull(text, &consumed);
    if (consumed != text.size()) {
        throw std::invalid_argument("point count is not an integer");
    }
    if (parsed > std::numeric_limits<std::size_t>::max()) {
        throw std::invalid_argument("point count is too large");
    }
    return static_cast<std::size_t>(parsed);
}

std::uint64_t parse_seed(const char* value) {
    std::size_t consumed = 0;
    const std::string text = value;
    const std::uint64_t seed = std::stoull(text, &consumed);
    if (consumed != text.size()) {
        throw std::invalid_argument("seed is not an integer");
    }
    return seed;
}

}  // namespace

Options parse_options(int argc, char** argv) {
    Options options;
    if (argc >= 2 && std::string(argv[1]) == "--input") {
        options.input_mode = true;
        for (int i = 1; i < argc; ++i) {
            const std::string argument = argv[i];
            if (argument == "--input" && i + 1 < argc) {
                options.input_path = argv[++i];
            } else if (argument == "--output" && i + 1 < argc) {
                options.output_path = argv[++i];
            } else {
                throw std::invalid_argument(
                    "unknown or incomplete option: " + argument);
            }
        }
        if (options.input_path.empty()) {
            throw std::invalid_argument("--input requires a JSON path");
        }
        return options;
    }

    if (argc > 4) {
        throw std::invalid_argument("too many positional arguments");
    }
    options.point_count = argc >= 2 ? parse_point_count(argv[1]) : 1000;
    options.output_path = argc >= 3 ? argv[2] : "mesh.svg";
    options.seed = argc >= 4 ? parse_seed(argv[3]) : 1;
    return options;
}

void print_usage(const char* executable) {
    std::cerr
        << "Usage:\n"
        << "  " << executable
        << " [point-count] [output.svg] [seed]\n"
        << "  " << executable
        << " --input geometry.json [--output output.svg]\n";
}

delaunay32::extras::Geometry load_geometry(const Options& options) {
    if (options.input_mode) {
        return delaunay32::extras::read_geometry_json(options.input_path);
    }
    delaunay32::extras::UniformIntOptions sampling;
    sampling.point_count = options.point_count;
    sampling.bounds = {0, 999, 0, 999};
    sampling.seed = options.seed;

    delaunay32::extras::Geometry geometry;
    geometry.points =
        delaunay32::extras::generate_uniform_int_points(sampling);
    return geometry;
}

}  // namespace delaunay32_example
