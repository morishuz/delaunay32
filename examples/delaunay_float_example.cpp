// SPDX-License-Identifier: MIT

#include "delaunay32/delaunay.hpp"
#include "svg_io.hpp"

#include <cstddef>
#include <cstdint>
#include <exception>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr std::size_t kPointCount = 5000;
constexpr std::uint64_t kSeed = 42;
constexpr float kMinX = -123.75F;
constexpr float kMaxX = 876.25F;
constexpr float kMinY = 40.125F;
constexpr float kMaxY = 640.625F;

std::vector<delaunay32::FloatPoint> make_random_points() {
    // Include the domain corners so the SVG bounds are stable, then fill the
    // interior with deterministic random float coordinates.
    std::vector<delaunay32::FloatPoint> points = {
        {kMinX, kMinY},
        {kMaxX, kMinY},
        {kMaxX, kMaxY},
        {kMinX, kMaxY},
    };
    points.reserve(kPointCount);

    std::mt19937_64 random(kSeed);
    std::uniform_real_distribution<float> x_coordinate(kMinX, kMaxX);
    std::uniform_real_distribution<float> y_coordinate(kMinY, kMaxY);
    while (points.size() < kPointCount) {
        points.push_back({
            x_coordinate(random),
            y_coordinate(random),
        });
    }
    return points;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc > 2) {
            throw std::invalid_argument(
                "expected at most one argument: output.svg");
        }
        const std::string output_path =
            argc == 2 ? argv[1] : "delaunay32_float.svg";

        // These are the caller-owned source coordinates. triangulate_float()
        // reads them but never changes or replaces them.
        const std::vector<delaunay32::FloatPoint> points =
            make_random_points();

        delaunay32::Triangulator triangulator(0);
        const delaunay32::TriangulationResult result =
            triangulator.triangulate_float_full(points);
        const std::vector<delaunay32::Triangle>& triangles =
            result.triangles;
        const delaunay32::QuantizationReport& report =
            result.quantization;

        // Triangle stores only indices, so rendering uses the original floats.
        delaunay32_example::write_svg(
            output_path, points, triangles, report);

        std::cout << std::setprecision(9)
                  << "wrote " << output_path << '\n'
                  << "source points: " << points.size() << '\n'
                  << "triangles: " << triangles.size() << '\n'
                  << "quantization origin: (" << report.origin_x << ", "
                  << report.origin_y << ")\n"
                  << "integer units per source unit: " << report.scale
                  << '\n'
                  << "source units per grid step: " << report.grid_step
                  << '\n'
                  << "maximum coordinate error: "
                  << report.max_coordinate_error << '\n'
                  << "unique grid points: " << report.unique_points << '\n'
                  << "collapsed points: " << report.collapsed_points << '\n';

        if (!triangles.empty()) {
            const delaunay32::Triangle& first = triangles.front();
            std::cout << "first triangle indexes original points: "
                      << first.i0 << ", " << first.i1 << ", " << first.i2
                      << '\n';
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Float example error: " << error.what() << '\n'
                  << "Usage: " << argv[0] << " [output.svg]\n";
        return 1;
    }
}
