// SPDX-License-Identifier: MIT

#include "delaunay32/delaunay.hpp"
#include "delaunay32/extras/svg.hpp"
#include "delaunay32/quantization.hpp"

#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void print_result(
    const std::string& name,
    const std::vector<delaunay32::FloatPoint>& source,
    const delaunay32::QuantizationResult& result) {
    const delaunay32::QuantizationReport& report = result.report;
    std::cout << "\n" << name << '\n'
              << "  origin: (" << std::fixed << std::setprecision(6)
              << report.origin_x << ", " << report.origin_y << ")\n"
              << "  scale: " << std::scientific << report.scale << '\n'
              << "  grid step: " << report.grid_step << '\n'
              << "  maximum error: " << report.max_coordinate_error << '\n'
              << "  unique points: " << report.unique_points << '\n'
              << "  collapsed points: " << report.collapsed_points << '\n';

    for (std::size_t i = 0; i < source.size(); ++i) {
        std::cout << "  [" << i << "] (" << std::fixed
                  << std::setprecision(6) << source[i].x << ", "
                  << source[i].y << ") -> (" << result.points[i].x
                  << ", " << result.points[i].y << ")\n";
    }
}

void show_rejection(
    const std::string& name,
    const std::vector<delaunay32::FloatPoint>& points,
    const delaunay32::QuantizationOptions& options) {
    try {
        static_cast<void>(delaunay32::quantize(points, options));
        std::cout << "\n" << name << ": unexpectedly accepted\n";
    } catch (const std::invalid_argument& error) {
        std::cout << "\n" << name << ": " << error.what() << '\n';
    }
}

}  // namespace

int main() {
    // Delaunay32 triangulates signed 32-bit integer coordinates, but many
    // applications begin with floating-point world coordinates. These points
    // are deliberately near 10^12 while being only thousandths apart. Their
    // absolute values are far outside int32, but their small local shape still
    // needs to be preserved.
    //
    // quantize() solves this by moving the points near an origin and scaling
    // their small differences onto a safe integer grid:
    //
    //     integer = round((source - origin) * scale)
    //
    // It returns new integer points for triangulation without modifying the
    // original floating-point coordinates used for rendering or other work.
    constexpr double base_x = 1'000'000'000'000.0;
    constexpr double base_y = -2'000'000'000'000.0;
    const std::vector<delaunay32::FloatPoint> source_points = {
        {base_x + 0.000, base_y + 0.000},
        {base_x + 0.012, base_y + 0.001},
        {base_x + 0.015, base_y + 0.010},
        {base_x + 0.010, base_y + 0.017},
        {base_x + 0.001, base_y + 0.015},
        {base_x - 0.004, base_y + 0.008},
        {base_x + 0.001, base_y + 0.001},
        {base_x + 0.009, base_y + 0.004},
        {base_x + 0.012, base_y + 0.009},
        {base_x + 0.008, base_y + 0.013},
        {base_x + 0.003, base_y + 0.011},
        {base_x + 0.007, base_y + 0.008},
    };

    // Automatic mode chooses an origin near the data and the finest safe
    // integer grid that fits this complete batch.
    const delaunay32::QuantizationResult automatic =
        delaunay32::quantize(source_points);
    print_result("Automatic", source_points, automatic);

    // GridStep lets us choose the spacing between integer grid lines in the
    // original coordinate units. This deliberately coarse grid snaps two
    // nearby source points to the same integer coordinate.
    delaunay32::QuantizationOptions grid_options;
    grid_options.mode = delaunay32::QuantizationMode::GridStep;
    grid_options.grid_step = 0.005;
    const delaunay32::QuantizationResult grid =
        delaunay32::quantize(source_points, grid_options);
    print_result("Grid step (collisions allowed)", source_points, grid);

    // FixedScale specifies the complete mapping. Reusing the same origin and
    // scale gives separate batches the same integer coordinate system.
    delaunay32::QuantizationOptions fixed_options;
    fixed_options.mode = delaunay32::QuantizationMode::FixedScale;
    fixed_options.origin_x = base_x;
    fixed_options.origin_y = base_y;
    fixed_options.scale = 1000.0;
    const delaunay32::QuantizationResult fixed =
        delaunay32::quantize(source_points, fixed_options);
    print_result("Fixed scale", source_points, fixed);

    // By default collisions are reported but allowed. Applications that need
    // every source point to remain distinct can reject them instead.
    delaunay32::QuantizationOptions reject_collisions = grid_options;
    reject_collisions.collision_policy =
        delaunay32::QuantizationCollisionPolicy::Reject;
    show_rejection(
        "Collision policy", source_points, reject_collisions);

    // Quantization rounds coordinates to grid lines. An error limit rejects a
    // mapping when that movement is too large in the original coordinate units.
    delaunay32::QuantizationOptions strict_error = grid_options;
    strict_error.max_coordinate_error = 0.0001;
    show_rejection("Error limit", source_points, strict_error);

    // The triangulator uses the converted integers, but its triangle indices
    // still address source_points. We can therefore render the original values.
    delaunay32::Triangulator triangulator;
    triangulator.set_points(automatic.points);
    const delaunay32::TriangulationResult mesh =
        triangulator.triangulate();
    delaunay32::extras::write_mesh_svg(
        "delaunay32_quantization.svg",
        source_points,
        mesh.triangles,
        automatic.report);
    std::cout << "\nwrote delaunay32_quantization.svg\n";
    return 0;
}
