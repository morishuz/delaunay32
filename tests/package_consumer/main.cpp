// SPDX-License-Identifier: MIT

#include <delaunay32/delaunay.hpp>

#include <vector>

int main() {
    const std::vector<delaunay32::Point> integer_points = {
        {0, 0},
        {100, 0},
        {100, 100},
        {0, 100},
        {50, 50},
    };
    delaunay32::Triangulator triangulator(1);
    if (triangulator.triangulate_int(integer_points).size() != 4) {
        return 1;
    }
    const delaunay32::TriangulationResult integer_result =
        triangulator.triangulate_int_full(integer_points);
    if (integer_result.triangles.size() != 4 ||
        integer_result.halfedges.size() != 12 ||
        integer_result.hull.size() != 4 ||
        integer_result.representatives.size() != integer_points.size()) {
        return 2;
    }

    const std::vector<delaunay32::FloatPoint> float_points = {
        {0.0F, 0.0F},
        {1.0F, 0.0F},
        {1.0F, 1.0F},
        {0.0F, 1.0F},
        {0.5F, 0.5F},
    };
    delaunay32::QuantizationOptions grid_options;
    grid_options.mode = delaunay32::QuantizationMode::GridStep;
    grid_options.grid_step = 0.125;
    if (triangulator.triangulate_float(
            float_points, grid_options).size() != 4) {
        return 3;
    }
    const delaunay32::TriangulationResult float_result =
        triangulator.triangulate_float_full(float_points);
    if (float_result.triangles.size() != 4 ||
        float_result.quantization.unique_points != float_points.size() ||
        float_result.quantization.collapsed_points != 0) {
        return 4;
    }
    return 0;
}
