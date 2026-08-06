// SPDX-License-Identifier: MIT

#include <delaunay32/delaunay.hpp>
#include <delaunay32/extras/json.hpp>
#include <delaunay32/extras/sampling.hpp>
#include <delaunay32/extras/svg.hpp>
#include <delaunay32/quantization.hpp>

#include <string>
#include <vector>

int main() {
    delaunay32::extras::PointSampler sampler;
    sampler.set_bounds({0.0, 100.0, 0.0, 100.0});
    delaunay32::extras::UniformSamplingOptions sample_options;
    sample_options.point_count = 8;
    const std::vector<delaunay32::FloatPoint> sampled_points =
        sampler.generate_uniform(sample_options);
    if (sampled_points.size() != sample_options.point_count) {
        return 7;
    }

    const std::vector<delaunay32::Point> integer_points = {
        {0, 0},
        {100, 0},
        {100, 100},
        {0, 100},
        {50, 50},
    };
    delaunay32::extras::Geometry json_geometry;
    json_geometry.points = integer_points;
    if (delaunay32::extras::geometry_to_json(json_geometry).find(
            "\"points\"") == std::string::npos) {
        return 9;
    }
    delaunay32::Triangulator triangulator;
    triangulator.set_points(integer_points);
    if (triangulator.triangulate().triangles.size() != 4) {
        return 1;
    }
    delaunay32::TriangulationOptions full_options;
    full_options.result_detail = delaunay32::ResultDetail::Full;
    triangulator.set_options(full_options);
    triangulator.set_points(integer_points);
    const delaunay32::TriangulationResult integer_result =
        triangulator.triangulate();
    if (integer_result.triangles.size() != 4 ||
        integer_result.halfedges.size() != 12 ||
        integer_result.hull.size() != 4 ||
        integer_result.representatives.size() != integer_points.size()) {
        return 2;
    }
    delaunay32::extras::Svg svg(320.0, 240.0);
    svg.draw_triangles(integer_points, integer_result.triangles);
    svg.draw_colored_triangles(integer_points, integer_result.triangles);
    svg.draw_points(integer_points);
    const std::string svg_text = svg.to_svg();
    if (svg_text.find("<svg") == std::string::npos ||
        svg_text.find("<path") == std::string::npos) {
        return 8;
    }
    const std::vector<delaunay32::Constraint> constraints = {{0, 2}};
    triangulator.set_options({});
    triangulator.set_points(integer_points);
    triangulator.set_constraints(constraints);
    if (triangulator.triangulate().triangles.size() != 4) {
        return 5;
    }
    triangulator.set_points(integer_points);
    triangulator.set_polygons({{{0, 1, 2, 3}, {}}});
    if (triangulator.triangulate().triangles.size() != 4) {
        return 6;
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
    const delaunay32::QuantizationResult grid =
        delaunay32::quantize(float_points, grid_options);
    triangulator.set_points(grid.points);
    if (triangulator.triangulate().triangles.size() != 4) {
        return 3;
    }
    const delaunay32::QuantizationResult automatic =
        delaunay32::quantize(float_points);
    triangulator.set_options(full_options);
    triangulator.set_points(automatic.points);
    const delaunay32::TriangulationResult float_result =
        triangulator.triangulate();
    if (float_result.triangles.size() != 4 ||
        automatic.report.unique_points != float_points.size() ||
        automatic.report.collapsed_points != 0) {
        return 4;
    }
    return 0;
}
