// SPDX-License-Identifier: MIT

#include "delaunay32/delaunay.hpp"
#include "delaunay32/extras/json.hpp"
#include "delaunay32/extras/sampling.hpp"
#include "delaunay32/extras/svg.hpp"
#include "delaunay32/quantization.hpp"

#include <cstddef>
#include <cstdint>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr std::size_t kBlueNoisePointCount = 625;
constexpr std::size_t kBestCandidateCount = 16;
constexpr std::uint64_t kSeed = 1;

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc > 3) {
            throw std::invalid_argument(
                "expected an optional outline.json and output.svg");
        }
        const std::string input_path =
            argc >= 2 ? argv[1] : DELAUNAY32_LOGO_INPUT_PATH;
        const std::string output_path =
            argc == 3 ? argv[2] : "delaunay32_logo_polygon.svg";

        // 1. Load the indexed letter outlines and holes from JSON.
        delaunay32::extras::Geometry geometry =
            delaunay32::extras::read_geometry_json(input_path);
        if (!geometry.constraints.empty() || geometry.polygon.has_value() ||
            geometry.polygons.empty()) {
            throw std::invalid_argument(
                "logo example expects points plus a polygons array");
        }

        const std::size_t outline_point_count = geometry.points.size();
        const std::size_t polygon_count = geometry.polygons.size();

        // 2. Generate well-spaced points inside the letters but outside holes.
        delaunay32::extras::PointSampler sampler;
        sampler.set_polygon_interiors(geometry.points, geometry.polygons);
        delaunay32::extras::BlueNoiseSamplingOptions sampling_options;
        sampling_options.point_count = kBlueNoisePointCount;
        sampling_options.candidates_per_point = kBestCandidateCount;
        sampling_options.seed = kSeed;
        const std::vector<delaunay32::FloatPoint> interior_points =
            sampler.generate_blue_noise(sampling_options);

        // 3. Quantize outlines and samples together so they share one mapping.
        // Keeping the same order preserves every polygon and triangle index.
        std::vector<delaunay32::FloatPoint> source_points;
        source_points.reserve(
            geometry.points.size() + interior_points.size());
        for (const delaunay32::Point& point : geometry.points) {
            source_points.push_back({
                static_cast<double>(point.x),
                static_cast<double>(point.y),
            });
        }
        source_points.insert(
            source_points.end(),
            interior_points.begin(),
            interior_points.end());
        const delaunay32::QuantizationResult quantized =
            delaunay32::quantize(source_points);

        // 4. Recover the letter boundaries and keep only triangles inside them.
        delaunay32::TriangulationOptions triangulation_options;
        triangulation_options.thread_count = 0;  // Let Delaunay32 choose.

        delaunay32::Triangulator triangulator;
        triangulator.set_options(triangulation_options);
        triangulator.set_points(quantized.points);
        triangulator.set_polygons(geometry.polygons);
        const delaunay32::TriangulationResult result =
            triangulator.triangulate();

        // 5. Draw with the original coordinates, using the computed indices.
        delaunay32::extras::Svg svg(1200.0, 440.0);
        svg.set_background("#ffffff00");  // Transparent background.
        svg.set_auto_fit({28.0, 68.0, 28.0, 28.0});

        // Edge-adjacent triangles are assigned different palette colors.
        delaunay32::extras::SvgTriangleColorStyle triangle_style;
        triangle_style.stroke = "#35566a";
        triangle_style.stroke_width = 0.3;
        triangle_style.palette = {
            "#ffefc4",
            "#aee3ff",
            "#58d5d5",
            "#6599af",
            "#4d466e",
        };
        svg.draw_colored_triangles(
            source_points, result.triangles, triangle_style);

        // Redraw the polygon rings so the letter edges and holes stay crisp.
        delaunay32::extras::SvgShapeStyle outline_style;
        outline_style.fill = "none";
        outline_style.stroke = "#475558";
        outline_style.stroke_width = 0.7;
        for (const delaunay32::PolygonDomain& domain : geometry.polygons) {
            svg.draw_polygon(source_points, domain, outline_style);
        }

        // Show the generated interior samples over the mesh.
        delaunay32::extras::SvgPointStyle interior_point_style;
        interior_point_style.fill = "#2479a6";
        interior_point_style.stroke = "#ffffff";
        interior_point_style.radius = 1.35;
        interior_point_style.stroke_width = 0.35;
        svg.draw_points(interior_points, interior_point_style);

        // Text uses canvas coordinates and is not included in auto-fitting.
        delaunay32::extras::SvgTextStyle title_style;
        title_style.fill = "#172326";
        title_style.font_size = 19.0;
        title_style.font_weight = "650";
        svg.draw_text(
            "Delaunay32 | polygon domains with blue-noise samples",
            28.0,
            34.0,
            title_style);
        svg.render_to_svg(output_path);

        std::cout << "wrote " << output_path << ": "
                  << interior_points.size()
                  << " blue-noise interior points, "
                  << outline_point_count << " outline points, "
                  << polygon_count << " polygon domains, "
                  << result.triangles.size() << " domain triangles, "
                  << result.report.actual_thread_count << " worker threads\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Logo polygon example error: " << error.what() << '\n'
                  << "Usage: " << argv[0]
                  << " [outline.json] [output.svg]\n";
        return 1;
    }
}
