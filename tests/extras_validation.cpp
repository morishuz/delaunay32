// SPDX-License-Identifier: MIT

#include "delaunay32/delaunay.hpp"
#include "delaunay32/extras/geometry.hpp"
#include "delaunay32/extras/json.hpp"
#include "delaunay32/extras/sampling.hpp"
#include "delaunay32/extras/svg.hpp"

#include <cstddef>
#include <cstdint>
#include <exception>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using delaunay32::Point;
using delaunay32::Triangle;
using delaunay32::extras::Geometry;
using delaunay32::PolygonDomain;

void expect(bool condition, const char* message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

template <typename Function>
void expect_invalid_argument(Function function, const char* message) {
    bool rejected = false;
    try {
        function();
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    expect(rejected, message);
}

template <typename Function>
void expect_logic_error(Function function, const char* message) {
    bool rejected = false;
    try {
        function();
    } catch (const std::logic_error&) {
        rejected = true;
    }
    expect(rejected, message);
}

template <typename Function>
void expect_runtime_error(Function function, const char* message) {
    bool rejected = false;
    try {
        function();
    } catch (const std::runtime_error&) {
        rejected = true;
    }
    expect(rejected, message);
}

std::size_t count_substring(
    const std::string& text,
    const std::string& substring) {
    std::size_t count = 0;
    std::size_t position = 0;
    while ((position = text.find(substring, position)) != std::string::npos) {
        ++count;
        position += substring.size();
    }
    return count;
}

void validate_sampling() {
    delaunay32::extras::PointSampler sampler;
    expect_logic_error(
        [&] { sampler.generate_uniform(); },
        "sampling without a region was accepted");
    expect_logic_error(
        [&] { sampler.generate_jittered_grid(); },
        "jittered-grid sampling without a region was accepted");
    expect_invalid_argument(
        [&] {
            sampler.set_bounds({
                std::numeric_limits<double>::quiet_NaN(),
                1.0,
                0.0,
                1.0,
            });
        },
        "sampling accepted non-finite bounds");
    sampler.set_bounds({-50.0, 70.0, -20.0, 90.0});
    delaunay32::extras::UniformSamplingOptions uniform;
    uniform.point_count = 128;
    uniform.seed = 42;
    const std::vector<delaunay32::FloatPoint> first =
        sampler.generate_uniform(uniform);
    const std::vector<delaunay32::FloatPoint> repeated =
        sampler.generate_uniform(uniform);
    expect(first.size() == 128, "wrong uniform sample count");
    expect(
        first.front().x != -50.0 || first.front().y != -20.0,
        "uniform sampling injected a corner by default");
    for (std::size_t i = 0; i < first.size(); ++i) {
        expect(
            first[i].x == repeated[i].x && first[i].y == repeated[i].y,
            "uniform sampling is not deterministic");
        expect(
            first[i].x >= -50.0 && first[i].x <= 70.0 &&
                first[i].y >= -20.0 && first[i].y <= 90.0,
            "uniform sample is outside bounds");
    }

    delaunay32::extras::UniformSamplingOptions corners;
    corners.point_count = 4;
    corners.include_bounds_corners = true;
    const std::vector<delaunay32::FloatPoint> with_corners =
        sampler.generate_uniform(corners);
    expect(
        with_corners[0].x == -50.0 && with_corners[0].y == -20.0 &&
            with_corners[2].x == 70.0 && with_corners[2].y == 90.0,
        "explicit bounds corners were not preserved");

    delaunay32::extras::BlueNoiseSamplingOptions blue_noise;
    blue_noise.point_count = 32;
    blue_noise.candidates_per_point = 8;
    blue_noise.seed = 11;
    const std::vector<delaunay32::FloatPoint> blue_first =
        sampler.generate_blue_noise(blue_noise);
    const std::vector<delaunay32::FloatPoint> blue_repeated =
        sampler.generate_blue_noise(blue_noise);
    expect(
        blue_first.size() == 32,
        "wrong bounds blue-noise sample count");
    for (std::size_t i = 0; i < blue_first.size(); ++i) {
        expect(
            blue_first[i].x == blue_repeated[i].x &&
                blue_first[i].y == blue_repeated[i].y,
            "blue-noise sampling is not deterministic");
        expect(
            blue_first[i].x >= -50.0 && blue_first[i].x <= 70.0 &&
                blue_first[i].y >= -20.0 && blue_first[i].y <= 90.0,
            "blue-noise sample is outside bounds");
    }

    delaunay32::extras::JitteredGridSamplingOptions jittered;
    jittered.point_count = 32;
    jittered.jitter = 0.1;
    jittered.seed = 19;
    const std::vector<delaunay32::FloatPoint> jittered_first =
        sampler.generate_jittered_grid(jittered);
    const std::vector<delaunay32::FloatPoint> jittered_repeated =
        sampler.generate_jittered_grid(jittered);
    jittered.seed = 20;
    const std::vector<delaunay32::FloatPoint> jittered_other_seed =
        sampler.generate_jittered_grid(jittered);
    expect(
        jittered_first.size() == 32,
        "wrong bounds jittered-grid sample count");
    for (std::size_t i = 0; i < jittered_first.size(); ++i) {
        expect(
            jittered_first[i].x == jittered_repeated[i].x &&
                jittered_first[i].y == jittered_repeated[i].y,
            "jittered-grid sampling is not deterministic");
        expect(
            jittered_first[i].x >= -50.0 && jittered_first[i].x <= 70.0 &&
                jittered_first[i].y >= -20.0 && jittered_first[i].y <= 90.0,
            "jittered-grid sample is outside bounds");
    }
    expect(
        jittered_first.front().x != jittered_other_seed.front().x ||
            jittered_first.front().y != jittered_other_seed.front().y,
        "jittered-grid seed did not change rotation, phase, or jitter");
    for (const double jitter : {0.0, 1.0}) {
        for (const std::size_t point_count : {1U, 2U, 31U, 33U, 127U}) {
            jittered.point_count = point_count;
            jittered.jitter = jitter;
            expect(
                sampler.generate_jittered_grid(jittered).size() == point_count,
                "jittered-grid count reconciliation failed");
        }
    }
    jittered.point_count = 32;
    jittered.jitter = 0.1;
    jittered.attempts_per_point = 1;
    expect_runtime_error(
        [&] { sampler.generate_jittered_grid(jittered); },
        "jittered-grid sampling ignored its candidate work limit");
    jittered.jitter = 1.1;
    jittered.attempts_per_point = 10000;
    expect_invalid_argument(
        [&] { sampler.generate_jittered_grid(jittered); },
        "jittered-grid sampling accepted jitter above one");
    jittered.jitter = 0.1;
    jittered.attempts_per_point = 0;
    expect_invalid_argument(
        [&] { sampler.generate_jittered_grid(jittered); },
        "jittered-grid sampling accepted zero attempts");

    const std::vector<Point> polygon_points = {
        {0, 0},
        {10, 0},
        {10, 10},
        {0, 10},
    };
    const PolygonDomain polygon{{0, 1, 2, 3}, {}};
    sampler.set_polygon_interiors(polygon_points, {polygon});
    uniform.point_count = 8;
    uniform.include_bounds_corners = true;
    expect_invalid_argument(
        [&] { sampler.generate_uniform(uniform); },
        "polygon sampling accepted bounds corner injection");
    uniform.include_bounds_corners = false;
    uniform.attempts_per_point = 0;
    expect_invalid_argument(
        [&] { sampler.generate_uniform(uniform); },
        "uniform sampling accepted zero attempts");
    expect_invalid_argument(
        [&] {
            delaunay32::extras::PointSampler invalid;
            invalid.set_polygon_interiors(
                polygon_points, std::vector<PolygonDomain>{});
        },
        "polygon sampling accepted no domains");

    sampler.set_bounds({0.0, 1.0, 0.0, 1.0});
    uniform.attempts_per_point = 10000;
    uniform.point_count = 4;
    expect(
        sampler.generate_uniform(uniform).size() == 4,
        "setting bounds did not replace the polygon region");
}

void validate_extreme_domain_query() {
    const std::int32_t low = std::numeric_limits<std::int32_t>::min();
    const std::int32_t high = std::numeric_limits<std::int32_t>::max();
    const std::vector<Point> points = {
        {low, low},
        {high, low},
        {high, high},
        {low, high},
    };
    const PolygonDomain domain{{0, 1, 2, 3}, {}};
    expect(
        delaunay32::extras::point_is_strictly_inside_domain(
            {0, 0}, domain, points),
        "extreme signed domain rejected an interior point");
    expect(
        !delaunay32::extras::point_is_strictly_inside_domain(
            {low, 0}, domain, points),
        "extreme signed domain accepted a boundary point");

    std::vector<delaunay32::FloatPoint> converted;
    for (const Point& point : points) {
        converted.push_back({
            static_cast<double>(point.x),
            static_cast<double>(point.y),
        });
    }
    expect(
        converted[0].x == static_cast<double>(low) &&
            converted[2].x == static_cast<double>(high),
        "FloatPoint did not exactly retain int32 coordinates");
    delaunay32::extras::PointSampler sampler;
    sampler.set_polygon_interiors(points, {domain});
    delaunay32::extras::UniformSamplingOptions options;
    options.point_count = 4;
    options.seed = 17;
    for (const delaunay32::FloatPoint& sample :
         sampler.generate_uniform(options)) {
        expect(
            delaunay32::extras::point_is_strictly_inside_domain(
                sample, domain, converted),
            "integer polygon overload lost an extreme coordinate");
    }
}

Geometry make_geometry() {
    Geometry geometry;
    geometry.points = {
        {0, 0},
        {100, 0},
        {100, 100},
        {0, 100},
        {40, 40},
        {60, 40},
        {60, 60},
        {40, 60},
    };
    geometry.constraints = {{0, 1}};
    geometry.polygon = PolygonDomain{
        {0, 1, 2, 3},
        {{4, 5, 6, 7}},
    };
    return geometry;
}

void validate_json_round_trip(
    const std::string& json_path,
    const Geometry& geometry) {
    const std::string encoded =
        delaunay32::extras::geometry_to_json(geometry);
    expect(
        encoded.find("\"points\"") != std::string::npos &&
            encoded.find("\"constraints\"") != std::string::npos &&
            encoded.find("\"polygon\"") != std::string::npos,
        "in-memory JSON omitted geometry fields");
    delaunay32::extras::write_geometry_json(json_path, geometry);
    std::ifstream written_input(json_path, std::ios::binary);
    const std::string written{
        std::istreambuf_iterator<char>(written_input),
        std::istreambuf_iterator<char>(),
    };
    expect(
        written == encoded,
        "file and in-memory JSON serialization differ");
    const Geometry restored =
        delaunay32::extras::read_geometry_json(json_path);
    expect(
        restored.points.size() == geometry.points.size(),
        "JSON changed point count");
    expect(
        restored.constraints.size() == geometry.constraints.size(),
        "JSON changed constraint count");
    expect(restored.polygon.has_value(), "JSON lost polygon");
    expect(
        restored.polygon->outer_ring == geometry.polygon->outer_ring &&
            restored.polygon->holes == geometry.polygon->holes,
        "JSON changed polygon rings");

    Geometry multi_domain = geometry;
    multi_domain.polygons = {*multi_domain.polygon};
    multi_domain.polygon.reset();
    delaunay32::extras::write_geometry_json(json_path, multi_domain);
    const Geometry restored_multi =
        delaunay32::extras::read_geometry_json(json_path);
    expect(
        !restored_multi.polygon.has_value() &&
            restored_multi.polygons.size() == 1 &&
            restored_multi.polygons[0].holes ==
                multi_domain.polygons[0].holes,
        "JSON changed multiple polygon domains");
}

void validate_polygon_sampling(
    Geometry& geometry,
    std::size_t sample_count) {
    delaunay32::extras::PointSampler sampler;
    sampler.set_polygon_interiors(
        geometry.points, {*geometry.polygon});
    delaunay32::extras::BlueNoiseSamplingOptions options;
    options.point_count = sample_count;
    options.candidates_per_point = 8;
    options.seed = 7;
    const std::vector<delaunay32::FloatPoint> samples =
        sampler.generate_blue_noise(options);
    expect(samples.size() == sample_count, "wrong polygon sample count");

    delaunay32::extras::JitteredGridSamplingOptions jittered_options;
    jittered_options.point_count = sample_count;
    jittered_options.seed = 23;
    const std::vector<delaunay32::FloatPoint> jittered_samples =
        sampler.generate_jittered_grid(jittered_options);
    expect(
        jittered_samples.size() == sample_count,
        "wrong polygon jittered-grid sample count");

    std::vector<delaunay32::FloatPoint> source;
    source.reserve(geometry.points.size() + samples.size());
    for (const Point& point : geometry.points) {
        source.push_back({
            static_cast<double>(point.x),
            static_cast<double>(point.y),
        });
    }
    for (const delaunay32::FloatPoint& point : samples) {
        expect(
            delaunay32::extras::point_is_strictly_inside_domain(
                point, *geometry.polygon, source),
            "polygon sampler returned a point outside the domain");
    }
    for (const delaunay32::FloatPoint& point : jittered_samples) {
        expect(
            delaunay32::extras::point_is_strictly_inside_domain(
                point, *geometry.polygon, source),
            "polygon jittered-grid sampler returned a point outside the domain");
    }
    source.insert(source.end(), samples.begin(), samples.end());
    geometry.points = delaunay32::quantize(source).points;

    std::vector<delaunay32::FloatPoint> floating_boundary(
        source.begin(), source.begin() + 8);
    delaunay32::extras::PointSampler floating_sampler;
    floating_sampler.set_polygon_interiors(
        std::move(floating_boundary), {*geometry.polygon});
    delaunay32::extras::UniformSamplingOptions uniform;
    uniform.point_count = 16;
    uniform.seed = 13;
    const std::vector<delaunay32::FloatPoint> uniform_samples =
        floating_sampler.generate_uniform(uniform);
    expect(uniform_samples.size() == 16, "wrong polygon uniform count");
    for (const delaunay32::FloatPoint& point : uniform_samples) {
        expect(
            delaunay32::extras::point_is_strictly_inside_domain(
                point,
                *geometry.polygon,
                std::vector<delaunay32::FloatPoint>(
                    source.begin(), source.begin() + 8)),
            "uniform polygon sampler returned a point outside the domain");
    }
}

void expect_svg(const std::string& path) {
    std::ifstream input(path, std::ios::binary);
    expect(static_cast<bool>(input), "SVG was not created");
    const std::string text{
        std::istreambuf_iterator<char>(input),
        std::istreambuf_iterator<char>(),
    };
    expect(text.find("<svg") != std::string::npos, "SVG has no root node");
}

void validate_svg_builder(const Geometry& geometry) {
    delaunay32::extras::Svg svg(320.0, 240.0);
    svg.set_background("#ffffff");
    svg.set_auto_fit({20.0, 40.0, 20.0, 20.0});

    delaunay32::extras::SvgShapeStyle domain_style;
    domain_style.fill = "#ddeeff";
    domain_style.stroke = "#123456";
    domain_style.stroke_width = 2.5;
    svg.draw_polygon(geometry.points, *geometry.polygon, domain_style);
    svg.draw_triangles(
        geometry.points,
        std::vector<delaunay32::Triangle>{{0, 1, 2}},
        domain_style);

    delaunay32::extras::SvgLineStyle line_style;
    line_style.stroke = "red";
    line_style.stroke_width = 3.0;
    svg.draw_line(geometry.points[0], geometry.points[2], line_style);

    delaunay32::extras::SvgPointStyle point_style;
    point_style.fill = "blue";
    point_style.radius = 4.0;
    svg.draw_points(geometry.points, point_style);
    svg.draw_points(
        std::vector<delaunay32::FloatPoint>{{0.5F, 0.25F}}, point_style);
    svg.draw_text("mesh <&>", 12.0, 24.0);

    const std::string generated = svg.to_svg();
    expect(
        generated.find("fill-rule=\"evenodd\"") != std::string::npos,
        "SVG polygon did not use the even-odd fill rule");
    expect(
        generated.find("<line") != std::string::npos,
        "SVG line primitive is missing");
    expect(
        generated.find("mesh &lt;&amp;&gt;") != std::string::npos,
        "SVG text was not escaped");

    delaunay32::extras::Svg transformed(100.0, 100.0);
    transformed.set_transform(2.0, -3.0, 5.0, 90.0);
    transformed.draw_point(10.0, 20.0);
    const std::string exact = transformed.to_svg();
    expect(
        exact.find("cx=\"25.000\" cy=\"30.000\"") !=
            std::string::npos,
        "SVG explicit transform was not applied");

    bool invalid_index_rejected = false;
    try {
        svg.draw_polygon(geometry.points, {0, 1, 999});
    } catch (const std::invalid_argument&) {
        invalid_index_rejected = true;
    }
    expect(invalid_index_rejected, "SVG accepted an invalid polygon index");
}

void validate_colored_triangles() {
    const std::vector<Point> points = {
        {0, 0},
        {10, 0},
        {0, 10},
        {10, 10},
        {-10, 0},
        {0, -10},
    };
    const std::vector<Triangle> triangles = {
        {0, 1, 2},
        {1, 3, 2},
        {0, 4, 5},
    };
    delaunay32::extras::SvgTriangleColorStyle style;
    style.palette = {
        "#aa000080",
        "rgb(0 170 0 / 50%)",
    };
    style.stroke = "#0000ff80";
    style.stroke_width = 1.5;

    delaunay32::extras::Svg colored(200.0, 200.0);
    colored.set_background("#ffffff80");
    colored.set_transform(1.0, 1.0);
    colored.draw_colored_triangles(points, triangles, style);
    const std::string generated = colored.to_svg();
    expect(
        count_substring(generated, "fill=\"#aa000080\"") == 1 &&
            count_substring(
                generated, "fill=\"rgb(0 170 0 / 50%)\"") == 1,
        "SVG edge-neighbor coloring did not use both palette colors");
    expect(
        generated.find(
            "M 0.000 0.000 L 10.000 0.000 L 0.000 10.000 Z M "
            "0.000 0.000 L -10.000 0.000 L 0.000 -10.000 Z") !=
            std::string::npos,
        "SVG vertex-only neighbors did not retain their shared color");
    expect(
        generated.find("fill=\"#ffffff80\"") != std::string::npos &&
            generated.find("stroke=\"#0000ff80\"") !=
                std::string::npos,
        "SVG alpha-bearing color strings were not preserved");

    delaunay32::extras::Svg repeated(200.0, 200.0);
    repeated.set_background("#ffffff80");
    repeated.set_transform(1.0, 1.0);
    repeated.draw_colored_triangles(points, triangles, style);
    expect(
        repeated.to_svg() == generated,
        "SVG edge-neighbor coloring is not deterministic");

    std::vector<delaunay32::FloatPoint> float_points;
    float_points.reserve(points.size());
    for (const Point& point : points) {
        float_points.push_back({
            static_cast<double>(point.x),
            static_cast<double>(point.y),
        });
    }
    delaunay32::extras::Svg float_svg(200.0, 200.0);
    float_svg.draw_colored_triangles(float_points, triangles, style);
    expect(
        float_svg.to_svg().find("rgb(0 170 0 / 50%)") !=
            std::string::npos,
        "SVG colored triangle float overload lost the palette");

    const std::vector<Triangle> color_cycle = {
        {0, 1, 2},
        {0, 2, 3},
        {0, 3, 1},
    };
    delaunay32::extras::SvgTriangleColorStyle empty_palette;
    empty_palette.palette.clear();
    expect_invalid_argument(
        [&] {
            delaunay32::extras::Svg invalid(100.0, 100.0);
            invalid.draw_colored_triangles(
                points, color_cycle, empty_palette);
        },
        "SVG accepted an empty triangle color palette");

    delaunay32::extras::SvgTriangleColorStyle duplicate_palette;
    duplicate_palette.palette = {"red", "red"};
    expect_invalid_argument(
        [&] {
            delaunay32::extras::Svg invalid(100.0, 100.0);
            invalid.draw_colored_triangles(
                points, color_cycle, duplicate_palette);
        },
        "SVG accepted duplicate triangle palette colors");

    const std::vector<Triangle> four_clique = {
        {0, 1, 2},
        {0, 3, 1},
        {1, 3, 2},
        {0, 2, 3},
    };
    delaunay32::extras::SvgTriangleColorStyle insufficient_palette;
    insufficient_palette.palette = {"red", "blue", "green"};
    std::string insufficient_message;
    try {
        delaunay32::extras::Svg invalid(100.0, 100.0);
        invalid.draw_colored_triangles(
            points, four_clique, insufficient_palette);
    } catch (const std::invalid_argument& error) {
        insufficient_message = error.what();
    }
    expect(
        insufficient_message.find("contains 3 colors") !=
                std::string::npos &&
            insufficient_message.find("at least 4 unique colors") !=
                std::string::npos,
        "SVG insufficient palette error did not explain the four-color "
        "guarantee");

    const std::vector<Triangle> non_manifold = {
        {0, 1, 2},
        {1, 0, 3},
        {0, 1, 4},
    };
    expect_invalid_argument(
        [&] {
            delaunay32::extras::Svg invalid(100.0, 100.0);
            invalid.draw_colored_triangles(points, non_manifold);
        },
        "SVG accepted a non-manifold colored triangle edge");

    delaunay32::extras::Svg four_colors(100.0, 100.0);
    delaunay32::extras::SvgTriangleColorStyle guaranteed_palette;
    guaranteed_palette.palette = {"red", "blue", "green", "yellow"};
    four_colors.draw_colored_triangles(
        points, four_clique, guaranteed_palette);
    expect(
        four_colors.to_svg().find("<path") != std::string::npos,
        "SVG four-color triangle palette did not render");
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 4) {
            throw std::invalid_argument(
                "expected JSON, mesh SVG, and polygon SVG output paths");
        }
        validate_sampling();
        validate_extreme_domain_query();
        Geometry geometry = make_geometry();
        validate_json_round_trip(argv[1], geometry);
        validate_svg_builder(geometry);
        validate_colored_triangles();
        validate_polygon_sampling(geometry, 24);

        delaunay32::Triangulator triangulator;
        triangulator.set_points(geometry.points);
        const std::vector<delaunay32::Triangle> ordinary =
            triangulator.triangulate().triangles;
        delaunay32::extras::write_mesh_svg(
            argv[2], geometry.points, ordinary);
        expect_svg(argv[2]);

        triangulator.set_points(geometry.points);
        triangulator.set_polygons({*geometry.polygon});
        const std::vector<delaunay32::Triangle> polygon =
            triangulator.triangulate().triangles;
        delaunay32::extras::write_polygon_mesh_svg(
            argv[3], geometry.points, *geometry.polygon, polygon);
        expect_svg(argv[3]);
        std::cout << "extras validation passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "extras validation failed: " << error.what() << '\n';
        return 1;
    }
}
