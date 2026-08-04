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
#include <unordered_set>
#include <vector>

namespace {

using delaunay32::Point;
using delaunay32::extras::Geometry;
using delaunay32::extras::PolygonDomain;

void expect(bool condition, const char* message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

std::uint64_t point_key(const Point& point) {
    return (static_cast<std::uint64_t>(
                static_cast<std::uint32_t>(point.x))
            << 32U) |
           static_cast<std::uint32_t>(point.y);
}

void validate_sampling() {
    delaunay32::extras::UniformIntOptions integer_options;
    integer_options.point_count = 128;
    integer_options.bounds = {-50, 70, -20, 90};
    integer_options.seed = 42;
    const std::vector<Point> integer_points =
        delaunay32::extras::generate_uniform_int_points(integer_options);
    expect(integer_points.size() == 128, "wrong integer sample count");
    expect(
        integer_points[0].x == -50 && integer_points[0].y == -20 &&
            integer_points[2].x == 70 && integer_points[2].y == 90,
        "integer sampling did not preserve the bounds corners");
    std::unordered_set<std::uint64_t> unique;
    for (const Point& point : integer_points) {
        unique.insert(point_key(point));
    }
    expect(unique.size() == integer_points.size(), "integer samples repeat");

    delaunay32::extras::UniformIntOptions dense_options;
    dense_options.point_count = 9;
    dense_options.bounds = {0, 2, 0, 2};
    dense_options.seed = 3;
    const std::vector<Point> dense =
        delaunay32::extras::generate_uniform_int_points(dense_options);
    unique.clear();
    for (const Point& point : dense) {
        unique.insert(point_key(point));
    }
    expect(unique.size() == 9, "dense integer sampling lost coordinates");

    delaunay32::extras::UniformFloatOptions float_options;
    float_options.point_count = 64;
    float_options.bounds = {-1.5F, 8.0F, 2.25F, 7.75F};
    float_options.seed = 11;
    expect(
        delaunay32::extras::generate_uniform_float_points(float_options)
                .size() == 64,
        "wrong float sample count");
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
    delaunay32::extras::write_geometry_json(json_path, geometry);
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
    delaunay32::extras::BestCandidateOptions options;
    options.point_count = sample_count;
    options.candidates_per_point = 8;
    options.seed = 7;
    const std::vector<Point> samples =
        delaunay32::extras::sample_polygon_interiors(
            geometry.points, {*geometry.polygon}, options);
    expect(samples.size() == sample_count, "wrong polygon sample count");
    for (const Point& point : samples) {
        expect(
            delaunay32::extras::point_is_strictly_inside_domain(
                point, *geometry.polygon, geometry.points),
            "polygon sampler returned a point outside the domain");
    }
    geometry.points.insert(
        geometry.points.end(), samples.begin(), samples.end());
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
        validate_polygon_sampling(geometry, 24);

        delaunay32::Triangulator triangulator(1);
        const std::vector<delaunay32::Triangle> ordinary =
            triangulator.triangulate_int(geometry.points);
        delaunay32::extras::write_mesh_svg(
            argv[2], geometry.points, ordinary);
        expect_svg(argv[2]);

        const std::vector<delaunay32::Triangle> polygon =
            triangulator.triangulate_polygon_int(
                geometry.points,
                geometry.polygon->outer_ring,
                geometry.polygon->holes);
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
