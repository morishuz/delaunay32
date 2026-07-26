// SPDX-License-Identifier: MIT

#include "delaunay32/delaunay.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {

using delaunay32::Point;
using delaunay32::Triangle;
using delaunay32::Triangulator;

constexpr std::int32_t kDomainMaximum = 999;
constexpr double kCanvasSize = 1200.0;
constexpr double kHeaderHeight = 56.0;
constexpr double kPlotMargin = 30.0;
constexpr double kPlotSize =
    kCanvasSize - kHeaderHeight - 2.0 * kPlotMargin;
constexpr double kPlotLeft = (kCanvasSize - kPlotSize) / 2.0;
constexpr std::size_t kInteriorCapacity =
    static_cast<std::size_t>(kDomainMaximum - 1) *
    static_cast<std::size_t>(kDomainMaximum - 1);
constexpr std::size_t kMaximumPoints = kInteriorCapacity + 4;

std::int64_t orient(
    const Point& a,
    const Point& b,
    const Point& c) {
    return static_cast<std::int64_t>(b.x - a.x) * (c.y - a.y) -
           static_cast<std::int64_t>(b.y - a.y) * (c.x - a.x);
}

std::uint64_t point_key(const Point& point) {
    return (static_cast<std::uint64_t>(
                static_cast<std::uint32_t>(point.x))
            << 32U) |
           static_cast<std::uint32_t>(point.y);
}

std::uint64_t edge_key(std::uint32_t a, std::uint32_t b) {
    if (b < a) {
        std::swap(a, b);
    }
    return (static_cast<std::uint64_t>(a) << 32U) | b;
}

std::vector<Point> generate_points(
    std::size_t point_count,
    std::uint64_t seed) {
    if (point_count < 4 || point_count > kMaximumPoints) {
        throw std::invalid_argument(
            "point count must be between 4 and " +
            std::to_string(kMaximumPoints));
    }

    // These four sites make the convex hull exactly the requested square.
    std::vector<Point> points = {
        {0, 0},
        {kDomainMaximum, 0},
        {kDomainMaximum, kDomainMaximum},
        {0, kDomainMaximum},
    };
    points.reserve(point_count);

    std::unordered_set<std::uint64_t> occupied;
    occupied.reserve(point_count * 2);
    for (const Point& point : points) {
        occupied.insert(point_key(point));
    }

    std::mt19937_64 random(seed);
    std::uniform_int_distribution<std::int32_t> coordinate(
        1, kDomainMaximum - 1);
    while (points.size() < point_count) {
        const Point point = {coordinate(random), coordinate(random)};
        if (occupied.insert(point_key(point)).second) {
            points.push_back(point);
        }
    }
    return points;
}

void validate_square_mesh(
    const std::vector<Point>& points,
    const std::vector<Triangle>& triangles) {
    const std::size_t expected_triangles = 2 * points.size() - 6;
    if (triangles.size() != expected_triangles) {
        throw std::runtime_error(
            "mesh check failed: triangle count does not satisfy Euler's "
            "formula for a four-corner hull");
    }

    const std::array<std::uint64_t, 4> boundary_edges = {
        edge_key(0, 1),
        edge_key(1, 2),
        edge_key(2, 3),
        edge_key(3, 0),
    };
    std::unordered_map<std::uint64_t, unsigned> edge_counts;
    edge_counts.reserve(triangles.size() * 2);
    std::int64_t area_sum = 0;
    for (const Triangle& triangle : triangles) {
        if (triangle.i0 >= points.size() ||
            triangle.i1 >= points.size() ||
            triangle.i2 >= points.size()) {
            throw std::runtime_error(
                "mesh check failed: triangle index is out of range");
        }
        const std::int64_t area = orient(
            points[triangle.i0],
            points[triangle.i1],
            points[triangle.i2]);
        if (area <= 0) {
            throw std::runtime_error(
                "mesh check failed: triangle is degenerate or clockwise");
        }
        area_sum += area;
        ++edge_counts[edge_key(triangle.i0, triangle.i1)];
        ++edge_counts[edge_key(triangle.i1, triangle.i2)];
        ++edge_counts[edge_key(triangle.i2, triangle.i0)];
    }

    const std::int64_t expected_area =
        2LL * kDomainMaximum * kDomainMaximum;
    if (area_sum != expected_area) {
        throw std::runtime_error(
            "mesh check failed: triangle area does not fill the square");
    }

    for (const auto& [edge, count] : edge_counts) {
        const bool boundary =
            std::find(
                boundary_edges.begin(), boundary_edges.end(), edge) !=
            boundary_edges.end();
        if (count != (boundary ? 1U : 2U)) {
            throw std::runtime_error(
                "mesh check failed: edge incidence is not manifold");
        }
    }
    for (const std::uint64_t edge : boundary_edges) {
        if (edge_counts.find(edge) == edge_counts.end()) {
            throw std::runtime_error(
                "mesh check failed: a square boundary edge is absent");
        }
    }
}

double svg_x(std::int32_t x) {
    return kPlotLeft +
           static_cast<double>(x) * kPlotSize / kDomainMaximum;
}

double svg_y(std::int32_t y) {
    return kHeaderHeight + kPlotMargin +
           static_cast<double>(kDomainMaximum - y) * kPlotSize /
               kDomainMaximum;
}

void write_point(std::ostream& output, const Point& point) {
    output << std::fixed << std::setprecision(2)
           << svg_x(point.x) << ',' << svg_y(point.y);
}

void write_svg(
    const std::string& output_path,
    const std::vector<Point>& points,
    const std::vector<Triangle>& triangles) {
    std::ofstream output(output_path);
    if (!output) {
        throw std::runtime_error("could not create SVG: " + output_path);
    }

    static constexpr std::array<const char*, 7> colors = {
        "#dcefe8",
        "#dbe8f4",
        "#f3e3a6",
        "#e8dff0",
        "#cde8ee",
        "#e7ebc5",
        "#f1d9cf",
    };

    output
        << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
        << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"1200\" "
           "height=\"1200\" viewBox=\"0 0 1200 1200\">\n"
        << "<rect width=\"1200\" height=\"1200\" fill=\"#f7f7f5\"/>\n"
        << "<text x=\"30\" y=\"35\" font-family=\"system-ui, sans-serif\" "
           "font-size=\"18\" fill=\"#202426\">Delaunay32 | "
        << points.size() << " points | " << triangles.size()
        << " triangles | square coverage PASS</text>\n"
        << "<rect x=\"" << std::fixed << std::setprecision(2)
        << svg_x(0) << "\" y=\"" << svg_y(kDomainMaximum)
        << "\" width=\"" << kPlotSize << "\" height=\"" << kPlotSize
        << "\" fill=\"#e00024\"/>\n";

    for (const Triangle& triangle : triangles) {
        const std::size_t color =
            (static_cast<std::size_t>(triangle.i0) * 17 +
             static_cast<std::size_t>(triangle.i1) * 31 +
             static_cast<std::size_t>(triangle.i2) * 43) %
            colors.size();
        output << "<polygon fill=\"" << colors[color]
               << "\" stroke=\"#344044\" stroke-width=\"0.72\" "
                  "stroke-linejoin=\"round\" points=\"";
        write_point(output, points[triangle.i0]);
        output << ' ';
        write_point(output, points[triangle.i1]);
        output << ' ';
        write_point(output, points[triangle.i2]);
        output << "\"/>\n";
    }

    const double radius =
        points.size() <= 2000 ? 2.2 : (points.size() <= 20000 ? 1.1 : 0.55);
    for (const Point& point : points) {
        output << "<circle cx=\"" << std::fixed << std::setprecision(2)
               << svg_x(point.x) << "\" cy=\"" << svg_y(point.y)
               << "\" r=\"" << radius
               << "\" fill=\"#111719\" stroke=\"#ffffff\" "
                  "stroke-width=\"0.35\"/>\n";
    }
    output << "</svg>\n";
    if (!output) {
        throw std::runtime_error(
            "failed while writing SVG: " + output_path);
    }
}

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

int main(int argc, char** argv) {
    try {
        if (argc > 4) {
            std::cerr
                << "Usage: " << argv[0]
                << " [point-count] [output.svg] [seed]\n";
            return 2;
        }
        const std::size_t point_count =
            argc >= 2 ? parse_point_count(argv[1]) : 1000;
        const std::string output_path =
            argc >= 3 ? argv[2] : "mesh.svg";
        const std::uint64_t seed =
            argc >= 4 ? parse_seed(argv[3]) : 1;

        const std::vector<Point> points =
            generate_points(point_count, seed);
        Triangulator triangulator(0);
        const std::vector<Triangle> triangles =
            triangulator.triangulate(points);
        validate_square_mesh(points, triangles);
        write_svg(output_path, points, triangles);

        std::cout << "wrote " << output_path << ": " << points.size()
                  << " points, " << triangles.size()
                  << " triangles, exact square hull, validation passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "SVG example error: " << error.what() << '\n';
        return 1;
    }
}
