// SPDX-License-Identifier: MIT

#include "delaunay32/delaunay.hpp"
#include "geometry_io.hpp"
#include "svg_io.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

namespace {

using delaunay32::Point;
using delaunay32_example::PolygonDomain;

constexpr std::size_t blue_noise_point_count = 625;
constexpr std::size_t best_candidate_count = 16;

std::uint64_t point_key(const Point& point) {
    return (static_cast<std::uint64_t>(
                static_cast<std::uint32_t>(point.x))
            << 32U) |
           static_cast<std::uint32_t>(point.y);
}

long double squared_distance_to_segment(
    const Point& point,
    const Point& a,
    const Point& b) {
    const long double dx = static_cast<long double>(b.x) - a.x;
    const long double dy = static_cast<long double>(b.y) - a.y;
    const long double length_squared = dx * dx + dy * dy;
    if (length_squared == 0.0L) {
        const long double point_dx =
            static_cast<long double>(point.x) - a.x;
        const long double point_dy =
            static_cast<long double>(point.y) - a.y;
        return point_dx * point_dx + point_dy * point_dy;
    }
    const long double projection = std::clamp(
        ((static_cast<long double>(point.x) - a.x) * dx +
         (static_cast<long double>(point.y) - a.y) * dy) /
            length_squared,
        0.0L,
        1.0L);
    const long double nearest_x = a.x + projection * dx;
    const long double nearest_y = a.y + projection * dy;
    const long double point_dx = point.x - nearest_x;
    const long double point_dy = point.y - nearest_y;
    return point_dx * point_dx + point_dy * point_dy;
}

long double minimum_boundary_distance_squared(
    const Point& point,
    const PolygonDomain& domain,
    const std::vector<Point>& points) {
    long double minimum = std::numeric_limits<long double>::max();
    const auto inspect_ring = [&](const std::vector<std::uint32_t>& ring) {
        for (std::size_t i = 0; i < ring.size(); ++i) {
            minimum = std::min(
                minimum,
                squared_distance_to_segment(
                    point,
                    points[ring[i]],
                    points[ring[(i + 1) % ring.size()]]));
        }
    };
    inspect_ring(domain.outer_ring);
    for (const std::vector<std::uint32_t>& hole : domain.holes) {
        inspect_ring(hole);
    }
    return minimum;
}

struct Candidate {
    Point point;
    std::size_t domain = 0;
};

void append_blue_noise_interior_points(
    std::vector<Point>& points,
    const std::vector<PolygonDomain>& domains,
    std::size_t count) {
    if (points.empty() || domains.empty()) {
        throw std::invalid_argument(
            "logo geometry requires outline points and polygon domains");
    }

    std::int32_t minimum_x = std::numeric_limits<std::int32_t>::max();
    std::int32_t minimum_y = std::numeric_limits<std::int32_t>::max();
    std::int32_t maximum_x = std::numeric_limits<std::int32_t>::min();
    std::int32_t maximum_y = std::numeric_limits<std::int32_t>::min();
    std::unordered_set<std::uint64_t> occupied;
    occupied.reserve((points.size() + count) * 2);
    for (const Point& point : points) {
        minimum_x = std::min(minimum_x, point.x);
        minimum_y = std::min(minimum_y, point.y);
        maximum_x = std::max(maximum_x, point.x);
        maximum_y = std::max(maximum_y, point.y);
        occupied.insert(point_key(point));
    }

    std::random_device entropy;
    std::seed_seq seed = {
        entropy(), entropy(), entropy(), entropy(),
        entropy(), entropy(), entropy(), entropy(),
    };
    std::mt19937_64 random(seed);
    std::uniform_int_distribution<std::int32_t> x_coordinate(
        minimum_x, maximum_x);
    std::uniform_int_distribution<std::int32_t> y_coordinate(
        minimum_y, maximum_y);

    const std::size_t target_size = points.size() + count;
    points.reserve(target_size);
    std::vector<std::vector<Point>> accepted_by_domain(domains.size());
    while (points.size() < target_size) {
        Candidate best;
        long double best_distance_squared = -1.0L;
        for (std::size_t sample = 0;
             sample < best_candidate_count;
             ++sample) {
            Candidate candidate;
            bool found = false;
            for (std::size_t attempt = 0; attempt < 10000; ++attempt) {
                candidate.point = {
                    x_coordinate(random),
                    y_coordinate(random),
                };
                if (occupied.find(point_key(candidate.point)) !=
                    occupied.end()) {
                    continue;
                }
                for (std::size_t domain = 0;
                     domain < domains.size();
                     ++domain) {
                    if (delaunay32_example::point_is_strictly_inside_domain(
                            candidate.point, domains[domain], points)) {
                        candidate.domain = domain;
                        found = true;
                        break;
                    }
                }
                if (found) {
                    break;
                }
            }
            if (!found) {
                throw std::runtime_error(
                    "could not sample a point inside the logo domains");
            }

            long double minimum_distance_squared =
                minimum_boundary_distance_squared(
                    candidate.point,
                    domains[candidate.domain],
                    points);
            for (const Point& accepted :
                 accepted_by_domain[candidate.domain]) {
                const long double dx =
                    static_cast<long double>(candidate.point.x) - accepted.x;
                const long double dy =
                    static_cast<long double>(candidate.point.y) - accepted.y;
                minimum_distance_squared = std::min(
                    minimum_distance_squared,
                    dx * dx + dy * dy);
            }
            if (minimum_distance_squared > best_distance_squared) {
                best = candidate;
                best_distance_squared = minimum_distance_squared;
            }
        }

        occupied.insert(point_key(best.point));
        accepted_by_domain[best.domain].push_back(best.point);
        points.push_back(best.point);
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc < 2 || argc > 3) {
            throw std::invalid_argument(
                "expected outline.json and an optional output.svg");
        }
        const std::string input_path = argv[1];
        const std::string output_path =
            argc == 3 ? argv[2] : "delaunay32_logo_polygon.svg";
        delaunay32_example::GeometryInput geometry =
            delaunay32_example::read_geometry_json(input_path);
        if (!geometry.constraints.empty() ||
            !geometry.outer_ring.empty() || !geometry.holes.empty() ||
            geometry.polygons.empty()) {
            throw std::invalid_argument(
                "logo example expects points plus a polygons array");
        }

        const std::size_t outline_point_count = geometry.points.size();
        append_blue_noise_interior_points(
            geometry.points, geometry.polygons, blue_noise_point_count);

        delaunay32::Triangulator triangulator(0);
        std::vector<delaunay32::Triangle> triangles;
        for (const PolygonDomain& domain : geometry.polygons) {
            const std::vector<delaunay32::Triangle> part =
                triangulator.triangulate_polygon_int(
                    geometry.points,
                    domain.outer_ring,
                    domain.holes);
            triangles.insert(triangles.end(), part.begin(), part.end());
        }

        delaunay32_example::write_logo_polygon_svg(
            output_path,
            geometry.points,
            blue_noise_point_count,
            geometry.polygons,
            triangles);
        std::cout << "wrote " << output_path << ": "
                  << blue_noise_point_count
                  << " blue-noise interior points, "
                  << outline_point_count << " outline points, "
                  << geometry.polygons.size() << " polygon domains, "
                  << triangles.size() << " domain triangles\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Logo polygon example error: " << error.what() << '\n'
                  << "Usage: " << argv[0]
                  << " outline.json [output.svg]\n";
        return 1;
    }
}
