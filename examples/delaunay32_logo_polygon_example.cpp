// SPDX-License-Identifier: MIT

#include "delaunay32/delaunay.hpp"
#include "delaunay32/extras/json.hpp"
#include "delaunay32/extras/sampling.hpp"
#include "presentation_svg.hpp"

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
        if (argc < 2 || argc > 3) {
            throw std::invalid_argument(
                "expected outline.json and an optional output.svg");
        }
        const std::string input_path = argv[1];
        const std::string output_path =
            argc == 3 ? argv[2] : "delaunay32_logo_polygon.svg";
        delaunay32::extras::Geometry geometry =
            delaunay32::extras::read_geometry_json(input_path);
        if (!geometry.constraints.empty() || geometry.polygon.has_value() ||
            geometry.polygons.empty()) {
            throw std::invalid_argument(
                "logo example expects points plus a polygons array");
        }

        const std::size_t outline_point_count = geometry.points.size();
        delaunay32::extras::BestCandidateOptions sampling;
        sampling.point_count = kBlueNoisePointCount;
        sampling.candidates_per_point = kBestCandidateCount;
        sampling.seed = kSeed;
        const std::vector<delaunay32::Point> interior_points =
            delaunay32::extras::sample_polygon_interiors(
                geometry.points, geometry.polygons, sampling);
        geometry.points.insert(
            geometry.points.end(),
            interior_points.begin(),
            interior_points.end());

        delaunay32::Triangulator triangulator(0);
        std::vector<delaunay32::Triangle> triangles;
        for (const delaunay32::extras::PolygonDomain& domain :
             geometry.polygons) {
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
            interior_points.size(),
            geometry.polygons,
            triangles);
        std::cout << "wrote " << output_path << ": "
                  << interior_points.size()
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
