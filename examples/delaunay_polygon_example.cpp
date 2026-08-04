// SPDX-License-Identifier: MIT

#include "delaunay32/delaunay.hpp"
#include "delaunay32/extras/json.hpp"
#include "delaunay32/extras/svg.hpp"

#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    try {
        if (argc < 2 || argc > 3) {
            throw std::invalid_argument(
                "expected input.json and an optional output.svg");
        }
        const std::string input_path = argv[1];
        const std::string output_path =
            argc == 3 ? argv[2] : "delaunay32_polygon.svg";
        const delaunay32::extras::Geometry geometry =
            delaunay32::extras::read_geometry_json(input_path);
        if (!geometry.polygon.has_value()) {
            throw std::invalid_argument(
                "polygon example requires polygon.outer");
        }
        if (!geometry.constraints.empty()) {
            throw std::invalid_argument(
                "polygon example does not consume standalone constraints");
        }

        delaunay32::Triangulator triangulator(0);
        const std::vector<delaunay32::Triangle> triangles =
            triangulator.triangulate_polygon_int(
                geometry.points,
                geometry.polygon->outer_ring,
                geometry.polygon->holes);
        delaunay32::extras::write_polygon_mesh_svg(
            output_path,
            geometry.points,
            *geometry.polygon,
            triangles);
        std::cout << "wrote " << output_path << ": "
                  << geometry.points.size() << " points, "
                  << geometry.polygon->holes.size() << " holes, "
                  << triangles.size() << " domain triangles\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Polygon example error: " << error.what() << '\n'
                  << "Usage: " << argv[0]
                  << " input.json [output.svg]\n";
        return 1;
    }
}
