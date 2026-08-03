// SPDX-License-Identifier: MIT

#include "delaunay32/delaunay.hpp"
#include "geometry_io.hpp"
#include "svg_io.hpp"

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
        const delaunay32_example::GeometryInput geometry =
            delaunay32_example::read_geometry_json(input_path);
        if (geometry.outer_ring.empty()) {
            throw std::invalid_argument(
                "polygon example requires polygon.outer");
        }

        delaunay32::Triangulator triangulator(0);
        const std::vector<delaunay32::Triangle> triangles =
            triangulator.triangulate_polygon_int(
                geometry.points,
                geometry.outer_ring,
                geometry.holes);
        delaunay32_example::write_polygon_svg(
            output_path,
            geometry.points,
            geometry.outer_ring,
            geometry.holes,
            triangles);
        std::cout << "wrote " << output_path << ": "
                  << geometry.points.size() << " points, "
                  << geometry.holes.size() << " holes, "
                  << triangles.size() << " domain triangles\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Polygon example error: " << error.what() << '\n'
                  << "Usage: " << argv[0]
                  << " input.json [output.svg]\n";
        return 1;
    }
}
