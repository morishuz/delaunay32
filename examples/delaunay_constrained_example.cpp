// SPDX-License-Identifier: MIT

#include "delaunay32/delaunay.hpp"
#include "delaunay32/extras/json.hpp"
#include "presentation_svg.hpp"

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
            argc == 3 ? argv[2] : "delaunay32_constrained.svg";
        const delaunay32::extras::Geometry geometry =
            delaunay32::extras::read_geometry_json(input_path);
        if (geometry.constraints.empty()) {
            throw std::invalid_argument(
                "constrained example requires at least one constraint");
        }
        if (geometry.polygon.has_value() || !geometry.polygons.empty()) {
            throw std::invalid_argument(
                "constrained example does not consume polygon rings");
        }

        delaunay32::Triangulator triangulator(0);
        const std::vector<delaunay32::Triangle> ordinary =
            triangulator.triangulate_int(geometry.points);
        const std::vector<delaunay32::Triangle> constrained =
            triangulator.triangulate_constrained_int(
                geometry.points, geometry.constraints);

        delaunay32_example::write_constrained_comparison_svg(
            output_path,
            geometry.points,
            geometry.constraints,
            ordinary,
            constrained);
        std::cout << "wrote " << output_path << ": "
                  << geometry.points.size() << " points, "
                  << geometry.constraints.size() << " constraints, "
                  << constrained.size() << " constrained triangles\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Constrained example error: " << error.what() << '\n'
                  << "Usage: " << argv[0]
                  << " input.json [output.svg]\n";
        return 1;
    }
}
