// SPDX-License-Identifier: MIT

#include "delaunay32/delaunay.hpp"
#include "geometry_io.hpp"
#include "svg_io.hpp"

#include <exception>
#include <iostream>
#include <stdexcept>
#include <vector>

int main(int argc, char** argv) {
    try {
        const delaunay32_example::Options options =
            delaunay32_example::parse_options(argc, argv);
        const delaunay32_example::GeometryInput geometry =
            delaunay32_example::load_geometry(options);
        if (!geometry.constraints.empty() ||
            !geometry.outer_ring.empty() || !geometry.holes.empty() ||
            !geometry.polygons.empty()) {
            throw std::invalid_argument(
                "ordinary SVG example expects a points-only JSON file");
        }
        const std::vector<delaunay32::Point>& points = geometry.points;

        delaunay32::Triangulator triangulator(0);
        const std::vector<delaunay32::Triangle> triangles =
            triangulator.triangulate_int(points);

        delaunay32_example::write_svg(
            options.output_path, points, triangles);
        std::cout << "wrote " << options.output_path << ": "
                  << points.size() << " points, "
                  << triangles.size() << " triangles\n";
        return 0;
    } catch (const std::exception& error) {
        delaunay32_example::print_usage(argv[0]);
        std::cerr << "SVG example error: " << error.what() << '\n';
        return 1;
    }
}
