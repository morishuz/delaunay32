// SPDX-License-Identifier: MIT

#include "delaunay32/delaunay.hpp"
#include "delaunay32/extras/svg.hpp"

#include <iostream>
#include <vector>

int main() {
    const std::vector<delaunay32::Point> points = {
        {40, 40},
        {180, 30},
        {320, 60},
        {370, 180},
        {300, 310},
        {150, 340},
        {30, 260},
        {120, 120},
        {240, 110},
        {290, 220},
        {180, 270},
        {90, 210},
        {200, 190},
    };

    delaunay32::Triangulator triangulator;
    triangulator.set_points(points);
    const delaunay32::TriangulationResult result =
        triangulator.triangulate();

    const char* output_path = "delaunay32_hello_mesh.svg";
    delaunay32::extras::write_mesh_svg(
        output_path, points, result.triangles);

    std::cout << "wrote " << output_path << '\n';
    return 0;
}
