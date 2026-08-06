// SPDX-License-Identifier: MIT

#include "delaunay32/delaunay.hpp"
#include "delaunay32/extras/svg.hpp"

#include <iostream>
#include <vector>

int main() {
    // Delaunay32 triangulates integer points. The position of each point in
    // this vector is also its index: the first point is 0, the next is 1,
    // and so on. Constraints refer to points using these indices.
    const std::vector<delaunay32::Point> points = {
        {40, 40},    // 0
        {200, 25},   // 1
        {360, 50},   // 2
        {390, 160},  // 3
        {350, 290},  // 4
        {190, 330},  // 5
        {35, 270},   // 6
        {20, 150},   // 7
        {100, 90},   // 8  start of the constrained polyline
        {200, 80},   // 9
        {300, 105},  // 10
        {90, 210},   // 11
        {200, 190},  // 12
        {295, 230},  // 13 bend in the constrained polyline
    };

    // A constraint says that a particular line segment must appear in the
    // finished mesh. Here, 8 -> 13 -> 5 makes one bent polyline from two
    // constrained segments.
    //
    // Delaunay32 does not create new points where constraints cross. Two
    // constraints may meet at a shared endpoint, as they do here, but they
    // must not cross anywhere else.
    const std::vector<delaunay32::Constraint> constraints = {
        {8, 13},
        {13, 5},
    };

    // Configure one problem, then solve it. Apart from set_constraints(),
    // this is the same workflow as an ordinary Delaunay triangulation.
    delaunay32::Triangulator triangulator;
    triangulator.set_points(points);
    triangulator.set_constraints(constraints);
    const delaunay32::TriangulationResult result =
        triangulator.triangulate();

    // Draw the mesh first, then draw the constraints over it so the required
    // segments are easy to see in the SVG. Drawing does not change the mesh.
    delaunay32::extras::Svg svg(700.0, 520.0);
    svg.set_background("#f4f7f6");
    svg.set_auto_fit({40.0, 40.0, 40.0, 40.0});

    delaunay32::extras::SvgShapeStyle mesh_style;
    mesh_style.fill = "#dcebef";
    mesh_style.stroke = "#70878d";
    mesh_style.stroke_width = 1.0;
    svg.draw_triangles(points, result.triangles, mesh_style);

    delaunay32::extras::SvgLineStyle constraint_style;
    constraint_style.stroke = "#d85f45";
    constraint_style.stroke_width = 4.0;
    for (const delaunay32::Constraint& constraint : constraints) {
        svg.draw_line(
            points[constraint.i0], points[constraint.i1], constraint_style);
    }

    delaunay32::extras::SvgPointStyle point_style;
    point_style.fill = "#172b31";
    point_style.stroke = "#ffffff";
    point_style.radius = 4.0;
    point_style.stroke_width = 1.3;
    svg.draw_points(points, point_style);

    const char* output_path = "delaunay32_constraints.svg";
    svg.render_to_svg(output_path);

    std::cout << "wrote " << output_path << '\n';
    return 0;
}
