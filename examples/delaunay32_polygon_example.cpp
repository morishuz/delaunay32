// SPDX-License-Identifier: MIT

#include "delaunay32/delaunay.hpp"
#include "delaunay32/extras/svg.hpp"

#include <iostream>
#include <vector>

int main() {
    // Every polygon corner and every point that may be used by the mesh lives
    // in this one vector. Polygon rings refer to these points by index.
    const std::vector<delaunay32::Point> points = {
        {30, 50},    // 0  outer boundary starts here
        {360, 40},   // 1
        {410, 170},  // 2
        {350, 310},  // 3
        {190, 350},  // 4
        {35, 280},   // 5

        {145, 135},  // 6  irregular five-sided hole starts here
        {225, 105},  // 7
        {285, 165},  // 8
        {245, 235},  // 9
        {155, 220},  // 10

        // These points give the triangulation some detail around the hole.
        {85, 95},    // 11
        {205, 75},   // 12
        {325, 90},   // 13
        {90, 175},   // 14
        {330, 175},  // 15
        {90, 260},   // 16
        {210, 285},  // 17
        {315, 265},  // 18

        // These two points deliberately lie inside the hole. They remain
        // valid input points, but no output triangle will use them because
        // the hole is excluded from the polygon domain.
        {185, 165},  // 19
        {225, 195},  // 20

        // These points sit outside the outer boundary. Polygon clipping
        // excludes them for the same reason: only the area inside the outer
        // ring belongs to the domain. The SVG draws them hollow too.
        {-10, 120},  // 21
        {180, 10},   // 22
        {430, 80},   // 23
        {430, 270},  // 24
        {120, 370},  // 25
    };

    delaunay32::PolygonDomain domain;

    // A ring is just a list of point indices. The final edge back to the
    // first point is added automatically, so index 0 is not repeated here.
    domain.outer_ring = {0, 1, 2, 3, 4, 5};

    // A domain can contain zero or more holes. Each hole is another ring.
    // The clockwise/counterclockwise direction of either ring does not matter.
    domain.holes = {
        {6, 7, 8, 9, 10},
    };

    // set_polygons() recovers all ring edges, then keeps triangles inside the
    // outer ring and outside every hole. It accepts a vector because one run
    // can contain several separate polygon domains; this example needs one.
    delaunay32::Triangulator triangulator;
    triangulator.set_points(points);
    triangulator.set_polygons({domain});
    const delaunay32::TriangulationResult result =
        triangulator.triangulate();

    // This convenience exporter draws the triangles, both boundary rings,
    // and all input points. Points unused because they are inside the hole
    // are drawn hollow, making the clipping step visible.
    const char* output_path = "delaunay32_polygon.svg";
    delaunay32::extras::write_polygon_mesh_svg(
        output_path, points, domain, result.triangles);

    std::cout << "wrote " << output_path << '\n';
    return 0;
}
