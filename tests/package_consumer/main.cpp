// SPDX-License-Identifier: MIT

#include <delaunay32/delaunay.hpp>

#include <vector>

int main() {
    const std::vector<delaunay32::Point> integer_points = {
        {0, 0},
        {100, 0},
        {100, 100},
        {0, 100},
        {50, 50},
    };
    delaunay32::Triangulator triangulator(1);
    if (triangulator.triangulate_int(integer_points).size() != 4) {
        return 1;
    }

    const std::vector<delaunay32::FloatPoint> float_points = {
        {0.0F, 0.0F},
        {1.0F, 0.0F},
        {1.0F, 1.0F},
        {0.0F, 1.0F},
        {0.5F, 0.5F},
    };
    delaunay32::QuantizationReport report;
    if (triangulator.triangulate_float(float_points, report).size() != 4 ||
        report.unique_points != float_points.size() ||
        report.collapsed_points != 0) {
        return 2;
    }
    return 0;
}
