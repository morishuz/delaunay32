// SPDX-License-Identifier: MIT

#include "delaunator_adapter.hpp"

#include <delaunator.hpp>

#include <limits>
#include <stdexcept>

namespace delaunay32 {

std::vector<Triangle> DelaunatorBaseline::triangulate(
    const std::vector<Point>& points) {
    prepare(points);
    return triangulate_prepared();
}

void DelaunatorBaseline::prepare(
    const std::vector<Point>& points) {
    if (points.size() < 3) {
        throw std::invalid_argument("need at least 3 points");
    }
    if (points.size() > std::numeric_limits<std::uint32_t>::max()) {
        throw std::invalid_argument("point count exceeds uint32 triangle indices");
    }

    coords_.clear();
    coords_.reserve(points.size() * 2);
    for (const Point& p : points) {
        coords_.push_back(p.x);
        coords_.push_back(p.y);
    }
}

std::vector<Triangle>
DelaunatorBaseline::triangulate_prepared() const {
    if (coords_.size() < 6) {
        throw std::logic_error("Delaunator input has not been prepared");
    }
    delaunator::Delaunator delaunay(coords_);
    std::vector<Triangle> triangles;
    triangles.reserve(delaunay.triangles.size() / 3);
    for (std::size_t i = 0; i < delaunay.triangles.size(); i += 3) {
        triangles.push_back(Triangle{
            static_cast<std::uint32_t>(delaunay.triangles[i]),
            static_cast<std::uint32_t>(delaunay.triangles[i + 1]),
            static_cast<std::uint32_t>(delaunay.triangles[i + 2]),
        });
    }
    return triangles;
}

}  // namespace delaunay32
