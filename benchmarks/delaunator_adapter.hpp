// SPDX-License-Identifier: MIT

#pragma once

#include "delaunay32/delaunay.hpp"

#include <vector>

namespace delaunay32 {

class DelaunatorBaseline {
public:
    std::vector<Triangle> triangulate(const std::vector<Point>& points);
    void prepare(const std::vector<Point>& points);
    std::vector<Triangle> triangulate_prepared() const;

private:
    std::vector<double> coords_;
};

}  // namespace delaunay32
