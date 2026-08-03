// SPDX-License-Identifier: MIT

#pragma once

#include "delaunay32/delaunay.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace delaunay32::benchmark_support {

enum class Dataset {
    Uniform,
    Clustered,
    Diagonal,
};

inline const char* dataset_name(Dataset dataset) {
    switch (dataset) {
        case Dataset::Uniform:
            return "uniform";
        case Dataset::Clustered:
            return "clustered";
        case Dataset::Diagonal:
            return "diagonal";
    }
    return "unknown";
}

inline std::uint64_t point_key(std::int32_t x, std::int32_t y) {
    return (static_cast<std::uint64_t>(
                static_cast<std::uint32_t>(x))
            << 32U) |
           static_cast<std::uint32_t>(y);
}

inline std::vector<Point> generate_points(
    Dataset dataset,
    std::size_t count,
    std::uint64_t seed,
    std::int32_t domain = 100000) {
    if (domain < 2) {
        throw std::invalid_argument(
            "point count exceeds the integer domain");
    }
    const std::uint64_t unsigned_domain =
        static_cast<std::uint64_t>(domain);
    const std::uint64_t capacity = unsigned_domain * unsigned_domain;
    if (count > capacity) {
        throw std::invalid_argument(
            "point count exceeds the integer domain");
    }

    std::mt19937_64 random(seed);
    std::uniform_int_distribution<std::int32_t> coordinate(0, domain - 1);
    std::uniform_int_distribution<int> cluster_index(0, 7);
    std::uniform_int_distribution<std::int32_t> cluster_offset(
        -domain / 20,
        domain / 20);
    std::uniform_int_distribution<int> band_index(0, 3);
    std::uniform_int_distribution<std::int32_t> diagonal_jitter(
        -domain / 100,
        domain / 100);
    const std::array<std::array<std::int32_t, 2>, 8> centers = {{
        {domain / 8, domain / 8},
        {3 * domain / 8, domain / 8},
        {5 * domain / 8, domain / 8},
        {7 * domain / 8, domain / 8},
        {domain / 8, 7 * domain / 8},
        {3 * domain / 8, 7 * domain / 8},
        {5 * domain / 8, 7 * domain / 8},
        {7 * domain / 8, 7 * domain / 8},
    }};
    const std::array<std::int32_t, 4> diagonal_offsets = {
        domain / 9,
        3 * domain / 9,
        5 * domain / 9,
        7 * domain / 9,
    };

    std::vector<Point> points;
    points.reserve(count);
    std::unordered_set<std::uint64_t> occupied;
    occupied.reserve(count * 2);
    while (points.size() < count) {
        Point point;
        if (dataset == Dataset::Uniform) {
            point = {coordinate(random), coordinate(random)};
        } else if (dataset == Dataset::Clustered) {
            const auto& center = centers[
                static_cast<std::size_t>(cluster_index(random))];
            point.x = std::clamp(
                center[0] + cluster_offset(random),
                std::int32_t{0},
                static_cast<std::int32_t>(domain - 1));
            point.y = std::clamp(
                center[1] + cluster_offset(random),
                std::int32_t{0},
                static_cast<std::int32_t>(domain - 1));
        } else {
            point.x = coordinate(random);
            const std::int64_t raw_y =
                static_cast<std::int64_t>(point.x) +
                diagonal_offsets[
                    static_cast<std::size_t>(band_index(random))] +
                diagonal_jitter(random);
            point.y = static_cast<std::int32_t>(
                (raw_y % domain + domain) % domain);
        }
        if (occupied.insert(point_key(point.x, point.y)).second) {
            points.push_back(point);
        }
    }
    return points;
}

inline Triangle canonical_triangle(Triangle triangle) {
    if (triangle.i1 < triangle.i0) {
        std::swap(triangle.i0, triangle.i1);
    }
    if (triangle.i2 < triangle.i0) {
        std::swap(triangle.i0, triangle.i2);
    }
    if (triangle.i2 < triangle.i1) {
        std::swap(triangle.i1, triangle.i2);
    }
    return triangle;
}

inline bool meshes_equal(
    std::vector<Triangle> a,
    std::vector<Triangle> b) {
    const auto canonicalize = [](std::vector<Triangle>& triangles) {
        for (Triangle& triangle : triangles) {
            triangle = canonical_triangle(triangle);
        }
        std::sort(
            triangles.begin(),
            triangles.end(),
            [](const Triangle& lhs, const Triangle& rhs) {
                if (lhs.i0 != rhs.i0) {
                    return lhs.i0 < rhs.i0;
                }
                if (lhs.i1 != rhs.i1) {
                    return lhs.i1 < rhs.i1;
                }
                return lhs.i2 < rhs.i2;
            });
    };
    canonicalize(a);
    canonicalize(b);
    return a.size() == b.size() &&
           std::equal(
               a.begin(),
               a.end(),
               b.begin(),
               [](const Triangle& lhs, const Triangle& rhs) {
                   return lhs.i0 == rhs.i0 &&
                          lhs.i1 == rhs.i1 &&
                          lhs.i2 == rhs.i2;
               });
}

#if defined(__SIZEOF_INT128__)
using Exact = __int128_t;
#else
using Exact = long double;
#endif

inline Exact orient(
    const Point& a,
    const Point& b,
    const Point& c) {
    const Exact abx = static_cast<Exact>(b.x) - a.x;
    const Exact aby = static_cast<Exact>(b.y) - a.y;
    const Exact acx = static_cast<Exact>(c.x) - a.x;
    const Exact acy = static_cast<Exact>(c.y) - a.y;
    return abx * acy - aby * acx;
}

inline bool inside_circumcircle(
    const Point& a,
    const Point& b,
    const Point& c,
    const Point& d) {
    const Exact ax = static_cast<Exact>(a.x) - d.x;
    const Exact ay = static_cast<Exact>(a.y) - d.y;
    const Exact bx = static_cast<Exact>(b.x) - d.x;
    const Exact by = static_cast<Exact>(b.y) - d.y;
    const Exact cx = static_cast<Exact>(c.x) - d.x;
    const Exact cy = static_cast<Exact>(c.y) - d.y;
    const Exact a2 = ax * ax + ay * ay;
    const Exact b2 = bx * bx + by * by;
    const Exact c2 = cx * cx + cy * cy;
    const Exact determinant =
        ax * (by * c2 - b2 * cy) -
        ay * (bx * c2 - b2 * cx) +
        a2 * (bx * cy - by * cx);
    const Exact winding = orient(a, b, c);
    return winding > 0 ? determinant > 0 : determinant < 0;
}

struct EdgeRecord {
    std::uint32_t a = 0;
    std::uint32_t b = 0;
    std::uint32_t opposite = 0;
    unsigned count = 0;
};

inline std::uint64_t edge_key(std::uint32_t a, std::uint32_t b) {
    if (b < a) {
        std::swap(a, b);
    }
    return (static_cast<std::uint64_t>(a) << 32U) | b;
}

inline bool validate_mesh(
    const std::vector<Point>& points,
    const std::vector<Triangle>& triangles,
    std::string& error) {
    std::vector<bool> used(points.size(), false);
    std::vector<Triangle> canonical;
    canonical.reserve(triangles.size());
    std::unordered_map<std::uint64_t, EdgeRecord> edges;
    edges.reserve(triangles.size() * 2);

    const auto add_edge =
        [&](std::uint32_t a,
            std::uint32_t b,
            std::uint32_t opposite) -> bool {
            const std::uint64_t key = edge_key(a, b);
            auto [iterator, inserted] =
                edges.emplace(key, EdgeRecord{a, b, opposite, 1});
            if (inserted) {
                return true;
            }
            EdgeRecord& edge = iterator->second;
            if (edge.count != 1) {
                error = "mesh contains a non-manifold edge";
                return false;
            }
            edge.count = 2;
            if (inside_circumcircle(
                    points[edge.a],
                    points[edge.b],
                    points[edge.opposite],
                    points[opposite])) {
                error = "mesh contains a locally illegal edge";
                return false;
            }
            return true;
        };

    for (const Triangle& triangle : triangles) {
        if (triangle.i0 >= points.size() ||
            triangle.i1 >= points.size() ||
            triangle.i2 >= points.size()) {
            error = "triangle index is out of range";
            return false;
        }
        if (triangle.i0 == triangle.i1 ||
            triangle.i1 == triangle.i2 ||
            triangle.i2 == triangle.i0) {
            error = "triangle repeats a vertex";
            return false;
        }
        if (orient(
                points[triangle.i0],
                points[triangle.i1],
                points[triangle.i2]) == 0) {
            error = "triangle is degenerate";
            return false;
        }
        used[triangle.i0] = true;
        used[triangle.i1] = true;
        used[triangle.i2] = true;
        canonical.push_back(canonical_triangle(triangle));
        if (!add_edge(triangle.i0, triangle.i1, triangle.i2) ||
            !add_edge(triangle.i1, triangle.i2, triangle.i0) ||
            !add_edge(triangle.i2, triangle.i0, triangle.i1)) {
            return false;
        }
    }

    std::sort(
        canonical.begin(),
        canonical.end(),
        [](const Triangle& lhs, const Triangle& rhs) {
            if (lhs.i0 != rhs.i0) {
                return lhs.i0 < rhs.i0;
            }
            if (lhs.i1 != rhs.i1) {
                return lhs.i1 < rhs.i1;
            }
            return lhs.i2 < rhs.i2;
        });
    for (std::size_t i = 1; i < canonical.size(); ++i) {
        const Triangle& previous = canonical[i - 1];
        const Triangle& current = canonical[i];
        if (previous.i0 == current.i0 &&
            previous.i1 == current.i1 &&
            previous.i2 == current.i2) {
            error = "mesh contains a duplicate triangle";
            return false;
        }
    }
    if (std::find(used.begin(), used.end(), false) != used.end()) {
        error = "mesh does not reference every input point";
        return false;
    }
    return true;
}

inline bool validate_against_reference(
    const std::vector<Point>& points,
    const std::vector<Triangle>& candidate,
    const std::vector<Triangle>& reference,
    std::string& error) {
    if (meshes_equal(candidate, reference)) {
        return true;
    }
    if (candidate.size() != reference.size()) {
        std::ostringstream message;
        message << "triangle count differs: candidate=" << candidate.size()
                << " reference=" << reference.size();
        error = message.str();
        return false;
    }
    if (!validate_mesh(points, reference, error)) {
        error = "Delaunator reference failed validation: " + error;
        return false;
    }
    return validate_mesh(points, candidate, error);
}

}  // namespace delaunay32::benchmark_support
