// SPDX-License-Identifier: MIT

#include "delaunay32/delaunay.hpp"

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace delaunay32 {
namespace {

std::uint32_t triangle_vertex(
    const Triangle& triangle,
    std::size_t local_index) {
    switch (local_index) {
        case 0:
            return triangle.i0;
        case 1:
            return triangle.i1;
        default:
            return triangle.i2;
    }
}

std::uint64_t undirected_edge_key(
    std::uint32_t a,
    std::uint32_t b) {
    if (b < a) {
        std::swap(a, b);
    }
    return (static_cast<std::uint64_t>(a) << 32U) | b;
}

std::vector<std::int64_t> build_halfedges(
    const std::vector<Triangle>& triangles) {
    if (triangles.size() >
            std::numeric_limits<std::size_t>::max() / 3 ||
        static_cast<std::uint64_t>(triangles.size()) >
            static_cast<std::uint64_t>(
                std::numeric_limits<std::int64_t>::max() / 3)) {
        throw std::length_error("halfedge index range exceeded");
    }

    const std::size_t edge_count = triangles.size() * 3;
    std::vector<std::int64_t> halfedges(edge_count, -1);
    std::unordered_map<std::uint64_t, std::size_t> first_edge;
    first_edge.reserve(triangles.size() * 2);
    for (std::size_t edge = 0; edge < edge_count; ++edge) {
        const std::size_t triangle_index = edge / 3;
        const std::size_t local_index = edge % 3;
        const Triangle& triangle = triangles[triangle_index];
        const std::uint32_t a = triangle_vertex(triangle, local_index);
        const std::uint32_t b =
            triangle_vertex(triangle, (local_index + 1) % 3);
        const auto [iterator, inserted] = first_edge.emplace(
            undirected_edge_key(a, b), edge);
        if (inserted) {
            continue;
        }

        const std::size_t opposite = iterator->second;
        if (halfedges[opposite] != -1) {
            throw std::logic_error(
                "triangulation contains a non-manifold edge");
        }
        const Triangle& opposite_triangle = triangles[opposite / 3];
        const std::size_t opposite_local = opposite % 3;
        const std::uint32_t opposite_a =
            triangle_vertex(opposite_triangle, opposite_local);
        const std::uint32_t opposite_b = triangle_vertex(
            opposite_triangle,
            (opposite_local + 1) % 3);
        if (a != opposite_b || b != opposite_a) {
            throw std::logic_error(
                "adjacent triangle edges have inconsistent winding");
        }
        halfedges[edge] = static_cast<std::int64_t>(opposite);
        halfedges[opposite] = static_cast<std::int64_t>(edge);
    }
    return halfedges;
}

std::vector<std::uint32_t> build_hull(
    const std::vector<Triangle>& triangles,
    const std::vector<std::int64_t>& halfedges) {
    std::unordered_map<std::uint32_t, std::uint32_t> next_vertex;
    std::unordered_set<std::uint32_t> previous_vertex;
    std::uint32_t start = std::numeric_limits<std::uint32_t>::max();
    for (std::size_t edge = 0; edge < halfedges.size(); ++edge) {
        if (halfedges[edge] != -1) {
            continue;
        }
        const Triangle& triangle = triangles[edge / 3];
        const std::size_t local_index = edge % 3;
        const std::uint32_t a = triangle_vertex(triangle, local_index);
        const std::uint32_t b =
            triangle_vertex(triangle, (local_index + 1) % 3);
        if (!next_vertex.emplace(a, b).second ||
            !previous_vertex.insert(b).second) {
            throw std::logic_error(
                "convex hull contains a branching boundary");
        }
        start = std::min(start, a);
    }

    std::vector<std::uint32_t> hull;
    hull.reserve(next_vertex.size());
    std::uint32_t current = start;
    for (std::size_t i = 0; i < next_vertex.size(); ++i) {
        if (i != 0 && current == start) {
            throw std::logic_error(
                "convex hull contains multiple boundary loops");
        }
        hull.push_back(current);
        const auto iterator = next_vertex.find(current);
        if (iterator == next_vertex.end()) {
            throw std::logic_error("convex hull boundary is open");
        }
        current = iterator->second;
    }
    if (!hull.empty() && current != start) {
        throw std::logic_error("convex hull contains multiple boundary loops");
    }
    return hull;
}

}  // namespace

TriangulationResult Triangulator::make_result(
    const QuantizationReport& quantization,
    PredicateWidth predicate_width,
    std::vector<std::uint32_t>&& representatives) {
    TriangulationResult result;
    result.triangles = std::move(triangles_out_);
    result.quantization = quantization;
    result.predicate_width = predicate_width;
    result.actual_thread_count = active_thread_count_;
    result.representatives = std::move(representatives);

    result.halfedges = build_halfedges(result.triangles);
    if (result.triangles.empty()) {
        const auto endpoints = std::minmax_element(
            points_.begin(), points_.end(), SiteLessXY{});
        result.hull = {
            endpoints.first->original,
            endpoints.second->original,
        };
    } else {
        result.hull = build_hull(result.triangles, result.halfedges);
    }
    return result;
}

}  // namespace delaunay32
