// SPDX-License-Identifier: MIT

#include "delaunator_adapter.hpp"
#include "support.hpp"

#include "delaunay32/delaunay.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace delaunay32 {
namespace {

using benchmark_support::Dataset;

void require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

template <typename Operation>
void require_invalid(Operation&& operation, const char* label) {
    try {
        operation();
    } catch (const std::invalid_argument&) {
        return;
    }
    throw std::runtime_error(
        std::string("expected invalid_argument for ") + label);
}

void require_valid_mesh(
    const std::vector<Point>& points,
    const std::vector<Triangle>& triangles,
    const char* label) {
    std::string error;
    require(
        benchmark_support::validate_mesh(points, triangles, error),
        std::string(label) + ": " + error);
    for (const Triangle& triangle : triangles) {
        require(
            benchmark_support::orient(
                points[triangle.i0],
                points[triangle.i1],
                points[triangle.i2]) > 0,
            std::string(label) + ": output triangle is not counterclockwise");
    }
}

void require_reference_match(
    const std::vector<Point>& points,
    const std::vector<Triangle>& candidate,
    const std::vector<Triangle>& reference,
    const char* label) {
    require_valid_mesh(points, candidate, label);
    std::string error;
    require(
        benchmark_support::validate_against_reference(
            points, candidate, reference, error),
        std::string(label) + ": " + error);
}

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

bool mesh_has_edge(
    const std::vector<Triangle>& triangles,
    std::uint32_t a,
    std::uint32_t b) {
    for (const Triangle& triangle : triangles) {
        const std::array<std::uint32_t, 3> vertices = {
            triangle.i0,
            triangle.i1,
            triangle.i2,
        };
        for (std::size_t i = 0; i < vertices.size(); ++i) {
            const std::uint32_t u = vertices[i];
            const std::uint32_t v = vertices[(i + 1) % vertices.size()];
            if ((u == a && v == b) || (u == b && v == a)) {
                return true;
            }
        }
    }
    return false;
}

bool point_on_segment(
    const std::vector<Point>& points,
    std::uint32_t point,
    std::uint32_t a,
    std::uint32_t b) {
    if (benchmark_support::orient(
            points[a], points[b], points[point]) != 0) {
        return false;
    }
    return points[point].x >= std::min(points[a].x, points[b].x) &&
           points[point].x <= std::max(points[a].x, points[b].x) &&
           points[point].y >= std::min(points[a].y, points[b].y) &&
           points[point].y <= std::max(points[a].y, points[b].y);
}

bool triangle_centroid_in_ring(
    const std::vector<Point>& points,
    const Triangle& triangle,
    const std::vector<std::uint32_t>& input_ring) {
    std::vector<std::uint32_t> ring = input_ring;
    if (ring.size() > 1 && ring.front() == ring.back()) {
        ring.pop_back();
    }
    const std::int64_t x3 =
        static_cast<std::int64_t>(points[triangle.i0].x) +
        points[triangle.i1].x + points[triangle.i2].x;
    const std::int64_t y3 =
        static_cast<std::int64_t>(points[triangle.i0].y) +
        points[triangle.i1].y + points[triangle.i2].y;
    bool inside = false;
    for (std::size_t i = 0; i < ring.size(); ++i) {
        const Point& a = points[ring[i]];
        const Point& b = points[ring[(i + 1) % ring.size()]];
        const std::int64_t ay3 = static_cast<std::int64_t>(a.y) * 3;
        const std::int64_t by3 = static_cast<std::int64_t>(b.y) * 3;
        const std::int64_t orientation =
            (static_cast<std::int64_t>(b.x) - a.x) *
                (y3 - ay3) -
            (static_cast<std::int64_t>(b.y) - a.y) *
                (x3 - static_cast<std::int64_t>(a.x) * 3);
        if ((ay3 <= y3 && y3 < by3 && orientation > 0) ||
            (by3 <= y3 && y3 < ay3 && orientation < 0)) {
            inside = !inside;
        }
    }
    return inside;
}

std::int64_t twice_ring_area(
    const std::vector<Point>& points,
    const std::vector<std::uint32_t>& input_ring) {
    std::vector<std::uint32_t> ring = input_ring;
    if (ring.size() > 1 && ring.front() == ring.back()) {
        ring.pop_back();
    }
    std::int64_t area = 0;
    for (std::size_t i = 0; i < ring.size(); ++i) {
        const Point& a = points[ring[i]];
        const Point& b = points[ring[(i + 1) % ring.size()]];
        area += static_cast<std::int64_t>(a.x) * b.y -
                static_cast<std::int64_t>(a.y) * b.x;
    }
    return area < 0 ? -area : area;
}

void require_polygon_boundary_chain(
    const std::vector<Point>& points,
    const std::unordered_map<std::uint64_t, unsigned>& edges,
    std::uint32_t start,
    std::uint32_t finish,
    const std::string& label) {
    std::vector<std::vector<std::uint32_t>> adjacency(points.size());
    for (const auto& entry : edges) {
        const std::uint32_t a =
            static_cast<std::uint32_t>(entry.first >> 32U);
        const std::uint32_t b =
            static_cast<std::uint32_t>(entry.first);
        if (point_on_segment(points, a, start, finish) &&
            point_on_segment(points, b, start, finish)) {
            adjacency[a].push_back(b);
            adjacency[b].push_back(a);
        }
    }
    std::vector<std::uint32_t> pending = {start};
    std::vector<bool> reached(points.size(), false);
    reached[start] = true;
    while (!pending.empty()) {
        const std::uint32_t vertex = pending.back();
        pending.pop_back();
        for (const std::uint32_t neighbor : adjacency[vertex]) {
            if (!reached[neighbor]) {
                reached[neighbor] = true;
                pending.push_back(neighbor);
            }
        }
    }
    require(reached[finish], label + ": polygon boundary chain is missing");
}

void require_valid_polygon_mesh(
    const std::vector<Point>& points,
    const std::vector<Triangle>& triangles,
    const std::vector<std::uint32_t>& outer,
    const std::vector<std::vector<std::uint32_t>>& holes,
    const std::string& label) {
    std::unordered_map<std::uint64_t, unsigned> edges;
    edges.reserve(triangles.size() * 2);
    std::int64_t triangle_area = 0;
    for (const Triangle& triangle : triangles) {
        require(
            triangle.i0 < points.size() &&
                triangle.i1 < points.size() &&
                triangle.i2 < points.size(),
            label + ": triangle index is outside the input");
        const std::int64_t area = static_cast<std::int64_t>(
            benchmark_support::orient(
                points[triangle.i0],
                points[triangle.i1],
                points[triangle.i2]));
        require(area > 0, label + ": triangle is not counterclockwise");
        triangle_area += area;
        require(
            triangle_centroid_in_ring(points, triangle, outer),
            label + ": triangle lies outside the outer ring");
        for (const std::vector<std::uint32_t>& hole : holes) {
            require(
                !triangle_centroid_in_ring(points, triangle, hole),
                label + ": triangle lies inside a hole");
        }
        const std::array<std::uint32_t, 3> vertices = {
            triangle.i0, triangle.i1, triangle.i2};
        for (std::size_t i = 0; i < vertices.size(); ++i) {
            const std::uint64_t key = benchmark_support::edge_key(
                vertices[i], vertices[(i + 1) % vertices.size()]);
            const unsigned count = ++edges[key];
            require(count <= 2, label + ": mesh has a non-manifold edge");
        }
    }

    std::int64_t expected_area = twice_ring_area(points, outer);
    for (const std::vector<std::uint32_t>& hole : holes) {
        expected_area -= twice_ring_area(points, hole);
    }
    require(
        triangle_area == expected_area,
        label + ": triangle area does not match the polygon domain");

    const auto require_ring = [&](const std::vector<std::uint32_t>& input) {
        std::vector<std::uint32_t> ring = input;
        if (ring.size() > 1 && ring.front() == ring.back()) {
            ring.pop_back();
        }
        for (std::size_t i = 0; i < ring.size(); ++i) {
            require_polygon_boundary_chain(
                points,
                edges,
                ring[i],
                ring[(i + 1) % ring.size()],
                label);
        }
    };
    require_ring(outer);
    for (const std::vector<std::uint32_t>& hole : holes) {
        require_ring(hole);
    }
}

bool segments_properly_intersect(
    const Point& a,
    const Point& b,
    const Point& c,
    const Point& d) {
    const auto ab_c = benchmark_support::orient(a, b, c);
    const auto ab_d = benchmark_support::orient(a, b, d);
    const auto cd_a = benchmark_support::orient(c, d, a);
    const auto cd_b = benchmark_support::orient(c, d, b);
    const auto opposite = [](const auto& lhs, const auto& rhs) {
        return (lhs > 0 && rhs < 0) || (lhs < 0 && rhs > 0);
    };
    return opposite(ab_c, ab_d) && opposite(cd_a, cd_b);
}

void require_valid_constrained_mesh(
    const std::vector<Point>& points,
    const std::vector<Triangle>& triangles,
    const std::vector<Constraint>& constraints,
    const std::string& label) {
    struct Record {
        std::uint32_t a = 0;
        std::uint32_t b = 0;
        std::uint32_t opposite = 0;
        std::uint32_t neighbor_opposite = 0;
        unsigned count = 0;
    };
    std::unordered_map<std::uint64_t, Record> edges;
    edges.reserve(triangles.size() * 2);
    std::vector<bool> used(points.size(), false);
    const auto add_edge = [&](std::uint32_t a,
                              std::uint32_t b,
                              std::uint32_t opposite) {
        const std::uint64_t key = benchmark_support::edge_key(a, b);
        auto [iterator, inserted] =
            edges.emplace(key, Record{a, b, opposite, 0, 1});
        if (!inserted) {
            require(
                iterator->second.count == 1,
                label + ": mesh contains a non-manifold edge");
            iterator->second.neighbor_opposite = opposite;
            iterator->second.count = 2;
        }
    };
    for (const Triangle& triangle : triangles) {
        require(
            triangle.i0 < points.size() &&
                triangle.i1 < points.size() &&
                triangle.i2 < points.size(),
            label + ": triangle index is outside the input");
        require(
            benchmark_support::orient(
                points[triangle.i0],
                points[triangle.i1],
                points[triangle.i2]) > 0,
            label + ": triangle is not counterclockwise");
        used[triangle.i0] = true;
        used[triangle.i1] = true;
        used[triangle.i2] = true;
        add_edge(triangle.i0, triangle.i1, triangle.i2);
        add_edge(triangle.i1, triangle.i2, triangle.i0);
        add_edge(triangle.i2, triangle.i0, triangle.i1);
    }
    require(
        std::find(used.begin(), used.end(), false) == used.end(),
        label + ": mesh does not reference every input point");

    std::unordered_set<std::uint64_t> constrained_edges;
    for (const Constraint constraint : constraints) {
        require(
            constraint.i0 < points.size() && constraint.i1 < points.size(),
            label + ": test constraint endpoint is invalid");
        std::vector<std::vector<std::uint32_t>> adjacency(points.size());
        for (const auto& entry : edges) {
            const Record& edge = entry.second;
            if (point_on_segment(
                    points, edge.a, constraint.i0, constraint.i1) &&
                point_on_segment(
                    points, edge.b, constraint.i0, constraint.i1)) {
                adjacency[edge.a].push_back(edge.b);
                adjacency[edge.b].push_back(edge.a);
                constrained_edges.insert(entry.first);
            } else {
                require(
                    !segments_properly_intersect(
                        points[constraint.i0],
                        points[constraint.i1],
                        points[edge.a],
                        points[edge.b]),
                    label + ": mesh edge crosses a constraint");
            }
        }
        std::vector<std::uint32_t> pending = {constraint.i0};
        std::vector<bool> reached(points.size(), false);
        reached[constraint.i0] = true;
        while (!pending.empty()) {
            const std::uint32_t vertex = pending.back();
            pending.pop_back();
            for (const std::uint32_t neighbor : adjacency[vertex]) {
                if (!reached[neighbor]) {
                    reached[neighbor] = true;
                    pending.push_back(neighbor);
                }
            }
        }
        require(
            reached[constraint.i1],
            label + ": constraint is not represented by a mesh-edge chain");
    }

    for (const auto& entry : edges) {
        const Record& edge = entry.second;
        if (edge.count != 2 ||
            constrained_edges.count(entry.first) != 0) {
            continue;
        }
        require(
            !benchmark_support::inside_circumcircle(
                points[edge.a],
                points[edge.b],
                points[edge.opposite],
                points[edge.neighbor_opposite]),
            label + ": unconstrained edge is locally illegal");
    }
}

void require_full_topology(
    const std::vector<Point>& points,
    const TriangulationResult& result,
    const std::string& label) {
    const std::size_t edge_count = result.triangles.size() * 3;
    require(
        result.halfedges.size() == edge_count,
        label + ": halfedge count does not match triangles");

    std::unordered_map<std::uint32_t, std::uint32_t> boundary;
    for (std::size_t edge = 0; edge < edge_count; ++edge) {
        const Triangle& triangle = result.triangles[edge / 3];
        const std::size_t local = edge % 3;
        const std::uint32_t a = triangle_vertex(triangle, local);
        const std::uint32_t b =
            triangle_vertex(triangle, (local + 1) % 3);
        const std::int64_t opposite = result.halfedges[edge];
        if (opposite == -1) {
            require(
                boundary.emplace(a, b).second,
                label + ": boundary branches at a hull vertex");
            continue;
        }
        require(
            opposite >= 0 &&
                static_cast<std::uint64_t>(opposite) < edge_count,
            label + ": halfedge points outside the flattened edge array");
        const std::size_t opposite_edge =
            static_cast<std::size_t>(opposite);
        require(
            result.halfedges[opposite_edge] ==
                static_cast<std::int64_t>(edge),
            label + ": halfedge pairing is not reciprocal");
        const Triangle& neighbor = result.triangles[opposite_edge / 3];
        const std::size_t neighbor_local = opposite_edge % 3;
        require(
            a == triangle_vertex(
                     neighbor, (neighbor_local + 1) % 3) &&
                b == triangle_vertex(neighbor, neighbor_local),
            label + ": paired halfedges do not reverse endpoints");
    }

    require(
        result.hull.size() == boundary.size(),
        label + ": hull length differs from boundary edge count");
    require(!result.hull.empty(), label + ": nonempty mesh has no hull");
    require(
        result.hull.front() ==
            *std::min_element(result.hull.begin(), result.hull.end()),
        label + ": hull is not rotated to its lowest input index");
    long double twice_area = 0.0L;
    for (std::size_t i = 0; i < result.hull.size(); ++i) {
        const std::uint32_t a = result.hull[i];
        const std::uint32_t b = result.hull[(i + 1) % result.hull.size()];
        require(
            a < points.size() && b < points.size(),
            label + ": hull index is outside the input");
        const auto iterator = boundary.find(a);
        require(
            iterator != boundary.end() && iterator->second == b,
            label + ": hull order does not follow boundary halfedges");
        twice_area +=
            static_cast<long double>(points[a].x) * points[b].y -
            static_cast<long double>(points[a].y) * points[b].x;
    }
    require(twice_area > 0.0L, label + ": hull is not counterclockwise");
}

std::vector<Point> quantize_from_report(
    const std::vector<FloatPoint>& points,
    const QuantizationReport& report) {
    std::vector<Point> quantized;
    quantized.reserve(points.size());
    for (const FloatPoint& point : points) {
        quantized.push_back({
            static_cast<std::int32_t>(std::llround(
                (static_cast<double>(point.x) - report.origin_x) *
                report.scale)),
            static_cast<std::int32_t>(std::llround(
                (static_cast<double>(point.y) - report.origin_y) *
                report.scale)),
        });
    }
    return quantized;
}

#if defined(__SIZEOF_INT128__)
std::pair<std::uint64_t, std::uint64_t> coordinate_spans(
    const std::vector<Point>& points) {
    std::int32_t min_x = points.front().x;
    std::int32_t max_x = min_x;
    std::int32_t min_y = points.front().y;
    std::int32_t max_y = min_y;
    for (const Point& point : points) {
        min_x = std::min(min_x, point.x);
        max_x = std::max(max_x, point.x);
        min_y = std::min(min_y, point.y);
        max_y = std::max(max_y, point.y);
    }
    return {
        static_cast<std::uint64_t>(
            static_cast<std::int64_t>(max_x) - min_x),
        static_cast<std::uint64_t>(
            static_cast<std::int64_t>(max_y) - min_y),
    };
}
#endif

void test_predicate_selection() {
    const auto expect =
        [](std::uint64_t x_span,
           std::uint64_t y_span,
           PredicateWidth expected,
           const char* label) {
            require(
                Triangulator::predicate_width_for_spans(x_span, y_span) ==
                    expected,
                std::string("predicate-width selector failed for ") + label);
        };

    expect(29609, 29609, PredicateWidth::Int64, "int64 square boundary");
    expect(50000, 10000, PredicateWidth::Int64, "thin int64 domain");
#if defined(__SIZEOF_INT128__)
    expect(29610, 29610, PredicateWidth::Int128, "first wider square");
    expect(100000, 100000, PredicateWidth::Int128, "100k square");
    expect(
        1940470527,
        1940470527,
        PredicateWidth::Int128,
        "int128 square boundary");
    expect(
        1940470528,
        1940470528,
        PredicateWidth::Unsupported,
        "first unsupported square");
    expect(
        std::numeric_limits<std::uint32_t>::max(),
        std::numeric_limits<std::uint32_t>::max(),
        PredicateWidth::Unsupported,
        "full int32 square");
    require(
        Triangulator::int64_wide_intermediates_for_spans(100000, 100000),
        "100k mixed-width selector failed");
    require(
        !Triangulator::int64_wide_intermediates_for_spans(
            1940470527, 1940470527),
        "large int128 selector unexpectedly chose int64 intermediates");
#else
    expect(29610, 29610, PredicateWidth::Unsupported, "unsupported square");
#endif
}

void test_deterministic_cases() {
    constexpr std::array<std::size_t, 16> sizes = {
        3,
        4,
        5,
        15,
        16,
        17,
        127,
        998,
        999,
        1000,
        1001,
        9998,
        9999,
        10000,
        10001,
        50000,
    };
    constexpr std::array<Dataset, 3> datasets = {
        Dataset::Uniform,
        Dataset::Clustered,
        Dataset::Diagonal,
    };

    Triangulator serial(1);
    Triangulator parallel(2);
    DelaunatorBaseline delaunator;
    for (std::size_t i = 0; i < sizes.size(); ++i) {
        const Dataset dataset = datasets[i % datasets.size()];
        const std::vector<Point> points =
            benchmark_support::generate_points(
                dataset, sizes[i], 0x1000ULL + i, 20000);
        const std::vector<Triangle> reference =
            delaunator.triangulate(points);
        const std::vector<Triangle> candidate =
            serial.triangulate_int(points);
        require_reference_match(
            points, candidate, reference, "serial deterministic case");

        if (sizes[i] >= 50000) {
            const std::vector<Triangle> parallel_candidate =
                parallel.triangulate_int(points);
            require_reference_match(
                points,
                parallel_candidate,
                reference,
                "parallel deterministic case");
            require(
                benchmark_support::meshes_equal(
                    candidate, parallel_candidate),
                "serial and parallel triangle sets differ");
        }
    }
}

void test_full_integer_result() {
    const std::vector<Point> points = {
        {0, 0},
        {10, 0},
        {0, 0},
        {10, 10},
        {0, 10},
        {5, 5},
        {5, 5},
    };
    const std::vector<std::uint32_t> expected_representatives = {
        0, 1, 0, 3, 4, 5, 5,
    };

    Triangulator triangulator(4);
    const TriangulationResult result =
        triangulator.triangulate_int_full(points);
    require(!result.triangles.empty(), "full integer result is empty");
    require(
        result.representatives == expected_representatives,
        "full integer representative mapping is incorrect");
    require(
        result.predicate_width == PredicateWidth::Int64,
        "full integer result reported the wrong predicate width");
    require(
        result.actual_thread_count == 1,
        "small full integer result did not report serial execution");
    require(
        result.quantization.scale == 0.0 &&
            result.quantization.unique_points == 0,
        "integer result unexpectedly contains quantization metadata");
    require_full_topology(points, result, "full integer result");

    const std::vector<Triangle> triangle_only =
        triangulator.triangulate_int(points);
    require(
        benchmark_support::meshes_equal(
            result.triangles, triangle_only),
        "full and triangle-only integer APIs differ");

    const std::vector<Point> collinear = {
        {3, 7}, {-4, 7}, {9, 7}, {0, 7},
    };
    const TriangulationResult collinear_result =
        triangulator.triangulate_int_full(collinear);
    require(
        collinear_result.triangles.empty() &&
            collinear_result.halfedges.empty(),
        "collinear full result contains faces or halfedges");
    require(
        collinear_result.hull == std::vector<std::uint32_t>({1, 2}),
        "collinear full result does not contain its endpoints");
}

void test_float_api() {
    const std::vector<FloatPoint> points = {
        {0.125F, 0.25F},
        {5.5F, 0.1F},
        {6.0F, 4.5F},
        {-1.0F, 5.0F},
        {2.75F, 2.0F},
        {4.125F, 3.25F},
    };
    const std::vector<FloatPoint> original = points;

    Triangulator triangulator;
    const TriangulationResult result =
        triangulator.triangulate_float_full(points);
    const QuantizationReport& report = result.quantization;
    const std::vector<Triangle>& vector_mesh = result.triangles;
    require(!vector_mesh.empty(), "float input produced no triangles");
    require(report.scale > 0.0, "float quantization scale is not positive");
    require(report.grid_step > 0.0, "float grid step is not positive");
    require(
        std::abs(report.scale * report.grid_step - 1.0) < 1e-12,
        "float grid step is not the inverse scale");
    require(
        report.max_coordinate_error <= report.grid_step * 0.500001,
        "reported float quantization error exceeds half a grid step");
    require(
        report.unique_points == points.size() &&
            report.collapsed_points == 0,
        "distinct float input unexpectedly collapsed");
    for (std::size_t i = 0; i < points.size(); ++i) {
        require(
            points[i].x == original[i].x && points[i].y == original[i].y,
            "float input coordinates were modified");
    }

    const std::vector<Point> quantized =
        quantize_from_report(points, report);
    const std::vector<Triangle> integer_mesh =
        triangulator.triangulate_int(quantized);
    require(
        benchmark_support::meshes_equal(vector_mesh, integer_mesh),
        "float API differs from its reported integer quantization");

    const std::vector<Triangle> triangle_only_mesh =
        triangulator.triangulate_float(points);
    require(
        benchmark_support::meshes_equal(vector_mesh, triangle_only_mesh),
        "full and triangle-only float APIs differ");
}

void test_quantization_options() {
    const std::vector<FloatPoint> points = {
        {0.125F, 0.125F},
        {2.125F, 0.125F},
        {2.125F, 2.125F},
        {0.125F, 2.125F},
        {1.125F, 0.875F},
    };
    QuantizationOptions grid_options;
    grid_options.mode = QuantizationMode::GridStep;
    grid_options.grid_step = 0.25;

    Triangulator triangulator;
    const TriangulationResult grid_result =
        triangulator.triangulate_float_full(points, grid_options);
    require(
        grid_result.quantization.origin_x == 0.125 &&
            grid_result.quantization.origin_y == 0.125 &&
            grid_result.quantization.scale == 4.0 &&
            grid_result.quantization.grid_step == 0.25,
        "grid-step quantization report is incorrect");
    const std::vector<Point> grid_points =
        quantize_from_report(points, grid_result.quantization);
    require(
        benchmark_support::meshes_equal(
            grid_result.triangles,
            triangulator.triangulate_int(grid_points)),
        "grid-step float mesh differs from its integer mapping");

    QuantizationOptions fixed_options;
    fixed_options.mode = QuantizationMode::FixedScale;
    fixed_options.origin_x = -10.0;
    fixed_options.origin_y = -20.0;
    fixed_options.scale = 8.0;
    const TriangulationResult fixed_result =
        triangulator.triangulate_float_full(points, fixed_options);
    require(
        fixed_result.quantization.origin_x == -10.0 &&
            fixed_result.quantization.origin_y == -20.0 &&
            fixed_result.quantization.scale == 8.0,
        "fixed quantization report is incorrect");
    require(
        benchmark_support::meshes_equal(
            fixed_result.triangles,
            triangulator.triangulate_int(
                quantize_from_report(
                    points, fixed_result.quantization))),
        "fixed-scale float mesh differs from its integer mapping");

    QuantizationOptions strict_error = grid_options;
    strict_error.max_coordinate_error = 0.05;
    const std::vector<FloatPoint> inexact = {
        {0.0F, 0.0F},
        {2.0F, 0.0F},
        {0.0F, 2.0F},
        {0.12F, 0.12F},
    };
    require_invalid(
        [&] { triangulator.triangulate_float(inexact, strict_error); },
        "requested maximum quantization error");

    QuantizationOptions reject_collisions;
    reject_collisions.mode = QuantizationMode::GridStep;
    reject_collisions.grid_step = 1.0;
    reject_collisions.collision_policy =
        QuantizationCollisionPolicy::Reject;
    const std::vector<FloatPoint> colliding = {
        {0.0F, 0.0F},
        {0.2F, 0.2F},
        {2.0F, 0.0F},
        {0.0F, 2.0F},
    };
    require_invalid(
        [&] {
            triangulator.triangulate_float(
                colliding, reject_collisions);
        },
        "rejected quantization collision");
}

void test_float_quantization_collisions() {
    const float maximum = std::numeric_limits<float>::max();
    const std::vector<FloatPoint> points = {
        {-maximum, -maximum},
        {maximum, -maximum},
        {maximum, maximum},
        {-maximum, maximum},
        {0.0F, 0.0F},
        {1.0F, 1.0F},
    };

    Triangulator triangulator;
    const TriangulationResult result =
        triangulator.triangulate_float_full(points);
    const QuantizationReport& report = result.quantization;
    const std::vector<Triangle>& candidate = result.triangles;
    require(
        report.unique_points == 5 && report.collapsed_points == 1,
        "full-range float collision was not reported");
    for (const Triangle& triangle : candidate) {
        require(
            triangle.i0 != 5 && triangle.i1 != 5 && triangle.i2 != 5,
            "quantized duplicate did not retain its lowest input index");
    }

    const std::vector<Point> quantized =
        quantize_from_report(points, report);
    const std::vector<Triangle> expected =
        triangulator.triangulate_int(quantized);
    require(
        benchmark_support::meshes_equal(candidate, expected),
        "full-range float mesh differs from integer quantization");
}

void test_full_float_result() {
    const float maximum = std::numeric_limits<float>::max();
    const std::vector<FloatPoint> points = {
        {-maximum, -maximum},
        {maximum, -maximum},
        {maximum, maximum},
        {-maximum, maximum},
        {0.0F, 0.0F},
        {1.0F, 1.0F},
    };
    const std::vector<std::uint32_t> expected_representatives = {
        0, 1, 2, 3, 4, 4,
    };

    Triangulator triangulator;
    const TriangulationResult result =
        triangulator.triangulate_float_full(points);
    require(
        result.quantization.unique_points == 5 &&
            result.quantization.collapsed_points == 1,
        "full float result has incorrect quantization counts");
    require(
        result.representatives == expected_representatives,
        "full float representative mapping is incorrect");
    require(
        result.predicate_width != PredicateWidth::Unsupported,
        "full float result reported unsupported predicates");
    const std::vector<Point> quantized =
        quantize_from_report(points, result.quantization);
    require_full_topology(quantized, result, "full float result");

    const std::vector<Triangle> triangle_only =
        triangulator.triangulate_float(points);
    require(
        benchmark_support::meshes_equal(
            result.triangles, triangle_only),
        "full and triangle-only float APIs differ");
}

void test_float_parallel() {
    const std::vector<Point> integer_points =
        benchmark_support::generate_points(
            Dataset::Uniform, 60000, 0xf10a7ULL, 20000);
    std::vector<FloatPoint> points;
    points.reserve(integer_points.size());
    for (const Point& point : integer_points) {
        points.push_back({
            static_cast<float>(point.x) * 0.25F + 0.125F,
            static_cast<float>(point.y) * 0.25F - 0.375F,
        });
    }

    Triangulator serial(1);
    Triangulator parallel(2);
    const std::vector<Triangle> serial_mesh =
        serial.triangulate_float(points);
    const std::vector<Triangle> parallel_mesh =
        parallel.triangulate_float(points);
    require(
        benchmark_support::meshes_equal(serial_mesh, parallel_mesh),
        "serial and parallel float triangle sets differ");

    const TriangulationResult full =
        parallel.triangulate_float_full(points);
    require(
        benchmark_support::meshes_equal(serial_mesh, full.triangles),
        "parallel full float triangle set differs");
    require(
        full.quantization.unique_points == points.size() &&
            full.quantization.collapsed_points == 0,
        "parallel float quantization unexpectedly collapsed points");
    require(
        full.actual_thread_count == 2,
        "parallel full result reported the wrong actual thread count");
    require(
        full.representatives.size() == points.size(),
        "parallel full result omitted representative entries");
    for (std::size_t i = 0; i < full.representatives.size(); ++i) {
        require(
            full.representatives[i] == i,
            "parallel full result changed a unique representative");
    }
    require_full_topology(
        quantize_from_report(points, full.quantization),
        full,
        "parallel full float result");
}

void test_float_collinear_input() {
    Triangulator triangulator;
    const std::vector<FloatPoint> points = {
        {-4.5F, 7.25F},
        {-1.5F, 7.25F},
        {0.0F, 7.25F},
        {3.0F, 7.25F},
        {9.0F, 7.25F},
    };
    require(
        triangulator.triangulate_float(points).empty(),
        "axis-aligned collinear float input produced a triangle");
}

void test_duplicates() {
    const std::vector<Point> unique = {
        {0, 0},
        {10, 0},
        {0, 10},
        {10, 10},
        {4, 3},
    };
    const std::array<std::uint32_t, 5> representative = {0, 1, 3, 4, 6};
    const std::vector<Point> duplicated = {
        unique[0],
        unique[1],
        unique[0],
        unique[2],
        unique[3],
        unique[1],
        unique[4],
        unique[4],
    };

    Triangulator triangulator;
    std::vector<Triangle> expected = triangulator.triangulate_int(unique);
    for (Triangle& triangle : expected) {
        triangle.i0 = representative[triangle.i0];
        triangle.i1 = representative[triangle.i1];
        triangle.i2 = representative[triangle.i2];
    }

    const std::vector<Triangle> candidate =
        triangulator.triangulate_int(duplicated);
    require(
        benchmark_support::meshes_equal(expected, candidate),
        "duplicate compaction changed the triangle set or representatives");
}

void test_collinear_input() {
    Triangulator triangulator;
    const std::vector<Triangle> horizontal = triangulator.triangulate_int(
        {{-4, 7}, {-1, 7}, {0, 7}, {3, 7}, {9, 7}});
    require(horizontal.empty(), "collinear input produced a triangle");

    const std::vector<Triangle> diagonal = triangulator.triangulate_int(
        {{-4, -7}, {-1, -1}, {0, 1}, {3, 7}, {9, 19}});
    require(diagonal.empty(), "diagonal collinear input produced a triangle");
}

void test_wide_predicates() {
#if defined(__SIZEOF_INT128__)
    constexpr std::array<Dataset, 3> datasets = {
        Dataset::Uniform,
        Dataset::Clustered,
        Dataset::Diagonal,
    };
    Triangulator serial(1);
    DelaunatorBaseline delaunator;
    for (std::size_t i = 0; i < datasets.size(); ++i) {
        const std::vector<Point> points =
            benchmark_support::generate_points(
                datasets[i], 10000, 0x9000ULL + i, 100000);
        const auto [x_span, y_span] = coordinate_spans(points);
        require(
            Triangulator::predicate_width_for_spans(x_span, y_span) ==
                PredicateWidth::Int128,
            "wide-domain input did not select int128 predicates");
        const std::vector<Triangle> candidate =
            serial.triangulate_int(points);
        const std::vector<Triangle> reference =
            delaunator.triangulate(points);
        require_reference_match(
            points, candidate, reference, "wide-predicate case");
    }

    const std::vector<Point> parallel_points =
        benchmark_support::generate_points(
            Dataset::Uniform, 60000, 0xfeedULL, 100000);
    const std::vector<Triangle> serial_mesh =
        serial.triangulate_int(parallel_points);
    Triangulator parallel(2);
    const std::vector<Triangle> parallel_mesh =
        parallel.triangulate_int(parallel_points);
    require(
        benchmark_support::meshes_equal(serial_mesh, parallel_mesh),
        "wide serial and parallel triangle sets differ");
    require_valid_mesh(
        parallel_points, parallel_mesh, "wide parallel case");

    const std::vector<Point> extreme_strip = {
        {std::numeric_limits<std::int32_t>::min(), 0},
        {std::numeric_limits<std::int32_t>::max(), 0},
        {std::numeric_limits<std::int32_t>::min(), 1},
        {std::numeric_limits<std::int32_t>::max(), 1},
    };
    const std::vector<Triangle> strip_mesh =
        serial.triangulate_int(extreme_strip);
    require(
        Triangulator::predicate_width_for_spans(
            std::numeric_limits<std::uint32_t>::max(), 1) ==
            PredicateWidth::Int128,
        "full-width thin domain did not select int128 predicates");
    require(
        !Triangulator::int64_wide_intermediates_for_spans(
            std::numeric_limits<std::uint32_t>::max(), 1),
        "full-width thin domain selected mixed-width intermediates");
    require_valid_mesh(extreme_strip, strip_mesh, "full-width thin domain");
#endif
}

void test_single_constraint_recovery() {
    const std::vector<Point> points = {
        {0, 0},
        {100, 0},
        {100, 100},
        {0, 100},
    };
    Triangulator triangulator;
    const std::vector<Triangle> unconstrained =
        triangulator.triangulate_int(points);
    const Constraint missing = mesh_has_edge(unconstrained, 0, 2)
                                   ? Constraint{1, 3}
                                   : Constraint{0, 2};
    const std::vector<Triangle> constrained =
        triangulator.triangulate_constrained_int(points, {missing});
    require_valid_constrained_mesh(
        points, constrained, {missing}, "single constraint recovery");
    require(
        mesh_has_edge(constrained, missing.i0, missing.i1),
        "single constraint was not recovered");
}

void test_constrained_random_cases() {
    for (std::uint64_t seed = 0; seed < 300; ++seed) {
        const std::size_t point_count = 20 + seed % 80;
        const std::vector<Point> points =
            benchmark_support::generate_points(
                Dataset::Uniform, point_count, 0xCD7000ULL + seed, 2000);
        const Constraint constraint = {
            static_cast<std::uint32_t>(seed % point_count),
            static_cast<std::uint32_t>(
                (seed * 17 + point_count / 2) % point_count),
        };
        if (constraint.i0 == constraint.i1) {
            continue;
        }
        Triangulator triangulator(seed % 2 == 0 ? 1 : 2);
        const std::vector<Triangle> triangles =
            triangulator.triangulate_constrained_int(
                points, {constraint});
        require_valid_constrained_mesh(
            points,
            triangles,
            {constraint},
            "random single-constraint case");
    }
}

void test_multiple_constraints() {
    const std::vector<Point> points =
        benchmark_support::generate_points(
            Dataset::Uniform, 80, 0xC07A1AULL, 5000);
    std::vector<Constraint> constraints;
    for (std::uint32_t vertex = 1; vertex < 16; ++vertex) {
        constraints.push_back({0, vertex});
    }
    constraints.push_back(constraints[3]);
    constraints.push_back({constraints[5].i1, constraints[5].i0});

    Triangulator triangulator;
    const std::vector<Triangle> triangles =
        triangulator.triangulate_constrained_int(points, constraints);
    require_valid_constrained_mesh(
        points, triangles, constraints, "multiple constraints");
}

void test_collinear_constraint_chain() {
    const std::vector<Point> points = {
        {0, 0},
        {10, 0},
        {20, 0},
        {30, 0},
        {0, 20},
        {30, 20},
        {15, 10},
    };
    const std::vector<Constraint> constraints = {
        {0, 3},
        {1, 2},
        {3, 0},
    };
    Triangulator triangulator;
    const std::vector<Triangle> triangles =
        triangulator.triangulate_constrained_int(points, constraints);
    require_valid_constrained_mesh(
        points, triangles, constraints, "collinear constraint chain");
    require(mesh_has_edge(triangles, 0, 1), "constraint chain misses 0-1");
    require(mesh_has_edge(triangles, 1, 2), "constraint chain misses 1-2");
    require(mesh_has_edge(triangles, 2, 3), "constraint chain misses 2-3");
}

void test_constraint_cycle() {
    constexpr std::size_t ring_size = 24;
    constexpr double pi = 3.14159265358979323846;
    std::vector<Point> points;
    points.reserve(ring_size + 100);
    for (std::size_t i = 0; i < ring_size; ++i) {
        const double angle =
            2.0 * pi * static_cast<double>(i) /
            static_cast<double>(ring_size);
        const double radius = i % 2 == 0 ? 900.0 : 450.0;
        points.push_back({
            static_cast<std::int32_t>(
                std::llround(radius * std::cos(angle))),
            static_cast<std::int32_t>(
                std::llround(radius * std::sin(angle))),
        });
    }
    std::vector<Point> background =
        benchmark_support::generate_points(
            Dataset::Uniform, 100, 0xC1C1EULL, 1800);
    for (Point& point : background) {
        point.x -= 900;
        point.y -= 900;
    }
    points.insert(points.end(), background.begin(), background.end());

    std::vector<Constraint> constraints;
    constraints.reserve(ring_size);
    for (std::uint32_t i = 0; i < ring_size; ++i) {
        constraints.push_back({
            i,
            static_cast<std::uint32_t>((i + 1) % ring_size),
        });
    }
    Triangulator triangulator;
    const std::vector<Triangle> triangles =
        triangulator.triangulate_constrained_int(points, constraints);
    require_valid_constrained_mesh(
        points, triangles, constraints, "nonconvex constraint cycle");
}

void test_constraint_split_at_unconnected_site() {
    const std::vector<Point> points = {
        {0, 0},
        {10, 0},
        {20, 0},
        {5, 1},
        {5, -1},
        {15, 1},
        {15, -1},
        {-5, -10},
        {25, -10},
        {25, 10},
        {-5, 10},
    };
    Triangulator triangulator;
    const std::vector<Triangle> ordinary =
        triangulator.triangulate_int(points);
    require(
        !mesh_has_edge(ordinary, 0, 1) &&
            !mesh_has_edge(ordinary, 1, 2),
        "split-site fixture unexpectedly starts with constraint subedges");

    const std::vector<Constraint> constraints = {{0, 2}};
    const std::vector<Triangle> constrained =
        triangulator.triangulate_constrained_int(points, constraints);
    require_valid_constrained_mesh(
        points,
        constrained,
        constraints,
        "constraint split at an initially unconnected site");
    require(mesh_has_edge(constrained, 0, 1), "split recovery misses 0-1");
    require(mesh_has_edge(constrained, 1, 2), "split recovery misses 1-2");
}

void test_constraints_meet_at_existing_site() {
    const std::vector<Point> points = {
        {-20, 0},
        {20, 0},
        {0, -20},
        {0, 20},
        {0, 0},
        {-30, -30},
        {30, -30},
        {30, 30},
        {-30, 30},
    };
    const std::vector<Constraint> constraints = {{0, 1}, {2, 3}};
    Triangulator triangulator;
    const std::vector<Triangle> triangles =
        triangulator.triangulate_constrained_int(points, constraints);
    require_valid_constrained_mesh(
        points,
        triangles,
        constraints,
        "constraints meeting at an existing site");
    require(mesh_has_edge(triangles, 0, 4), "horizontal chain misses 0-4");
    require(mesh_has_edge(triangles, 4, 1), "horizontal chain misses 4-1");
    require(mesh_has_edge(triangles, 2, 4), "vertical chain misses 2-4");
    require(mesh_has_edge(triangles, 4, 3), "vertical chain misses 4-3");
}

void test_constrained_duplicates_and_empty_input() {
    const std::vector<Point> unique = {
        {0, 0},
        {100, 0},
        {0, 100},
        {100, 100},
        {35, 45},
    };
    const std::array<std::uint32_t, 5> representative = {0, 1, 3, 4, 6};
    const std::vector<Point> duplicated = {
        unique[0],
        unique[1],
        unique[0],
        unique[2],
        unique[3],
        unique[1],
        unique[4],
        unique[4],
    };

    Triangulator triangulator;
    std::vector<Triangle> expected =
        triangulator.triangulate_constrained_int(unique, {{0, 4}});
    for (Triangle& triangle : expected) {
        triangle.i0 = representative[triangle.i0];
        triangle.i1 = representative[triangle.i1];
        triangle.i2 = representative[triangle.i2];
    }
    const std::vector<Triangle> candidate =
        triangulator.triangulate_constrained_int(duplicated, {{2, 7}});
    require(
        benchmark_support::meshes_equal(expected, candidate),
        "duplicate constraint endpoints did not map to representatives");
    require(
        mesh_has_edge(candidate, representative[0], representative[4]),
        "constraint using duplicate endpoints was not recovered");

    const std::vector<Triangle> ordinary =
        triangulator.triangulate_int(unique);
    const std::vector<Triangle> empty =
        triangulator.triangulate_constrained_int(unique, {});
    require(
        benchmark_support::meshes_equal(ordinary, empty),
        "empty constraint input changed the ordinary triangulation");

    require_invalid(
        [&] {
            triangulator.triangulate_constrained_int(
                duplicated, {{0, 2}});
        },
        "distinct constraint indices at the same coordinate");
}

void test_collinear_constrained_input() {
    const std::vector<Point> points = {
        {-20, -39},
        {-10, -19},
        {0, 1},
        {10, 21},
        {20, 41},
    };
    Triangulator triangulator;
    require(
        triangulator.triangulate_constrained_int(points, {{0, 4}}).empty(),
        "collinear constrained input produced a triangle");
}

void test_parallel_constraints() {
    std::vector<Point> grid;
    grid.reserve(50000);
    for (std::int32_t y = 0; y < 200; ++y) {
        for (std::int32_t x = 0; x < 250; ++x) {
            grid.push_back({x, y});
        }
    }
    const std::vector<Constraint> grid_constraints = {{0, 49999}};
    Triangulator grid_serial(1);
    Triangulator grid_parallel(4);
    const std::vector<Triangle> grid_serial_mesh =
        grid_serial.triangulate_constrained_int(grid, grid_constraints);
    const std::vector<Triangle> grid_parallel_mesh =
        grid_parallel.triangulate_constrained_int(grid, grid_constraints);
    require(
        benchmark_support::meshes_equal(
            grid_serial_mesh, grid_parallel_mesh),
        "cocircular constrained serial and parallel triangle sets differ");
    require_valid_constrained_mesh(
        grid,
        grid_parallel_mesh,
        grid_constraints,
        "cocircular parallel constrained case");

#if defined(__SIZEOF_INT128__)
    const std::vector<Point> points =
        benchmark_support::generate_points(
            Dataset::Uniform, 60000, 0xCD7128ULL, 100000);
    const std::vector<Constraint> constraints = {{17, 42371}};
    Triangulator serial(1);
    Triangulator parallel(4);
    const std::vector<Triangle> serial_mesh =
        serial.triangulate_constrained_int(points, constraints);
    const std::vector<Triangle> parallel_mesh =
        parallel.triangulate_constrained_int(points, constraints);
    require(
        benchmark_support::meshes_equal(serial_mesh, parallel_mesh),
        "wide constrained serial and parallel triangle sets differ");
    require_valid_constrained_mesh(
        points,
        parallel_mesh,
        constraints,
        "wide parallel constrained case");
#endif
}

void test_invalid_constraints() {
    const std::vector<Point> square = {
        {0, 0},
        {100, 0},
        {100, 100},
        {0, 100},
    };
    Triangulator triangulator;
    require_invalid(
        [&] {
            triangulator.triangulate_constrained_int(
                square, {{0, 4}});
        },
        "constraint endpoint outside point array");
    require_invalid(
        [&] {
            triangulator.triangulate_constrained_int(
                square, {{2, 2}});
        },
        "zero-length constraint");
    require_invalid(
        [&] {
            triangulator.triangulate_constrained_int(
                square, {{0, 2}, {1, 3}});
        },
        "properly intersecting constraints");

    const std::vector<Triangle> recovered =
        triangulator.triangulate_constrained_int(square, {{0, 2}});
    require_valid_constrained_mesh(
        square,
        recovered,
        {{0, 2}},
        "reuse after invalid constrained input");
}

void test_polygon_outer_domains() {
    {
        const std::vector<Point> points = {
            {0, 0}, {100, 0}, {100, 100}, {0, 100},
            {20, 20}, {80, 20}, {50, 75},
            {-40, 50}, {140, 50}, {50, -40}, {50, 140},
        };
        const std::vector<std::uint32_t> outer = {0, 3, 2, 1, 0};
        Triangulator triangulator;
        const std::vector<Triangle> triangles =
            triangulator.triangulate_polygon_int(points, outer);
        require_valid_polygon_mesh(
            points, triangles, outer, {}, "clockwise square polygon");
        std::vector<bool> used(points.size(), false);
        for (const Triangle& triangle : triangles) {
            used[triangle.i0] = true;
            used[triangle.i1] = true;
            used[triangle.i2] = true;
        }
        for (std::size_t outside = 7; outside < points.size(); ++outside) {
            require(
                !used[outside],
                "square polygon retained a point outside its domain");
        }
    }

    {
        const std::vector<Point> points = {
            {0, 0}, {120, 0}, {120, 40},
            {50, 40}, {50, 120}, {0, 120},
            {20, 20}, {80, 20}, {20, 80},
            {80, 80}, {-20, 60},
        };
        const std::vector<std::uint32_t> outer = {0, 1, 2, 3, 4, 5};
        Triangulator triangulator;
        const std::vector<Triangle> triangles =
            triangulator.triangulate_polygon_int(points, outer, {});
        require_valid_polygon_mesh(
            points, triangles, outer, {}, "concave polygon");
        for (const Triangle& triangle : triangles) {
            require(
                triangle.i0 != 9 && triangle.i1 != 9 && triangle.i2 != 9 &&
                    triangle.i0 != 10 && triangle.i1 != 10 &&
                    triangle.i2 != 10,
                "concave polygon retained an excluded point");
        }
    }
}

void test_polygon_holes_and_boundary_chains() {
    const std::vector<Point> points = {
        {0, 0}, {100, 0}, {100, 100}, {0, 100},
        {20, 20}, {40, 20}, {40, 40}, {20, 40},
        {60, 55}, {85, 55}, {85, 80}, {60, 80},
        {50, 20}, {50, 85}, {30, 30}, {70, 65},
        {-20, 50}, {50, 0}, {20, 30},
    };
    const std::vector<std::uint32_t> outer = {0, 3, 2, 1, 0};
    const std::vector<std::vector<std::uint32_t>> holes = {
        {4, 5, 6, 7},
        {8, 11, 10, 9, 8},
    };
    Triangulator triangulator;
    const std::vector<Triangle> triangles =
        triangulator.triangulate_polygon_int(points, outer, holes);
    require_valid_polygon_mesh(
        points,
        triangles,
        outer,
        holes,
        "polygon with two holes and boundary-chain sites");

    std::vector<bool> used(points.size(), false);
    for (const Triangle& triangle : triangles) {
        used[triangle.i0] = true;
        used[triangle.i1] = true;
        used[triangle.i2] = true;
    }
    require(!used[14], "polygon retained a point inside its first hole");
    require(!used[15], "polygon retained a point inside its second hole");
    require(!used[16], "polygon retained a point outside its outer ring");
    require(used[17], "polygon omitted an unlisted outer-boundary point");
    require(used[18], "polygon omitted an unlisted hole-boundary point");
}

void test_polygon_duplicate_representatives() {
    const std::vector<Point> points = {
        {0, 0}, {100, 0}, {100, 100}, {0, 100}, {50, 50},
        {0, 0}, {100, 0}, {100, 100}, {0, 100},
    };
    Triangulator triangulator;
    const std::vector<Triangle> expected =
        triangulator.triangulate_polygon_int(points, {0, 1, 2, 3});
    const std::vector<Triangle> candidate =
        triangulator.triangulate_polygon_int(points, {5, 6, 7, 8, 5});
    require(
        benchmark_support::meshes_equal(expected, candidate),
        "polygon ring indices did not resolve to duplicate representatives");
}

void test_polygon_radial_cases() {
    constexpr double pi = 3.14159265358979323846;
    for (std::uint64_t seed = 0; seed < 120; ++seed) {
        const std::size_t outer_size = 8 + seed % 13;
        const std::size_t hole_size = 5 + seed % 7;
        std::vector<Point> points;
        std::vector<std::uint32_t> outer;
        std::vector<std::uint32_t> hole;
        points.reserve(outer_size + hole_size + 40);
        for (std::size_t i = 0; i < outer_size; ++i) {
            const double angle =
                2.0 * pi * static_cast<double>(i) /
                static_cast<double>(outer_size);
            const double radius =
                700.0 + static_cast<double>((i * 97 + seed * 31) % 180);
            outer.push_back(static_cast<std::uint32_t>(points.size()));
            points.push_back({
                static_cast<std::int32_t>(
                    std::llround(radius * std::cos(angle))),
                static_cast<std::int32_t>(
                    std::llround(radius * std::sin(angle))),
            });
        }
        for (std::size_t i = 0; i < hole_size; ++i) {
            const double angle =
                2.0 * pi * static_cast<double>(i) /
                static_cast<double>(hole_size);
            const double radius =
                130.0 + static_cast<double>((i * 41 + seed * 17) % 45);
            hole.push_back(static_cast<std::uint32_t>(points.size()));
            points.push_back({
                static_cast<std::int32_t>(
                    std::llround(70.0 + radius * std::cos(angle))),
                static_cast<std::int32_t>(
                    std::llround(-40.0 + radius * std::sin(angle))),
            });
        }
        std::vector<Point> background =
            benchmark_support::generate_points(
                Dataset::Uniform, 40, 0xB01700ULL + seed, 2200);
        for (Point& point : background) {
            point.x -= 1100;
            point.y -= 1100;
        }
        points.insert(points.end(), background.begin(), background.end());

        if (seed % 2 == 0) {
            std::reverse(outer.begin(), outer.end());
        }
        if (seed % 3 == 0) {
            std::reverse(hole.begin(), hole.end());
        }
        if (seed % 5 == 0) {
            outer.push_back(outer.front());
            hole.push_back(hole.front());
        }
        const std::vector<std::vector<std::uint32_t>> holes = {hole};
        Triangulator triangulator(seed % 2 == 0 ? 1 : 2);
        const std::vector<Triangle> triangles =
            triangulator.triangulate_polygon_int(points, outer, holes);
        require_valid_polygon_mesh(
            points,
            triangles,
            outer,
            holes,
            "radial polygon case " + std::to_string(seed));
    }
}

void test_parallel_polygon() {
    constexpr std::int32_t width = 250;
    constexpr std::int32_t height = 200;
    std::vector<Point> points;
    points.reserve(static_cast<std::size_t>(width * height));
    for (std::int32_t y = 0; y < height; ++y) {
        for (std::int32_t x = 0; x < width; ++x) {
            points.push_back({x, y});
        }
    }
    const auto index = [](std::int32_t x, std::int32_t y) {
        return static_cast<std::uint32_t>(y * width + x);
    };
    const std::vector<std::uint32_t> outer = {
        index(0, 0),
        index(width - 1, 0),
        index(width - 1, height - 1),
        index(0, height - 1),
    };
    const std::vector<std::vector<std::uint32_t>> holes = {{
        index(80, 70),
        index(170, 70),
        index(170, 130),
        index(80, 130),
    }};

    Triangulator serial(1);
    Triangulator parallel(4);
    const std::vector<Triangle> expected =
        serial.triangulate_polygon_int(points, outer, holes);
    const std::vector<Triangle> candidate =
        parallel.triangulate_polygon_int(points, outer, holes);
    require(
        benchmark_support::meshes_equal(expected, candidate),
        "polygon serial and parallel triangle sets differ");
    require_valid_polygon_mesh(
        points, candidate, outer, holes, "parallel polygon grid");
}

void test_invalid_polygons() {
    const std::vector<Point> square = {
        {0, 0}, {100, 0}, {100, 100}, {0, 100},
        {25, 25}, {75, 25}, {75, 75}, {25, 75},
        {35, 35}, {65, 35}, {65, 65}, {35, 65},
        {-20, 20}, {-5, 20}, {-5, 40}, {-20, 40},
    };
    Triangulator triangulator;
    require_invalid(
        [&] { triangulator.triangulate_polygon_int(square, {0, 1}); },
        "polygon ring with fewer than three indices");
    require_invalid(
        [&] {
            triangulator.triangulate_polygon_int(square, {0, 1, 16});
        },
        "polygon ring index outside the point array");
    require_invalid(
        [&] {
            triangulator.triangulate_polygon_int(square, {0, 1, 2, 1, 3});
        },
        "polygon ring with a repeated coordinate");
    require_invalid(
        [&] {
            triangulator.triangulate_polygon_int(square, {0, 2, 1, 3});
        },
        "self-intersecting polygon ring");

    const std::vector<Point> backtracking = {
        {0, 0}, {100, 0}, {50, 0}, {0, 100},
    };
    require_invalid(
        [&] {
            triangulator.triangulate_polygon_int(
                backtracking, {0, 1, 2, 3});
        },
        "polygon ring with overlapping adjacent edges");

    require_invalid(
        [&] {
            triangulator.triangulate_polygon_int(
                square, {0, 1, 2, 3}, {{12, 13, 14, 15}});
        },
        "hole outside the polygon");
    require_invalid(
        [&] {
            triangulator.triangulate_polygon_int(
                square, {0, 1, 2, 3}, {{0, 4, 5}});
        },
        "hole touching the outer ring");

    const std::vector<Point> crossing_outer = {
        {0, 0}, {100, 0}, {100, 100}, {0, 100},
        {-10, 25}, {25, 25}, {25, 75}, {-10, 75},
    };
    require_invalid(
        [&] {
            triangulator.triangulate_polygon_int(
                crossing_outer, {0, 1, 2, 3}, {{4, 5, 6, 7}});
        },
        "hole crossing the outer ring");
    require_invalid(
        [&] {
            triangulator.triangulate_polygon_int(
                square, {0, 1, 2, 3}, {{4, 5, 6, 7}, {8, 9, 10, 11}});
        },
        "nested polygon holes");

    const std::vector<Point> touching_holes = {
        {0, 0}, {100, 0}, {100, 100}, {0, 100},
        {10, 10}, {50, 10}, {50, 50}, {10, 50},
        {50, 20}, {80, 20}, {80, 60}, {50, 60},
    };
    require_invalid(
        [&] {
            triangulator.triangulate_polygon_int(
                touching_holes,
                {0, 1, 2, 3},
                {{4, 5, 6, 7}, {8, 9, 10, 11}});
        },
        "touching polygon holes");

    const std::vector<Point> crossing_holes = {
        {0, 0}, {100, 0}, {100, 100}, {0, 100},
        {10, 10}, {60, 10}, {60, 60}, {10, 60},
        {40, 40}, {80, 40}, {80, 80}, {40, 80},
    };
    require_invalid(
        [&] {
            triangulator.triangulate_polygon_int(
                crossing_holes,
                {0, 1, 2, 3},
                {{4, 5, 6, 7}, {8, 9, 10, 11}});
        },
        "crossing polygon holes");

    const std::vector<Point> duplicate_coordinate = {
        {0, 0}, {100, 0}, {100, 100}, {0, 100}, {100, 0},
    };
    require_invalid(
        [&] {
            triangulator.triangulate_polygon_int(
                duplicate_coordinate, {0, 1, 2, 4, 3});
        },
        "ring repeating a duplicate representative");

    const std::vector<Triangle> recovered =
        triangulator.triangulate_polygon_int(square, {0, 1, 2, 3});
    require_valid_polygon_mesh(
        square,
        recovered,
        {0, 1, 2, 3},
        {},
        "reuse after invalid polygon input");
}

void test_invalid_inputs() {
    Triangulator triangulator;
    const std::vector<FloatPoint> valid_float_triangle = {
        {0.0F, 0.0F},
        {1.0F, 0.0F},
        {0.0F, 1.0F},
    };
    const auto require_invalid_options =
        [&](const QuantizationOptions& options, const char* label) {
            require_invalid(
                [&] {
                    triangulator.triangulate_float(
                        valid_float_triangle, options);
                },
                label);
        };
    require_invalid(
        [&] {
            triangulator.triangulate_int({{0, 0}, {1, 1}});
        },
        "fewer than three points");
    require_invalid(
        [&] {
            triangulator.triangulate_int(
                {{0, 0}, {0, 0}, {1, 1}, {1, 1}});
        },
        "fewer than three unique points");
    require_invalid(
        [&] {
            triangulator.triangulate_float(
                std::vector<FloatPoint>{{0.0F, 0.0F}, {1.0F, 1.0F}});
        },
        "fewer than three float points");
    require_invalid(
        [&] {
            triangulator.triangulate_float(std::vector<FloatPoint>{
                {0.0F, 0.0F},
                {1.0F, 0.0F},
                {std::numeric_limits<float>::quiet_NaN(), 1.0F},
            });
        },
        "NaN float coordinate");
    require_invalid(
        [&] {
            triangulator.triangulate_float(std::vector<FloatPoint>{
                {0.0F, 0.0F},
                {1.0F, 0.0F},
                {0.0F, std::numeric_limits<float>::infinity()},
            });
        },
        "infinite float coordinate");
    QuantizationOptions invalid_grid;
    invalid_grid.mode = QuantizationMode::GridStep;
    require_invalid_options(invalid_grid, "zero quantization grid step");
    QuantizationOptions invalid_fixed;
    invalid_fixed.mode = QuantizationMode::FixedScale;
    invalid_fixed.origin_x = 0.0;
    invalid_fixed.origin_y = 0.0;
    invalid_fixed.scale = -1.0;
    require_invalid_options(
        invalid_fixed, "negative fixed quantization scale");
    QuantizationOptions invalid_error;
    invalid_error.max_coordinate_error =
        std::numeric_limits<double>::quiet_NaN();
    require_invalid_options(
        invalid_error, "non-finite maximum quantization error");
    QuantizationOptions invalid_mode;
    invalid_mode.mode = static_cast<QuantizationMode>(255);
    require_invalid_options(invalid_mode, "unknown quantization mode");
    QuantizationOptions invalid_collision_policy;
    invalid_collision_policy.collision_policy =
        static_cast<QuantizationCollisionPolicy>(255);
    require_invalid_options(
        invalid_collision_policy,
        "unknown quantization collision policy");
    const float maximum = std::numeric_limits<float>::max();
    require_invalid(
        [&] {
            triangulator.triangulate_float(std::vector<FloatPoint>{
                {-maximum, 0.0F},
                {0.0F, 0.0F},
                {1.0F, 0.0F},
                {2.0F, 0.0F},
            });
        },
        "fewer than three unique quantized float points");
#if defined(__SIZEOF_INT128__)
    require_invalid(
        [&] {
            triangulator.triangulate_int({
                {
                    std::numeric_limits<std::int32_t>::min(),
                    std::numeric_limits<std::int32_t>::min(),
                },
                {
                    std::numeric_limits<std::int32_t>::max(),
                    std::numeric_limits<std::int32_t>::min(),
                },
                {
                    std::numeric_limits<std::int32_t>::min(),
                    std::numeric_limits<std::int32_t>::max(),
                },
            });
        },
        "unsupported full int32 square");
#endif
}

}  // namespace
}  // namespace delaunay32

int main() {
    try {
        delaunay32::test_predicate_selection();
        delaunay32::test_deterministic_cases();
        delaunay32::test_full_integer_result();
        delaunay32::test_float_api();
        delaunay32::test_quantization_options();
        delaunay32::test_float_quantization_collisions();
        delaunay32::test_full_float_result();
        delaunay32::test_float_parallel();
        delaunay32::test_float_collinear_input();
        delaunay32::test_duplicates();
        delaunay32::test_collinear_input();
        delaunay32::test_wide_predicates();
        delaunay32::test_single_constraint_recovery();
        delaunay32::test_constrained_random_cases();
        delaunay32::test_multiple_constraints();
        delaunay32::test_collinear_constraint_chain();
        delaunay32::test_constraint_cycle();
        delaunay32::test_constraint_split_at_unconnected_site();
        delaunay32::test_constraints_meet_at_existing_site();
        delaunay32::test_constrained_duplicates_and_empty_input();
        delaunay32::test_collinear_constrained_input();
        delaunay32::test_parallel_constraints();
        delaunay32::test_invalid_constraints();
        delaunay32::test_polygon_outer_domains();
        delaunay32::test_polygon_holes_and_boundary_chains();
        delaunay32::test_polygon_duplicate_representatives();
        delaunay32::test_polygon_radial_cases();
        delaunay32::test_parallel_polygon();
        delaunay32::test_invalid_polygons();
        delaunay32::test_invalid_inputs();
        std::cout
            << "validation passed: deterministic, duplicate, collinear, "
               "integer/float/full API, quantization options, "
               "predicate-width, constrained, polygon, and parallel cases\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "validation failed: " << error.what() << '\n';
        return 1;
    }
}
