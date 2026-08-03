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
        delaunay32::test_invalid_inputs();
        std::cout
            << "validation passed: deterministic, duplicate, collinear, "
               "integer/float/full API, quantization options, "
               "predicate-width, and parallel cases\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "validation failed: " << error.what() << '\n';
        return 1;
    }
}
