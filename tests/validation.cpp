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

void test_struct_of_arrays_api() {
    const std::vector<Point> points =
        benchmark_support::generate_points(
            Dataset::Uniform, 4096, 0xa11ceULL, 20000);
    std::vector<std::int32_t> xs(points.size());
    std::vector<std::int32_t> ys(points.size());
    for (std::size_t i = 0; i < points.size(); ++i) {
        xs[i] = points[i].x;
        ys[i] = points[i].y;
    }

    Triangulator triangulator;
    const std::vector<Triangle> vector_mesh =
        triangulator.triangulate_int(points);
    const std::vector<Triangle> array_mesh =
        triangulator.triangulate_int(xs.data(), ys.data(), points.size());
    require(
        benchmark_support::meshes_equal(vector_mesh, array_mesh),
        "vector and struct-of-arrays APIs differ");
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
    QuantizationReport report;
    const std::vector<Triangle> vector_mesh =
        triangulator.triangulate_float(points, &report);
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

    std::vector<float> xs;
    std::vector<float> ys;
    xs.reserve(points.size());
    ys.reserve(points.size());
    for (const FloatPoint& point : points) {
        xs.push_back(point.x);
        ys.push_back(point.y);
    }
    QuantizationReport array_report;
    const std::vector<Triangle> array_mesh = triangulator.triangulate_float(
        xs.data(), ys.data(), xs.size(), &array_report);
    require(
        benchmark_support::meshes_equal(vector_mesh, array_mesh),
        "float vector and struct-of-arrays APIs differ");
    require(
        report.origin_x == array_report.origin_x &&
            report.origin_y == array_report.origin_y &&
            report.scale == array_report.scale &&
            report.unique_points == array_report.unique_points,
        "float overloads reported different quantization");

    const std::vector<Triangle> no_report_mesh =
        triangulator.triangulate_float(points);
    require(
        benchmark_support::meshes_equal(vector_mesh, no_report_mesh),
        "optional float report changed the mesh");
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
    QuantizationReport report;
    const std::vector<Triangle> candidate =
        triangulator.triangulate_float(points, &report);
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
    QuantizationReport report;
    const std::vector<Triangle> parallel_mesh =
        parallel.triangulate_float(points, &report);
    require(
        report.unique_points == points.size() &&
            report.collapsed_points == 0,
        "parallel float quantization unexpectedly collapsed points");
    require(
        benchmark_support::meshes_equal(serial_mesh, parallel_mesh),
        "serial and parallel float triangle sets differ");
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

    std::vector<std::int32_t> xs;
    std::vector<std::int32_t> ys;
    xs.reserve(duplicated.size());
    ys.reserve(duplicated.size());
    for (const Point& point : duplicated) {
        xs.push_back(point.x);
        ys.push_back(point.y);
    }
    const std::vector<Triangle> array_candidate =
        triangulator.triangulate_int(xs.data(), ys.data(), xs.size());
    require(
        benchmark_support::meshes_equal(expected, array_candidate),
        "struct-of-arrays duplicate compaction differs");
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
            triangulator.triangulate_int(nullptr, nullptr, 3);
        },
        "null coordinate arrays");
    require_invalid(
        [&] {
            triangulator.triangulate_float(
                std::vector<FloatPoint>{{0.0F, 0.0F}, {1.0F, 1.0F}});
        },
        "fewer than three float points");
    require_invalid(
        [&] {
            triangulator.triangulate_float(nullptr, nullptr, 3);
        },
        "null float coordinate arrays");
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
        delaunay32::test_struct_of_arrays_api();
        delaunay32::test_float_api();
        delaunay32::test_float_quantization_collisions();
        delaunay32::test_float_parallel();
        delaunay32::test_float_collinear_input();
        delaunay32::test_duplicates();
        delaunay32::test_collinear_input();
        delaunay32::test_wide_predicates();
        delaunay32::test_invalid_inputs();
        std::cout
            << "validation passed: deterministic, duplicate, collinear, "
               "integer/float API, predicate-width, and parallel cases\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "validation failed: " << error.what() << '\n';
        return 1;
    }
}
