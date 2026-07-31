// SPDX-License-Identifier: MIT

#include "delaunator_adapter.hpp"
#include "support.hpp"

#include "delaunay32/delaunay.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
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
            serial.triangulate(points);
        require_reference_match(
            points, candidate, reference, "serial deterministic case");

        if (sizes[i] >= 50000) {
            const std::vector<Triangle> parallel_candidate =
                parallel.triangulate(points);
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
        triangulator.triangulate(points);
    const std::vector<Triangle> array_mesh =
        triangulator.triangulate_int(xs.data(), ys.data(), points.size());
    require(
        benchmark_support::meshes_equal(vector_mesh, array_mesh),
        "vector and struct-of-arrays APIs differ");
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
    std::vector<Triangle> expected = triangulator.triangulate(unique);
    for (Triangle& triangle : expected) {
        triangle.i0 = representative[triangle.i0];
        triangle.i1 = representative[triangle.i1];
        triangle.i2 = representative[triangle.i2];
    }

    const std::vector<Triangle> candidate =
        triangulator.triangulate(duplicated);
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
    const std::vector<Triangle> horizontal = triangulator.triangulate(
        {{-4, 7}, {-1, 7}, {0, 7}, {3, 7}, {9, 7}});
    require(horizontal.empty(), "collinear input produced a triangle");

    const std::vector<Triangle> diagonal = triangulator.triangulate(
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
        const std::vector<Triangle> candidate =
            serial.triangulate(points);
        const std::vector<Triangle> reference =
            delaunator.triangulate(points);
        require_reference_match(
            points, candidate, reference, "wide-predicate case");
    }

    const std::vector<Point> parallel_points =
        benchmark_support::generate_points(
            Dataset::Uniform, 60000, 0xfeedULL, 100000);
    const std::vector<Triangle> serial_mesh =
        serial.triangulate(parallel_points);
    Triangulator parallel(2);
    const std::vector<Triangle> parallel_mesh =
        parallel.triangulate(parallel_points);
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
        serial.triangulate(extreme_strip);
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
            triangulator.triangulate({{0, 0}, {1, 1}});
        },
        "fewer than three points");
    require_invalid(
        [&] {
            triangulator.triangulate(
                {{0, 0}, {0, 0}, {1, 1}, {1, 1}});
        },
        "fewer than three unique points");
    require_invalid(
        [&] {
            triangulator.triangulate_int(nullptr, nullptr, 3);
        },
        "null coordinate arrays");
#if defined(__SIZEOF_INT128__)
    require_invalid(
        [&] {
            triangulator.triangulate({
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
        delaunay32::test_duplicates();
        delaunay32::test_collinear_input();
        delaunay32::test_wide_predicates();
        delaunay32::test_invalid_inputs();
        std::cout
            << "validation passed: deterministic, duplicate, collinear, API, "
               "predicate-width, and parallel cases\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "validation failed: " << error.what() << '\n';
        return 1;
    }
}
