// SPDX-License-Identifier: MIT

#include "support.hpp"

#include "delaunay32/delaunay.hpp"
#include "delaunay32/quantization.hpp"

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
#include <thread>
#include <unordered_set>
#include <utility>
#include <vector>

namespace delaunay32 {
namespace {

namespace support = benchmark_support;

void require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void configure(
    Triangulator& triangulator,
    std::size_t thread_count,
    ResultDetail detail = ResultDetail::Triangles) {
    TriangulationOptions options;
    options.thread_count = thread_count;
    options.result_detail = detail;
    triangulator.set_options(options);
}

std::vector<Triangle> triangulate_points(
    Triangulator& triangulator,
    const std::vector<Point>& points) {
    triangulator.set_points(points);
    return triangulator.triangulate().triangles;
}

TriangulationResult triangulate_full(
    Triangulator& triangulator,
    const std::vector<Point>& points,
    std::size_t thread_count) {
    configure(triangulator, thread_count, ResultDetail::Full);
    triangulator.set_points(points);
    return triangulator.triangulate();
}

template <typename Operation>
void require_invalid(Operation&& operation, const std::string& label) {
    try {
        operation();
    } catch (const std::invalid_argument&) {
        return;
    }
    throw std::runtime_error(label + ": expected invalid_argument");
}

bool ordered_meshes_equal(
    const std::vector<Triangle>& a,
    const std::vector<Triangle>& b) {
    if (a.size() != b.size()) {
        return false;
    }
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i].i0 != b[i].i0 ||
            a[i].i1 != b[i].i1 ||
            a[i].i2 != b[i].i2) {
            return false;
        }
    }
    return true;
}

bool all_collinear(const std::vector<Point>& points) {
    for (std::size_t i = 2; i < points.size(); ++i) {
        if (support::orient(points[0], points[1], points[i]) != 0) {
            return false;
        }
    }
    return true;
}

void require_valid_mesh(
    const std::vector<Point>& points,
    const std::vector<Triangle>& triangles,
    const std::string& label) {
    if (triangles.empty()) {
        require(all_collinear(points), label + ": unexpected empty mesh");
        return;
    }

    std::string error;
    require(
        support::validate_mesh(points, triangles, error),
        label + ": " + error);
    for (const Triangle& triangle : triangles) {
        require(
            support::orient(
                points[triangle.i0],
                points[triangle.i1],
                points[triangle.i2]) > 0,
            label + ": triangle is not counterclockwise");
    }
}

void require_same_mesh(
    const std::vector<Triangle>& expected,
    const std::vector<Triangle>& candidate,
    const std::string& label) {
    require(
        support::meshes_equal(expected, candidate),
        label + ": triangle set differs from expected mesh");
}

void test_seeded_serial_properties() {
    constexpr std::array<std::size_t, 13> sizes = {
        3, 4, 5, 15, 16, 17, 31, 32, 33, 127, 255, 1023, 4095,
    };
    constexpr std::array<support::Dataset, 3> datasets = {
        support::Dataset::Uniform,
        support::Dataset::Clustered,
        support::Dataset::Diagonal,
    };

    Triangulator triangulator;
    configure(triangulator, 1);
    for (std::size_t dataset_index = 0;
         dataset_index < datasets.size();
         ++dataset_index) {
        for (std::size_t size_index = 0;
             size_index < sizes.size();
             ++size_index) {
            for (std::size_t seed_index = 0; seed_index < 6; ++seed_index) {
                const std::uint64_t seed =
                    0x51eed000ULL + dataset_index * 0x10000ULL +
                    size_index * 0x100ULL + seed_index;
                std::vector<Point> points = support::generate_points(
                    datasets[dataset_index],
                    sizes[size_index],
                    seed,
                    1000);
                if (seed_index == 1) {
                    for (Point& point : points) {
                        point.x -= 1000000000;
                        point.y += 1000000000;
                    }
                }

                const std::string label =
                    std::string(support::dataset_name(
                        datasets[dataset_index])) +
                    " size " + std::to_string(sizes[size_index]) +
                    " seed " + std::to_string(seed_index);
                const std::vector<Triangle> first =
                    triangulate_points(triangulator, points);
                require_valid_mesh(points, first, label);
                const std::vector<Triangle> repeated =
                    triangulate_points(triangulator, points);
                require(
                    ordered_meshes_equal(first, repeated),
                    label + ": repeated serial output order changed");
            }
        }
    }
}

void test_parallel_thresholds() {
    const std::vector<Point> master = support::generate_points(
        support::Dataset::Uniform,
        100001,
        0x7a125eedULL,
        1000);
    constexpr std::array<std::size_t, 5> sizes = {
        49999, 50000, 99999, 100000, 100001,
    };

    Triangulator serial;
    configure(serial, 1);
    Triangulator parallel;
    configure(parallel, 2);
    for (const std::size_t size : sizes) {
        const std::vector<Point> points(
            master.begin(),
            master.begin() + static_cast<std::ptrdiff_t>(size));
        const std::vector<Triangle> expected =
            triangulate_points(serial, points);
        const std::vector<Triangle> candidate =
            triangulate_points(parallel, points);
        const std::string label =
            "parallel threshold size " + std::to_string(size);
        require_same_mesh(expected, candidate, label);
        require_valid_mesh(points, candidate, label);
        // Parallel edge blocks are acquired according to worker scheduling,
        // so repeatability here applies to topology rather than vector order.
        require_same_mesh(
            candidate,
            triangulate_points(parallel, points),
            label + " repeat");
    }
}

void test_thread_counts_and_reconfiguration() {
    const std::vector<Point> points = support::generate_points(
        support::Dataset::Clustered,
        60000,
        0x7a12c0deULL,
        2000);
    Triangulator serial;
    configure(serial, 1);
    const std::vector<Triangle> expected = triangulate_points(serial, points);
    require_valid_mesh(points, expected, "thread-count serial baseline");

    constexpr std::array<std::size_t, 5> thread_counts = {1, 2, 4, 8, 0};
    for (const std::size_t thread_count : thread_counts) {
        Triangulator triangulator;
        configure(triangulator, thread_count);
        const std::vector<Triangle> candidate =
            triangulate_points(triangulator, points);
        require_same_mesh(
            expected,
            candidate,
            "fresh requested thread count " +
                std::to_string(thread_count));
        require_same_mesh(
            candidate,
            triangulate_points(triangulator, points),
            "fresh requested thread count " +
                std::to_string(thread_count) + " repeat");
    }

    Triangulator reused;
    configure(reused, 1);
    for (const std::size_t thread_count : thread_counts) {
        configure(reused, thread_count);
        const std::vector<Triangle> candidate =
            triangulate_points(reused, points);
        require_same_mesh(
            expected,
            candidate,
            "reconfigured thread count " +
                std::to_string(thread_count));
        require_same_mesh(
            candidate,
            triangulate_points(reused, points),
            "reconfigured thread count " +
                std::to_string(thread_count) + " repeat");
    }
}

void test_grid_and_duplicates() {
    std::vector<Point> grid;
    grid.reserve(60000);
    for (std::int32_t y = 0; y < 60; ++y) {
        for (std::int32_t x = 0; x < 1000; ++x) {
            grid.push_back({x, y});
        }
    }
    Triangulator serial;
    configure(serial, 1);
    Triangulator parallel;
    configure(parallel, 4);
    const std::vector<Triangle> grid_expected =
        triangulate_points(serial, grid);
    const std::vector<Triangle> grid_candidate =
        triangulate_points(parallel, grid);
    require_same_mesh(
        grid_expected,
        grid_candidate,
        "highly skewed cocircular grid");
    require_valid_mesh(
        grid,
        grid_candidate,
        "highly skewed cocircular grid");
    require_same_mesh(
        grid_candidate,
        triangulate_points(parallel, grid),
        "highly skewed cocircular grid repeat");

    const std::vector<Point> unique = support::generate_points(
        support::Dataset::Uniform,
        20000,
        0xd0011ca7eULL,
        1000);
    std::vector<Point> duplicated = unique;
    duplicated.insert(duplicated.end(), unique.rbegin(), unique.rend());
    duplicated.insert(duplicated.end(), unique.begin(), unique.end());

    const std::vector<Triangle> unique_mesh =
        triangulate_points(serial, unique);
    const std::vector<Triangle> duplicate_mesh =
        triangulate_points(parallel, duplicated);
    require_same_mesh(
        unique_mesh,
        duplicate_mesh,
        "duplicate-heavy parallel compaction");
    require_valid_mesh(
        unique,
        duplicate_mesh,
        "duplicate-heavy parallel compaction");
    require_same_mesh(
        duplicate_mesh,
        triangulate_points(parallel, duplicated),
        "duplicate-heavy parallel compaction repeat");
}

void test_ring_and_nearly_collinear_inputs() {
    constexpr std::size_t ring_candidate_count = 2048;
    constexpr double tau = 6.283185307179586476925286766559;
    std::vector<Point> ring;
    ring.reserve(ring_candidate_count);
    std::unordered_set<std::uint64_t> occupied;
    occupied.reserve(ring_candidate_count * 2);
    for (std::size_t i = 0; i < ring_candidate_count; ++i) {
        const double angle = tau * static_cast<double>(i) /
                             static_cast<double>(ring_candidate_count);
        const Point point = {
            static_cast<std::int32_t>(
                std::llround(1000.0 * std::cos(angle))),
            static_cast<std::int32_t>(
                std::llround(1000.0 * std::sin(angle))),
        };
        if (occupied.insert(support::point_key(point.x, point.y)).second) {
            ring.push_back(point);
        }
    }
    Triangulator serial;
    configure(serial, 1);
    const std::vector<Triangle> ring_mesh = triangulate_points(serial, ring);
    require_valid_mesh(ring, ring_mesh, "ring-like input");

    std::vector<Point> nearly_collinear;
    nearly_collinear.reserve(60000);
    for (std::int32_t band = 0; band < 60; ++band) {
        for (std::int32_t x = 0; x < 1000; ++x) {
            nearly_collinear.push_back({x, x + band});
        }
    }
    Triangulator parallel;
    configure(parallel, 4);
    const std::vector<Triangle> expected =
        triangulate_points(serial, nearly_collinear);
    const std::vector<Triangle> candidate =
        triangulate_points(parallel, nearly_collinear);
    require_same_mesh(
        expected,
        candidate,
        "nearly collinear parallel bands");
    require_valid_mesh(
        nearly_collinear,
        candidate,
        "nearly collinear parallel bands");
    require_same_mesh(
        candidate,
        triangulate_points(parallel, nearly_collinear),
        "nearly collinear parallel bands repeat");
}

void test_moves_and_exception_recovery() {
    const std::vector<Point> points = support::generate_points(
        support::Dataset::Diagonal,
        60000,
        0xc0de5eedULL,
        1000);
    Triangulator serial;
    configure(serial, 1);
    const std::vector<Triangle> expected = triangulate_points(serial, points);

    Triangulator source;
    configure(source, 2);
    require_same_mesh(
        expected,
        triangulate_points(source, points),
        "pre-move triangulator");
    Triangulator moved(std::move(source));
    require_same_mesh(
        expected,
        triangulate_points(moved, points),
        "move-constructed triangulator");

    Triangulator assigned;
    configure(assigned, 1);
    assigned = std::move(moved);
    require_same_mesh(
        expected,
        triangulate_points(assigned, points),
        "move-assigned triangulator");

    require_invalid(
        [&] { triangulate_points(assigned, {{0, 0}, {1, 1}}); },
        "short integer input");
    require_same_mesh(
        expected,
        triangulate_points(assigned, points),
        "reuse after short-input exception");

    require_invalid(
        [&] {
            quantize(std::vector<FloatPoint>{
                {0.0F, 0.0F},
                {1.0F, 0.0F},
                {std::numeric_limits<float>::quiet_NaN(), 1.0F},
            });
        },
        "NaN float input");
    require_same_mesh(
        expected,
        triangulate_points(assigned, points),
        "reuse after float-input exception");
}

void test_independent_concurrent_instances() {
    const std::vector<Point> points = support::generate_points(
        support::Dataset::Uniform,
        60000,
        0xc0ac077eULL,
        1000);
    Triangulator serial;
    configure(serial, 1);
    const std::vector<Triangle> expected = triangulate_points(serial, points);

    constexpr std::array<std::size_t, 3> thread_counts = {1, 2, 4};
    std::array<std::vector<Triangle>, thread_counts.size()> outputs;
    std::array<std::exception_ptr, thread_counts.size()> exceptions;
    std::vector<std::thread> callers;
    callers.reserve(thread_counts.size());
    for (std::size_t i = 0; i < thread_counts.size(); ++i) {
        callers.emplace_back([&, i] {
            try {
                Triangulator triangulator;
                configure(triangulator, thread_counts[i]);
                if (i + 1 == thread_counts.size()) {
                    TriangulationResult result =
                        triangulate_full(
                            triangulator, points, thread_counts[i]);
                    require(
                        result.halfedges.size() ==
                                result.triangles.size() * 3 &&
                            !result.hull.empty() &&
                            result.representatives.size() == points.size(),
                        "concurrent full result omitted auxiliary data");
                    outputs[i] = std::move(result.triangles);
                } else {
                    outputs[i] = triangulate_points(triangulator, points);
                }
            } catch (...) {
                exceptions[i] = std::current_exception();
            }
        });
    }
    for (std::thread& caller : callers) {
        caller.join();
    }

    for (std::size_t i = 0; i < thread_counts.size(); ++i) {
        if (exceptions[i] != nullptr) {
            std::rethrow_exception(exceptions[i]);
        }
        require_same_mesh(
            expected,
            outputs[i],
            "independent concurrent instance " + std::to_string(i));
    }
}

}  // namespace
}  // namespace delaunay32

int main() {
    try {
        delaunay32::test_seeded_serial_properties();
        delaunay32::test_parallel_thresholds();
        delaunay32::test_thread_counts_and_reconfiguration();
        delaunay32::test_grid_and_duplicates();
        delaunay32::test_ring_and_nearly_collinear_inputs();
        delaunay32::test_moves_and_exception_recovery();
        delaunay32::test_independent_concurrent_instances();
        std::cout
            << "stress validation passed: seeded, threshold, thread-count, "
               "grid, ring, nearly-collinear, duplicate, move, recovery, "
               "and concurrent cases\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "stress validation failed: " << error.what() << '\n';
        return 1;
    }
}
