// SPDX-License-Identifier: MIT

#include "support.hpp"

#include "delaunay32/delaunay.hpp"

#include <atomic>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <new>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::atomic<std::size_t> failure_bytes{0};
std::atomic<unsigned> allocations_remaining{0};

bool fail_allocation(std::size_t bytes) noexcept {
    const std::size_t selected =
        failure_bytes.load(std::memory_order_relaxed);
    if (selected == 0 || selected != bytes) {
        return false;
    }
    if (allocations_remaining.fetch_sub(1, std::memory_order_relaxed) != 1) {
        return false;
    }
    failure_bytes.store(0, std::memory_order_relaxed);
    return true;
}

void require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

std::vector<delaunay32::Point> make_points(std::size_t count) {
    std::vector<delaunay32::Point> points;
    points.reserve(count);
    for (std::size_t i = 0; i < count; ++i) {
        points.push_back({
            static_cast<std::int32_t>(i * 7),
            static_cast<std::int32_t>((i * i * 13) % 121),
        });
    }
    return points;
}

void require_recovery(
    const std::vector<delaunay32::Point>& failed_points,
    std::size_t failed_allocation_bytes,
    unsigned failed_allocation,
    const std::vector<delaunay32::Point>& retry_points,
    bool warm,
    const std::string& label) {
    using namespace delaunay32;
    namespace support = benchmark_support;
    const TriangulationOptions options{1, ResultDetail::Full};
    Triangulator reference;
    reference.set_options(options);
    reference.set_points(retry_points);
    const TriangulationResult expected = reference.triangulate();

    Triangulator candidate;
    candidate.set_options(options);
    if (warm) {
        candidate.set_points(make_points(11));
        (void)candidate.triangulate();
    }
    candidate.set_points(failed_points);
    allocations_remaining.store(failed_allocation, std::memory_order_relaxed);
    failure_bytes.store(failed_allocation_bytes, std::memory_order_relaxed);
    bool caught = false;
    try {
        (void)candidate.triangulate();
    } catch (const std::bad_alloc&) {
        caught = true;
    } catch (...) {
        failure_bytes.store(0, std::memory_order_relaxed);
        throw;
    }
    failure_bytes.store(0, std::memory_order_relaxed);
    require(caught, label + ": allocation failure was not injected");
    require(
        allocations_remaining.load(std::memory_order_relaxed) == 0,
        label + ": allocation failed before the intended arena vector");

    for (unsigned repeat = 0; repeat < 2; ++repeat) {
        candidate.set_points(retry_points);
        const TriangulationResult actual = candidate.triangulate();
        require(
            support::meshes_equal(expected.triangles, actual.triangles),
            label + ": recovered triangle set differs");
        require(
            actual.halfedges == expected.halfedges &&
                actual.hull == expected.hull &&
                actual.representatives == expected.representatives,
            label + ": recovered full result differs");
        require(
            actual.report.unique_points == retry_points.size() &&
                actual.report.actual_thread_count == 1,
            label + ": recovered report differs");
    }
}

}  // namespace

// Intercept only a selected allocation size and ordinal. The tests use one
// worker and arm injection immediately before triangulate(), keeping unrelated
// setup allocations out of the failure sequence. No production test API is
// needed, and every case verifies that the intended allocation actually failed.
void* operator new(std::size_t bytes) {
    if (fail_allocation(bytes)) {
        throw std::bad_alloc();
    }
    if (void* memory = std::malloc(bytes == 0 ? 1 : bytes)) {
        return memory;
    }
    throw std::bad_alloc();
}

void operator delete(void* memory) noexcept {
    std::free(memory);
}

void operator delete(void* memory, std::size_t) noexcept {
    std::free(memory);
}

int main() {
    try {
#if defined(DELAUNAY32_TEST_PARALLEL_EDGE_ARENA_DART_LIMIT)
        // The private variant starts below this input's required arena size,
        // so the matching allocation occurs in make_edge()'s serial growth.
        constexpr std::size_t limit =
            DELAUNAY32_TEST_PARALLEL_EDGE_ARENA_DART_LIMIT;
        const auto points = delaunay32::benchmark_support::generate_points(
            delaunay32::benchmark_support::Dataset::Uniform,
            60000,
            0xa110ca7eULL,
            20000);
        for (unsigned failed : {1U, 2U, 3U}) {
            for (std::size_t retry_count : {17U, 60000U, 70000U}) {
                const auto retry =
                    delaunay32::benchmark_support::generate_points(
                        delaunay32::benchmark_support::Dataset::Uniform,
                        retry_count,
                        0xa110ca7eULL,
                        20000);
                require_recovery(
                    points,
                    limit * 2 * sizeof(std::uint32_t),
                    failed,
                    retry,
                    false,
                    "serial arena growth allocation " +
                        std::to_string(failed) + " retry " +
                        std::to_string(retry_count));
            }
        }
#else
        const auto points = make_points(50);
        for (unsigned failed : {1U, 2U, 3U}) {
            for (bool warm : {false, true}) {
                for (std::size_t retry_count : {17U, 50U, 120U}) {
                    require_recovery(
                        points,
                        points.size() * 9 * sizeof(std::uint32_t),
                        failed,
                        make_points(retry_count),
                        warm,
                        "initial arena allocation " +
                            std::to_string(failed) + " retry " +
                            std::to_string(retry_count) +
                            (warm ? " reused" : " fresh"));
                }
            }
        }
#endif
        std::cout << "edge arena allocation recovery passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "edge arena allocation recovery failed: "
                  << error.what() << '\n';
        return 1;
    }
}
