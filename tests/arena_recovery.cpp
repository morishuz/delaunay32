// SPDX-License-Identifier: MIT

#include "support.hpp"

#include "delaunay32/delaunay.hpp"

#include <array>
#include <atomic>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <new>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

std::atomic<std::size_t> failure_bytes{0};
std::atomic<unsigned> allocations_remaining{0};
std::atomic<bool> record_allocation{false};
std::atomic<std::size_t> recorded_bytes{0};

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

template <typename Element>
std::size_t vector_allocation_bytes(std::size_t count) {
    // Measure the request made by this standard library: large vector
    // allocations can include alignment padding (for example, on MSVC).
    std::vector<Element> probe;
    recorded_bytes.store(0, std::memory_order_relaxed);
    record_allocation.store(true, std::memory_order_relaxed);
    try {
        probe.reserve(count);
    } catch (...) {
        record_allocation.store(false, std::memory_order_relaxed);
        throw;
    }
    record_allocation.store(false, std::memory_order_relaxed);
    const std::size_t bytes = recorded_bytes.load(std::memory_order_relaxed);
    require(bytes >= count * sizeof(Element),
            "could not measure the vector allocation size");
    return bytes;
}

std::size_t arena_allocation_bytes(std::size_t dart_count) {
    return vector_allocation_bytes<std::uint32_t>(dart_count);
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
    const std::string& label,
    std::size_t thread_count = 1) {
    using namespace delaunay32;
    namespace support = benchmark_support;
    const TriangulationOptions options{thread_count, ResultDetail::Full};
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
        label + ": allocation failed before the intended vector");

    for (unsigned repeat = 0; repeat < 2; ++repeat) {
        candidate.set_points(retry_points);
        const TriangulationResult actual = candidate.triangulate();
        require(
            support::meshes_equal(expected.triangles, actual.triangles),
            label + ": recovered triangle set differs");
        // Parallel arena allocation can change triangle order between runs.
        // Reconstruct adjacency in the recovered result's own output order.
        std::vector<std::int64_t> adjacency(actual.triangles.size() * 3, -1);
        std::unordered_map<std::uint64_t, std::size_t> first_edges;
        for (std::size_t i = 0; i < actual.triangles.size(); ++i) {
            const Triangle triangle = actual.triangles[i];
            const std::array<std::uint32_t, 3> vertices = {
                triangle.i0, triangle.i1, triangle.i2,
            };
            for (std::size_t local = 0; local < 3; ++local) {
                const std::size_t edge = i * 3 + local;
                const auto found = first_edges.emplace(
                    support::edge_key(
                        vertices[local], vertices[(local + 1) % 3]),
                    edge);
                if (!found.second) {
                    const std::size_t opposite = found.first->second;
                    adjacency[edge] = static_cast<std::int64_t>(opposite);
                    adjacency[opposite] = static_cast<std::int64_t>(edge);
                }
            }
        }
        require(
            actual.halfedges == adjacency &&
                actual.hull == expected.hull &&
                actual.representatives == expected.representatives,
            label + ": recovered full result differs");
        require(
            actual.report.unique_points == retry_points.size() &&
                actual.report.actual_thread_count ==
                    expected.report.actual_thread_count,
            label + ": recovered report differs");
    }
}

}  // namespace

// Intercept only a selected allocation size and ordinal. The tests use one
// worker for arena failures and also cover parallel output allocation. Arm
// injection immediately before triangulate(), keeping setup allocations out
// of the failure sequence. No production test API is needed, and every case
// verifies that the intended allocation actually failed.
void* operator new(std::size_t bytes) {
    if (record_allocation.load(std::memory_order_relaxed)) {
        recorded_bytes.store(bytes, std::memory_order_relaxed);
    }
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
        const std::size_t growth_bytes = arena_allocation_bytes(limit * 2);
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
                    growth_bytes,
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
        const std::size_t initial_bytes =
            arena_allocation_bytes(points.size() * 9);
        for (unsigned failed : {1U, 2U, 3U}) {
            for (bool warm : {false, true}) {
                for (std::size_t retry_count : {17U, 50U, 120U}) {
                    require_recovery(
                        points,
                        initial_bytes,
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

        // Halfedge allocation happens after serial face discovery, or while
        // parallel workers wait for their output slices. Failure must discard
        // any temporary face IDs and release every waiting worker.
        const auto export_points =
            delaunay32::benchmark_support::generate_points(
                delaunay32::benchmark_support::Dataset::Uniform,
                65536, 0xe87047ULL, 20000);
        delaunay32::Triangulator reference;
        reference.set_options({1, delaunay32::ResultDetail::Full});
        reference.set_points(export_points);
        const auto full = reference.triangulate();
        const std::size_t output_bytes =
            vector_allocation_bytes<std::int64_t>(full.halfedges.size());
        for (const std::size_t threads : {1U, 8U}) {
            for (bool warm : {false, true}) {
                require_recovery(
                    export_points, output_bytes, 1, export_points, warm,
                    "full export allocation with " +
                        std::to_string(threads) + " threads",
                    threads);
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
