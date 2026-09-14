// SPDX-License-Identifier: MIT

#pragma once

#include "internal.hpp"

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace delaunay32::detail {

// Append-only primal dart storage. Deleted darts are not recycled during
// construction. Serial allocation grows; parallel allocation uses fixed blocks
// and signals exhaustion so the caller can discard topology and retry serially.
struct EdgeArena {
    // Production benchmarks use about 8.1 darts/point. Nine provides measured
    // headroom, not a geometric upper bound.
    static constexpr std::size_t kDartsPerPoint = 9;
    static constexpr std::uint32_t kBlockDarts = 4096;

    struct Range {
        std::uint32_t first = 0;
        std::uint32_t last = 0;
    };

    // Parallel jobs update cursors independently. Separate adjacent cursors
    // even on machines with 128-byte cache lines.
    struct alignas(128) Cursor {
        std::uint32_t next = 0;
        std::uint32_t end = 0;
        std::uint32_t range_first = 0;
        std::atomic<std::size_t>* block_counter = nullptr;
        std::vector<Range> ranges;
    };

    std::vector<std::uint32_t> origin;
    std::vector<std::uint32_t> next;
    std::vector<std::uint32_t> prev;
    std::vector<Range> ranges;
    std::size_t count = 0;
    std::size_t capacity_limit = 0;

    void resize(std::size_t dart_count);
    static void finish_cursor(Cursor& cursor);
    template <bool ParallelAllocation = false>
    std::uint32_t make_edge(
        std::uint32_t origin_site,
        std::uint32_t destination,
        Cursor* cursor = nullptr);

private:
    void acquire_block(Cursor& cursor);
};

inline void EdgeArena::resize(std::size_t dart_count) {
    // Reserve every backing allocation before changing any logical size.
    // If allocation fails, all three sizes still agree, including when the
    // caller starts a smaller problem after the failed run. Resizing uint32_t
    // elements within reserved capacity cannot allocate or throw.
    origin.reserve(dart_count);
    next.reserve(dart_count);
    prev.reserve(dart_count);
    origin.resize(dart_count);
    next.resize(dart_count);
    prev.resize(dart_count);
}

inline void EdgeArena::acquire_block(Cursor& cursor) {
    finish_cursor(cursor);
    if (cursor.block_counter == nullptr) {
        throw std::logic_error("parallel edge cursor has no block counter");
    }
    std::size_t limit = capacity_limit;
#if defined(DELAUNAY32_TEST_PARALLEL_EDGE_ARENA_DART_LIMIT)
    limit = std::min(
        limit,
        static_cast<std::size_t>(
            DELAUNAY32_TEST_PARALLEL_EDGE_ARENA_DART_LIMIT));
#endif
    const std::size_t first =
        cursor.block_counter->fetch_add(
            kBlockDarts, std::memory_order_relaxed);
    if (first >= limit || limit - first < 2) {
        throw ParallelEdgeArenaExhausted{};
    }
    const std::size_t last =
        std::min(first + kBlockDarts, limit);
    cursor.next = static_cast<std::uint32_t>(first);
    cursor.end = static_cast<std::uint32_t>(last);
    cursor.range_first = cursor.next;
}

inline void EdgeArena::finish_cursor(Cursor& cursor) {
    if (cursor.range_first < cursor.next) {
        cursor.ranges.push_back({cursor.range_first, cursor.next});
    }
    cursor.range_first = cursor.next;
}

template <bool ParallelAllocation>
std::uint32_t EdgeArena::make_edge(
    std::uint32_t origin_site,
    std::uint32_t destination,
    Cursor* cursor) {
    std::uint32_t edge = 0;
    if constexpr (ParallelAllocation) {
        if (cursor == nullptr) {
            throw std::logic_error("parallel edge allocation has no cursor");
        }
        if (cursor->next + 2 > cursor->end) {
            acquire_block(*cursor);
        }
        edge = cursor->next;
        cursor->next += 2;
    } else {
        if (count + 2 > origin.size()) {
            const std::size_t new_size =
                std::max(origin.size() * 2, count + 2);
            resize(new_size);
        }
        edge = static_cast<std::uint32_t>(count);
        count += 2;
    }
    origin[edge] = origin_site;
    origin[edge + 1] = destination;
    next[edge] = edge;
    next[edge + 1] = edge + 1;
    prev[edge] = edge;
    prev[edge + 1] = edge + 1;
    return edge;
}

}  // namespace delaunay32::detail
