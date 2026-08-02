// SPDX-License-Identifier: MIT

#include "delaunay32/delaunay.hpp"
#include "internal.hpp"

#include <algorithm>
#include <atomic>
#include <exception>
#include <mutex>
#include <stdexcept>
#include <vector>

#if defined(_MSC_VER)
#define DELAUNAY32_ALWAYS_INLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__)
#define DELAUNAY32_ALWAYS_INLINE inline __attribute__((always_inline))
#else
#define DELAUNAY32_ALWAYS_INLINE inline
#endif

namespace delaunay32 {
using detail::ThreadBarrier;

// Keep topology and predicates in one translation unit so the recursive
// small-leaf merge can inline the exact predicate specialization.
template <
    bool WidePredicates,
    bool ParallelAllocation>
Triangulator::HullEdges Triangulator::build_range(
    std::size_t first,
    std::size_t last,
    EdgeCursor* cursor) {
    const std::size_t count = last - first;
    if (count == 2) {
        const std::uint32_t edge =
            make_edge<ParallelAllocation>(
            static_cast<std::uint32_t>(first),
            static_cast<std::uint32_t>(first + 1),
            cursor);
        return {edge, sym(edge)};
    }

    if (count == 3) {
        const std::uint32_t a =
            make_edge<ParallelAllocation>(
            static_cast<std::uint32_t>(first),
            static_cast<std::uint32_t>(first + 1),
            cursor);
        const std::uint32_t b =
            make_edge<ParallelAllocation>(
            static_cast<std::uint32_t>(first + 1),
            static_cast<std::uint32_t>(first + 2),
            cursor);
        splice(sym(a), b);

        const std::int64_t winding = orient(
            static_cast<std::uint32_t>(first),
            static_cast<std::uint32_t>(first + 1),
            static_cast<std::uint32_t>(first + 2));
        if (winding > 0) {
            connect<ParallelAllocation>(b, a, cursor);
            return {a, sym(b)};
        }
        if (winding < 0) {
            const std::uint32_t c =
                connect<ParallelAllocation>(b, a, cursor);
            return {sym(c), c};
        }
        return {a, sym(b)};
    }
    const std::size_t middle = first + count / 2;
    HullEdges left =
        build_range<WidePredicates, ParallelAllocation>(
            first, middle, cursor);
    HullEdges right =
        build_range<WidePredicates, ParallelAllocation>(
            middle, last, cursor);
    return merge_hulls_inline<WidePredicates, ParallelAllocation>(
        left, right, cursor);
}

Triangulator::MortonSplit
Triangulator::find_morton_split(
    std::size_t first,
    std::size_t last) const {
    const std::uint32_t differing =
        morton_keys_[first] ^ morton_keys_[last - 1];
    if (differing == 0) {
        return {};
    }
    unsigned split_bit = 0;
    for (std::uint32_t bits = differing; (bits >>= 1U) != 0;) {
        ++split_bit;
    }
    const std::uint32_t mask = 1U << split_bit;
    std::size_t low = first;
    std::size_t high = last;
    while (low < high) {
        const std::size_t middle = low + (high - low) / 2;
        if ((morton_keys_[middle] & mask) == 0) {
            low = middle + 1;
        } else {
            high = middle;
        }
    }
    if (low - first < 2 || last - low < 2) {
        return {};
    }
    return {low, split_bit, true};
}

template <
    bool WidePredicates,
    bool ParallelAllocation>
Triangulator::DirectionalHulls
Triangulator::build_morton_range(
    std::size_t first,
    std::size_t last,
    EdgeCursor* cursor) {
    const auto build_leaf = [&] {
        std::sort(
            points_.begin() + static_cast<std::ptrdiff_t>(first),
            points_.begin() + static_cast<std::ptrdiff_t>(last),
            SiteLessXY{});
        bool leaf_uses_int64 = false;
        if constexpr (WidePredicates) {
            // Leaf ranges are disjoint under parallel allocation, so local
            // int64 lifts are safe to compute concurrently. Upper merges keep
            // using wide_lifts_ with the global predicate width.
            std::int32_t leaf_min_y = points_[first].y;
            std::int32_t leaf_max_y = leaf_min_y;
            for (std::size_t i = first + 1; i < last; ++i) {
                leaf_min_y = std::min(leaf_min_y, points_[i].y);
                leaf_max_y = std::max(leaf_max_y, points_[i].y);
            }
            const std::int64_t leaf_min_x = points_[first].x;
            const std::int64_t leaf_max_x = points_[last - 1].x;
            leaf_uses_int64 =
                leaf_max_x - leaf_min_x <= kFastCoordinateSpan &&
                static_cast<std::int64_t>(leaf_max_y) - leaf_min_y <=
                    kFastCoordinateSpan;

            for (std::size_t i = first; i < last; ++i) {
                const std::uint64_t x = static_cast<std::uint64_t>(
                    static_cast<std::int64_t>(points_[i].x) - min_x_);
                const std::uint64_t y = static_cast<std::uint64_t>(
                    static_cast<std::int64_t>(points_[i].y) - min_y_);
                wide_lifts_[i] = x * x + y * y;
                if (leaf_uses_int64) {
                    const std::uint32_t local_x = static_cast<std::uint32_t>(
                        static_cast<std::int64_t>(points_[i].x) - leaf_min_x);
                    const std::uint32_t local_y = static_cast<std::uint32_t>(
                        static_cast<std::int64_t>(points_[i].y) - leaf_min_y);
                    points_[i].lift = local_x * local_x + local_y * local_y;
                }
            }
        }
        const HullEdges hull = [&] {
            if constexpr (WidePredicates) {
                if (leaf_uses_int64) {
                    return build_range<false, ParallelAllocation>(
                        first, last, cursor);
                }
            }
            return build_range<WidePredicates, ParallelAllocation>(
                first, last, cursor);
        }();
        return scan_directional_hulls(sym(hull.left_outer));
    };

    const std::size_t count = last - first;
    if (count <= kMortonLeafSize) {
        return build_leaf();
    }

    const MortonSplit split = find_morton_split(first, last);
    if (!split.valid) {
        return build_leaf();
    }
    const std::size_t middle = split.middle;

    const DirectionalHulls left =
        build_morton_range<WidePredicates, ParallelAllocation>(
            first, middle, cursor);
    const DirectionalHulls right =
        build_morton_range<WidePredicates, ParallelAllocation>(
            middle, last, cursor);
    const bool horizontal = (split.split_bit & 1U) != 0;
    const HullEdges merged =
        merge_hulls<WidePredicates, ParallelAllocation>(
            horizontal ? left.y : left.x,
            horizontal ? right.y : right.x,
            cursor);
    if (horizontal) {
        return scan_merged_hulls<true>(sym(merged.left_outer), merged);
    }
    return scan_merged_hulls<false>(sym(merged.left_outer), merged);
}

template <
    bool WidePredicates,
    bool ParallelAllocation>
DELAUNAY32_ALWAYS_INLINE
Triangulator::HullEdges
Triangulator::merge_hulls_inline(
    HullEdges left,
    HullEdges right,
    EdgeCursor* cursor) {
    std::uint32_t ldi = left.right_outer;
    std::uint32_t rdi = right.left_outer;

    while (true) {
        if (left_of(org(rdi), ldi)) {
            ldi = lnext(ldi);
        } else if (right_of(org(ldi), rdi)) {
            rdi = onext(sym(rdi));
        } else {
            break;
        }
    }

    std::uint32_t base =
        connect<ParallelAllocation>(sym(rdi), ldi, cursor);
    if (org(ldi) == org(left.left_outer)) {
        left.left_outer = sym(base);
    }
    if (org(rdi) == org(right.right_outer)) {
        right.right_outer = base;
    }

    // Left/right prune loops stay expanded: factoring them into lambdas
    // measurably regresses the leaf hot path on Apple Silicon.
    while (true) {
        std::uint32_t lcand = onext(sym(base));
        bool valid_l = right_of(dest(lcand), base);
        if (valid_l) {
            bool deleted = false;
            while (in_circle<WidePredicates>(
                dest(base),
                org(base),
                dest(lcand),
                dest(onext(lcand)))) {
                const std::uint32_t next = onext(lcand);
                delete_edge(lcand);
                lcand = next;
                deleted = true;
            }
            if (deleted) {
                valid_l = right_of(dest(lcand), base);
            }
        }

        std::uint32_t rcand = oprev(base);
        bool valid_r = right_of(dest(rcand), base);
        if (valid_r) {
            bool deleted = false;
            while (in_circle<WidePredicates>(
                dest(base),
                org(base),
                dest(rcand),
                dest(oprev(rcand)))) {
                const std::uint32_t next = oprev(rcand);
                delete_edge(rcand);
                rcand = next;
                deleted = true;
            }
            if (deleted) {
                valid_r = right_of(dest(rcand), base);
            }
        }

        if (!valid_l && !valid_r) {
            break;
        }
        if (!valid_l ||
            (valid_r &&
             in_circle<WidePredicates>(
                 dest(lcand),
                 org(lcand),
                 org(rcand),
                 dest(rcand)))) {
            base = connect<ParallelAllocation>(rcand, sym(base), cursor);
        } else {
            base = connect<ParallelAllocation>(
                sym(base), sym(lcand), cursor);
        }
    }
    return {left.left_outer, right.right_outer};
}

template <
    bool WidePredicates,
    bool ParallelAllocation>
Triangulator::HullEdges
Triangulator::merge_hulls(
    HullEdges left,
    HullEdges right,
    EdgeCursor* cursor) {
    // Leaf recursion expands the body through merge_hulls_inline(). Keeping
    // this wrapper out of other callers avoids duplicating the large kernel
    // throughout the Morton and parallel merge machinery.
    return merge_hulls_inline<WidePredicates, ParallelAllocation>(
        left, right, cursor);
}

std::size_t Triangulator::add_parallel_node(
    std::size_t first,
    std::size_t last,
    std::size_t target_size,
    std::vector<ParallelNode>& nodes,
    std::vector<std::size_t>& leaves) const {
    const std::size_t index = nodes.size();
    ParallelNode new_node;
    new_node.first = first;
    new_node.last = last;
    nodes.push_back(new_node);
    if (last - first <= target_size) {
        leaves.push_back(index);
        return index;
    }

    const MortonSplit split = find_morton_split(first, last);
    if (!split.valid) {
        leaves.push_back(index);
        return index;
    }

    const std::size_t left = add_parallel_node(
        first,
        split.middle,
        target_size,
        nodes,
        leaves);
    const std::size_t right = add_parallel_node(
        split.middle,
        last,
        target_size,
        nodes,
        leaves);
    ParallelNode& node = nodes[index];
    node.left = left;
    node.right = right;
    node.split_bit = split.split_bit;
    node.leaf = false;
    return index;
}

template <bool WidePredicates>
Triangulator::DirectionalHulls
Triangulator::build_parallel(
    std::size_t thread_count,
    detail::WorkerTeam& workers) {
    const std::size_t target_jobs =
        std::max<std::size_t>(2, thread_count * kParallelJobsPerThread);
    const std::size_t target_size = std::max(
        kParallelMinJobPoints,
        (points_.size() + target_jobs - 1) / target_jobs);

    std::vector<ParallelNode> nodes;
    std::vector<std::size_t> leaves;
    nodes.reserve(target_jobs * 2);
    leaves.reserve(target_jobs);
    add_parallel_node(0, points_.size(), target_size, nodes, leaves);
    if (leaves.size() < 2) {
        active_thread_count_ = 1;
        const DirectionalHulls hull =
            build_morton_range<WidePredicates>(0, points_.size());
        edge_ranges_.push_back(
            {0, static_cast<std::uint32_t>(edge_count_)});
        return hull;
    }

    std::atomic<std::size_t> next_edge_block{0};
    std::vector<EdgeCursor> leaf_cursors(leaves.size());
    std::vector<EdgeCursor> node_cursors(nodes.size());
    for (EdgeCursor& cursor : leaf_cursors) {
        cursor.block_counter = &next_edge_block;
    }
    for (EdgeCursor& cursor : node_cursors) {
        cursor.block_counter = &next_edge_block;
    }
    std::vector<std::size_t> node_heights(nodes.size(), 0);
    std::vector<std::vector<std::size_t>> merge_levels(1);
    for (std::size_t i = nodes.size(); i-- > 0;) {
        const ParallelNode& node = nodes[i];
        if (node.leaf) {
            continue;
        }
        const std::size_t height =
            1 + std::max(
                    node_heights[node.left],
                    node_heights[node.right]);
        node_heights[i] = height;
        if (merge_levels.size() <= height) {
            merge_levels.resize(height + 1);
        }
        merge_levels[height].push_back(i);
    }
    std::atomic<std::size_t> next_job{0};
    std::atomic<bool> cancelled{false};
    std::exception_ptr first_exception;
    std::mutex exception_mutex;
    const std::size_t worker_count =
        std::min(thread_count, leaves.size());
    ThreadBarrier barrier(worker_count);
    const auto record_exception = [&] {
        cancelled.store(true, std::memory_order_relaxed);
        std::lock_guard<std::mutex> lock(exception_mutex);
        if (first_exception == nullptr) {
            first_exception = std::current_exception();
        }
    };
    const auto run_jobs = [&](std::size_t worker_index) {
        // Parallel topology proceeds bottom-up:
        //
        //   independent leaves -> barrier -> height-1 merges -> barrier -> ...
        //
        // Worker 0 resets the shared job cursor before each level. The first
        // level barrier publishes that reset; the second ensures every child
        // hull is complete before a parent level starts.
        while (!cancelled.load(std::memory_order_relaxed)) {
            const std::size_t job =
                next_job.fetch_add(1, std::memory_order_relaxed);
            if (job >= leaves.size()) {
                break;
            }
            EdgeCursor& cursor = leaf_cursors[job];
            try {
                ParallelNode& node = nodes[leaves[job]];
                node.hull =
                    build_morton_range<WidePredicates, true>(
                        node.first,
                        node.last,
                        &cursor);
                finish_edge_cursor(cursor);
            } catch (...) {
                finish_edge_cursor(cursor);
                record_exception();
                break;
            }
        }
        barrier.wait();

        for (std::size_t level = 1;
             level < merge_levels.size();
             ++level) {
            if (worker_index == 0) {
                next_job.store(0, std::memory_order_relaxed);
            }
            barrier.wait();
            while (!cancelled.load(std::memory_order_relaxed)) {
                const std::size_t job =
                    next_job.fetch_add(1, std::memory_order_relaxed);
                if (job >= merge_levels[level].size()) {
                    break;
                }
                const std::size_t node_index =
                    merge_levels[level][job];
                ParallelNode& node = nodes[node_index];
                EdgeCursor& cursor = node_cursors[node_index];
                try {
                    const DirectionalHulls& left =
                        nodes[node.left].hull;
                    const DirectionalHulls& right =
                        nodes[node.right].hull;
                    const bool horizontal =
                        (node.split_bit & 1U) != 0;
                    const HullEdges merged =
                        merge_hulls<WidePredicates, true>(
                            horizontal ? left.y : left.x,
                            horizontal ? right.y : right.x,
                            &cursor);
                    node.hull =
                        horizontal
                            ? scan_merged_hulls<true>(
                                  sym(merged.left_outer), merged)
                            : scan_merged_hulls<false>(
                                  sym(merged.left_outer), merged);
                    finish_edge_cursor(cursor);
                } catch (...) {
                    finish_edge_cursor(cursor);
                    record_exception();
                    break;
                }
            }
            barrier.wait();
        }
    };

    workers.run(worker_count, run_jobs);
    if (first_exception != nullptr) {
        std::rethrow_exception(first_exception);
    }
    edge_count_ = 0;
    const auto collect_ranges = [&](const EdgeCursor& cursor) {
        for (const EdgeRange range : cursor.ranges) {
            edge_ranges_.push_back(range);
            edge_count_ += range.last - range.first;
        }
    };
    for (const EdgeCursor& cursor : leaf_cursors) {
        collect_ranges(cursor);
    }
    for (const EdgeCursor& cursor : node_cursors) {
        collect_ranges(cursor);
    }
    std::sort(
        edge_ranges_.begin(),
        edge_ranges_.end(),
        [](const EdgeRange& a, const EdgeRange& b) {
            return a.first < b.first;
        });
    return nodes.front().hull;
}

Triangulator::DirectionalHulls
Triangulator::scan_directional_hulls(
    std::uint32_t outer_seed) const {
    DirectionalHulls result;
    bool have_x_left = false;
    bool have_x_right = false;
    bool have_y_left = false;
    bool have_y_right = false;
    std::int32_t x_left_x = 0;
    std::int32_t x_left_y = 0;
    std::int32_t x_right_x = 0;
    std::int32_t x_right_y = 0;
    std::int32_t y_left_y = 0;
    std::int32_t y_left_x = 0;
    std::int32_t y_right_y = 0;
    std::int32_t y_right_x = 0;

    std::uint32_t outer = outer_seed;
    do {
        const std::uint32_t outer_origin = org(outer);
        const std::uint32_t outer_destination = dest(outer);
        const Site& origin = points_[outer_origin];
        const Site& destination = points_[outer_destination];

        if (!have_x_left ||
            destination.x < x_left_x ||
            (destination.x == x_left_x && destination.y < x_left_y)) {
            have_x_left = true;
            x_left_x = destination.x;
            x_left_y = destination.y;
            result.x.left_outer = sym(outer);
        }
        if (!have_x_right ||
            origin.x > x_right_x ||
            (origin.x == x_right_x && origin.y < x_right_y)) {
            have_x_right = true;
            x_right_x = origin.x;
            x_right_y = origin.y;
            result.x.right_outer = outer;
        }

        if (!have_y_left ||
            destination.y < y_left_y ||
            (destination.y == y_left_y && destination.x > y_left_x)) {
            have_y_left = true;
            y_left_y = destination.y;
            y_left_x = destination.x;
            result.y.left_outer = sym(outer);
        }
        if (!have_y_right ||
            origin.y > y_right_y ||
            (origin.y == y_right_y && origin.x > y_right_x)) {
            have_y_right = true;
            y_right_y = origin.y;
            y_right_x = origin.x;
            result.y.right_outer = outer;
        }

        outer = lnext(outer);
    } while (outer != outer_seed);
    return result;
}

template <bool Horizontal>
Triangulator::DirectionalHulls
Triangulator::scan_merged_hulls(
    std::uint32_t outer_seed,
    HullEdges split_hull) const {
    DirectionalHulls result;
    if constexpr (Horizontal) {
        result.y = split_hull;
    } else {
        result.x = split_hull;
    }

    bool have_left = false;
    bool have_right = false;
    std::int32_t left_primary = 0;
    std::int32_t left_secondary = 0;
    std::int32_t right_primary = 0;
    std::int32_t right_secondary = 0;
    std::uint32_t outer = outer_seed;
    do {
        const Site& origin = points_[org(outer)];
        const Site& destination = points_[dest(outer)];
        if constexpr (Horizontal) {
            if (!have_left ||
                destination.x < left_primary ||
                (destination.x == left_primary &&
                 destination.y < left_secondary)) {
                have_left = true;
                left_primary = destination.x;
                left_secondary = destination.y;
                result.x.left_outer = sym(outer);
            }
            if (!have_right ||
                origin.x > right_primary ||
                (origin.x == right_primary &&
                 origin.y < right_secondary)) {
                have_right = true;
                right_primary = origin.x;
                right_secondary = origin.y;
                result.x.right_outer = outer;
            }
        } else {
            if (!have_left ||
                destination.y < left_primary ||
                (destination.y == left_primary &&
                 destination.x > left_secondary)) {
                have_left = true;
                left_primary = destination.y;
                left_secondary = destination.x;
                result.y.left_outer = sym(outer);
            }
            if (!have_right ||
                origin.y > right_primary ||
                (origin.y == right_primary &&
                 origin.x > right_secondary)) {
                have_right = true;
                right_primary = origin.y;
                right_secondary = origin.x;
                result.y.right_outer = outer;
            }
        }
        outer = lnext(outer);
    } while (outer != outer_seed);
    return result;
}

void Triangulator::acquire_edge_block(EdgeCursor& cursor) {
    finish_edge_cursor(cursor);
    if (cursor.block_counter == nullptr) {
        throw std::logic_error("parallel edge cursor has no block counter");
    }
    std::size_t capacity_limit = edge_capacity_limit_;
#if defined(DELAUNAY32_TEST_PARALLEL_EDGE_ARENA_DART_LIMIT)
    capacity_limit = std::min(
        capacity_limit,
        static_cast<std::size_t>(
            DELAUNAY32_TEST_PARALLEL_EDGE_ARENA_DART_LIMIT));
#endif
    const std::size_t first =
        cursor.block_counter->fetch_add(
            kEdgeBlockDarts, std::memory_order_relaxed);
    if (first >= capacity_limit) {
        throw detail::ParallelEdgeArenaExhausted{};
    }
    const std::size_t last =
        std::min(first + kEdgeBlockDarts, capacity_limit);
    cursor.next = static_cast<std::uint32_t>(first);
    cursor.end = static_cast<std::uint32_t>(last);
    cursor.range_first = cursor.next;
}

void Triangulator::finish_edge_cursor(EdgeCursor& cursor) {
    if (cursor.range_first < cursor.next) {
        cursor.ranges.push_back({cursor.range_first, cursor.next});
    }
    cursor.range_first = cursor.next;
}

template <bool ParallelAllocation>
std::uint32_t Triangulator::make_edge(
    std::uint32_t origin,
    std::uint32_t destination,
    EdgeCursor* cursor) {
    std::uint32_t edge = 0;
    if constexpr (ParallelAllocation) {
        if (cursor == nullptr) {
            throw std::logic_error("parallel edge allocation has no cursor");
        }
        if (cursor->next + 2 > cursor->end) {
            acquire_edge_block(*cursor);
        }
        edge = cursor->next;
        cursor->next += 2;
    } else {
        if (edge_count_ + 2 > edge_origin_.size()) {
            const std::size_t new_size =
                std::max(edge_origin_.size() * 2, edge_count_ + 2);
            edge_origin_.resize(new_size);
            edge_next_.resize(new_size);
            edge_prev_.resize(new_size);
        }
        edge = static_cast<std::uint32_t>(edge_count_);
        edge_count_ += 2;
    }
    edge_origin_[edge] = origin;
    edge_origin_[edge + 1] = destination;
    edge_next_[edge] = edge;
    edge_next_[edge + 1] = edge + 1;
    edge_prev_[edge] = edge;
    edge_prev_[edge + 1] = edge + 1;
    return edge;
}

void Triangulator::splice(std::uint32_t a, std::uint32_t b) {
    const std::uint32_t a_next = edge_next_[a];
    const std::uint32_t b_next = edge_next_[b];
    edge_next_[a] = b_next;
    edge_prev_[b_next] = a;
    edge_next_[b] = a_next;
    edge_prev_[a_next] = b;
}

template <bool ParallelAllocation>
std::uint32_t Triangulator::connect(
    std::uint32_t a,
    std::uint32_t b,
    EdgeCursor* cursor) {
    const std::uint32_t edge =
        make_edge<ParallelAllocation>(dest(a), org(b), cursor);
    splice(edge, lnext(a));
    splice(sym(edge), b);
    return edge;
}

void Triangulator::delete_edge(std::uint32_t edge) {
    splice(edge, oprev(edge));
    splice(sym(edge), oprev(sym(edge)));
    const std::uint32_t pair = edge & ~1U;
    edge_origin_[pair] = kDeletedEdge;
    edge_origin_[pair + 1] = kDeletedEdge;
}

std::int64_t Triangulator::orient(
    std::uint32_t a,
    std::uint32_t b,
    std::uint32_t c) const {
    const Site& pa = points_[a];
    const Site& pb = points_[b];
    const Site& pc = points_[c];
    return (static_cast<std::int64_t>(pb.x) - pa.x) *
               (static_cast<std::int64_t>(pc.y) - pa.y) -
           (static_cast<std::int64_t>(pb.y) - pa.y) *
               (static_cast<std::int64_t>(pc.x) - pa.x);
}

bool Triangulator::left_of(
    std::uint32_t point,
    std::uint32_t edge) const {
    return orient(org(edge), dest(edge), point) > 0;
}

bool Triangulator::right_of(
    std::uint32_t point,
    std::uint32_t edge) const {
    return orient(org(edge), dest(edge), point) < 0;
}

template <bool WidePredicates>
bool Triangulator::in_circle(
    std::uint32_t a,
    std::uint32_t b,
    std::uint32_t c,
    std::uint32_t d) const {
    const Site& pa = points_[a];
    const Site& pb = points_[b];
    const Site& pc = points_[c];
    const Site& pd = points_[d];
    if constexpr (WidePredicates) {
#if defined(__SIZEOF_INT128__)
        if (int64_wide_intermediates_) {
            // Inner products stay in int64; only the final three multiply
            // into the completed determinant need __int128.
            const std::int64_t ax =
                static_cast<std::int64_t>(pa.x) - pd.x;
            const std::int64_t ay =
                static_cast<std::int64_t>(pa.y) - pd.y;
            const std::int64_t bx =
                static_cast<std::int64_t>(pb.x) - pd.x;
            const std::int64_t by =
                static_cast<std::int64_t>(pb.y) - pd.y;
            const std::int64_t cx =
                static_cast<std::int64_t>(pc.x) - pd.x;
            const std::int64_t cy =
                static_cast<std::int64_t>(pc.y) - pd.y;
            const std::int64_t d_lift =
                static_cast<std::int64_t>(wide_lifts_[d]);
            const std::int64_t ap =
                static_cast<std::int64_t>(wide_lifts_[a]) - d_lift;
            const std::int64_t bp =
                static_cast<std::int64_t>(wide_lifts_[b]) - d_lift;
            const std::int64_t cp =
                static_cast<std::int64_t>(wide_lifts_[c]) - d_lift;
            const std::int64_t first = by * cp - bp * cy;
            const std::int64_t second = bx * cp - bp * cx;
            const std::int64_t third = bx * cy - by * cx;
            return static_cast<__int128_t>(ax) * first -
                       static_cast<__int128_t>(ay) * second +
                       static_cast<__int128_t>(ap) * third >
                   0;
        }
        const __int128_t ax = static_cast<__int128_t>(pa.x) - pd.x;
        const __int128_t ay = static_cast<__int128_t>(pa.y) - pd.y;
        const __int128_t bx = static_cast<__int128_t>(pb.x) - pd.x;
        const __int128_t by = static_cast<__int128_t>(pb.y) - pd.y;
        const __int128_t cx = static_cast<__int128_t>(pc.x) - pd.x;
        const __int128_t cy = static_cast<__int128_t>(pc.y) - pd.y;
        const __int128_t d_lift = wide_lifts_[d];
        const __int128_t ap =
            static_cast<__int128_t>(wide_lifts_[a]) - d_lift;
        const __int128_t bp =
            static_cast<__int128_t>(wide_lifts_[b]) - d_lift;
        const __int128_t cp =
            static_cast<__int128_t>(wide_lifts_[c]) - d_lift;
        return ax * (by * cp - bp * cy) -
                   ay * (bx * cp - bp * cx) +
                   ap * (bx * cy - by * cx) >
               0;
#else
        return false;
#endif
    } else {
        const std::int64_t ax = static_cast<std::int64_t>(pa.x) - pd.x;
        const std::int64_t ay = static_cast<std::int64_t>(pa.y) - pd.y;
        const std::int64_t bx = static_cast<std::int64_t>(pb.x) - pd.x;
        const std::int64_t by = static_cast<std::int64_t>(pb.y) - pd.y;
        const std::int64_t cx = static_cast<std::int64_t>(pc.x) - pd.x;
        const std::int64_t cy = static_cast<std::int64_t>(pc.y) - pd.y;
        const std::int64_t ap =
            static_cast<std::int64_t>(pa.lift) - pd.lift;
        const std::int64_t bp =
            static_cast<std::int64_t>(pb.lift) - pd.lift;
        const std::int64_t cp =
            static_cast<std::int64_t>(pc.lift) - pd.lift;
        return ax * (by * cp - bp * cy) -
                   ay * (bx * cp - bp * cx) +
                   ap * (bx * cy - by * cx) >
               0;
    }
}

// The API translation unit calls these specializations by declaration. Emit
// them here so all topology template expansion remains centralized.
template Triangulator::DirectionalHulls
Triangulator::build_morton_range<false, false>(
    std::size_t,
    std::size_t,
    EdgeCursor*);
template Triangulator::DirectionalHulls
Triangulator::build_morton_range<true, false>(
    std::size_t,
    std::size_t,
    EdgeCursor*);
template Triangulator::DirectionalHulls
Triangulator::build_parallel<false>(
    std::size_t,
    detail::WorkerTeam&);
template Triangulator::DirectionalHulls
Triangulator::build_parallel<true>(
    std::size_t,
    detail::WorkerTeam&);

}  // namespace delaunay32

#undef DELAUNAY32_ALWAYS_INLINE
