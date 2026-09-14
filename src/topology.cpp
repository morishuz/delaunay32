// SPDX-License-Identifier: MIT

#include "triangulator_impl.hpp"

#include <algorithm>
#include <atomic>
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
Triangulator::Impl::HullEdges Triangulator::Impl::build_range(
    std::size_t first,
    std::size_t last,
    EdgeCursor* cursor) {
    const std::size_t count = last - first;
    if (count == 2) {
        const std::uint32_t edge =
            arena_.make_edge<ParallelAllocation>(
            static_cast<std::uint32_t>(first),
            static_cast<std::uint32_t>(first + 1),
            cursor);
        return {edge, sym(edge)};
    }

    if (count == 3) {
        const std::uint32_t a =
            arena_.make_edge<ParallelAllocation>(
            static_cast<std::uint32_t>(first),
            static_cast<std::uint32_t>(first + 1),
            cursor);
        const std::uint32_t b =
            arena_.make_edge<ParallelAllocation>(
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

Triangulator::Impl::MortonSplit
Triangulator::Impl::find_morton_split(
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
Triangulator::Impl::DirectionalHulls
Triangulator::Impl::build_morton_range(
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
    return merge_directional_hulls<WidePredicates, ParallelAllocation>(
        left, right, horizontal, cursor);
}

template <
    bool WidePredicates,
    bool ParallelAllocation>
DELAUNAY32_ALWAYS_INLINE
Triangulator::Impl::HullEdges
Triangulator::Impl::merge_hulls_inline(
    HullEdges left,
    HullEdges right,
    EdgeCursor* cursor,
    OuterBridges* bridges) {
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
    if (bridges != nullptr) {
        bridges->first_outer = base;
    }
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
    if (bridges != nullptr) {
        bridges->last_outer = sym(base);
    }
    return {left.left_outer, right.right_outer};
}

template <
    bool WidePredicates,
    bool ParallelAllocation>
Triangulator::Impl::HullEdges
Triangulator::Impl::merge_hulls(
    HullEdges left,
    HullEdges right,
    EdgeCursor* cursor,
    OuterBridges* bridges) {
    // Leaf recursion expands the body through merge_hulls_inline(). Keeping
    // this wrapper out of other callers avoids duplicating the large kernel
    // throughout the Morton and parallel merge machinery.
    return merge_hulls_inline<WidePredicates, ParallelAllocation>(
        left, right, cursor, bridges);
}

template <
    bool WidePredicates,
    bool ParallelAllocation>
Triangulator::Impl::DirectionalHulls
Triangulator::Impl::merge_directional_hulls(
    const DirectionalHulls& left,
    const DirectionalHulls& right,
    bool horizontal,
    EdgeCursor* cursor) {
    // Extreme vertices of the union are extreme vertices of a child. Cache
    // their indices before merging, because discarded darts lose their origins.
    const HullEdges& a = horizontal ? left.x : left.y;
    const HullEdges& b = horizontal ? right.x : right.y;
    std::uint32_t left_vertex = org(a.left_outer);
    std::uint32_t right_vertex = org(a.right_outer);
    const std::uint32_t b_left_vertex = org(b.left_outer);
    const std::uint32_t b_right_vertex = org(b.right_outer);
    const Site& a_left = points_[left_vertex];
    const Site& a_right = points_[right_vertex];
    const Site& b_left = points_[b_left_vertex];
    const Site& b_right = points_[b_right_vertex];
    HullEdges perpendicular = a;
    if (horizontal) {
        if (b_left.x < a_left.x ||
            (b_left.x == a_left.x && b_left.y < a_left.y)) {
            perpendicular.left_outer = b.left_outer;
            left_vertex = b_left_vertex;
        }
        if (b_right.x > a_right.x ||
            (b_right.x == a_right.x && b_right.y < a_right.y)) {
            perpendicular.right_outer = b.right_outer;
            right_vertex = b_right_vertex;
        }
    } else {
        if (b_left.y < a_left.y ||
            (b_left.y == a_left.y && b_left.x > a_left.x)) {
            perpendicular.left_outer = b.left_outer;
            left_vertex = b_left_vertex;
        }
        if (b_right.y > a_right.y ||
            (b_right.y == a_right.y && b_right.x > a_right.x)) {
            perpendicular.right_outer = b.right_outer;
            right_vertex = b_right_vertex;
        }
    }

    OuterBridges bridges;
    const HullEdges merged =
        merge_hulls<WidePredicates, ParallelAllocation>(
            horizontal ? left.y : left.x,
            horizontal ? right.y : right.x,
            cursor,
            &bridges);
    // The merged boundary consists of retained child chains and these two
    // bridges. Only a bridge can replace an extreme vertex's boundary dart.
    const auto correct = [&](std::uint32_t bridge) {
        if (dest(bridge) == left_vertex) {
            perpendicular.left_outer = sym(bridge);
        }
        if (org(bridge) == right_vertex) {
            perpendicular.right_outer = bridge;
        }
    };
    correct(bridges.first_outer);
    correct(bridges.last_outer);

    DirectionalHulls result;
    if (horizontal) {
        result.x = perpendicular;
        result.y = merged;
    } else {
        result.x = merged;
        result.y = perpendicular;
    }
    return result;
}

std::size_t Triangulator::Impl::add_parallel_node(
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
Triangulator::Impl::DirectionalHulls
Triangulator::Impl::build_parallel(
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
        arena_.ranges.push_back(
            {0, static_cast<std::uint32_t>(arena_.count)});
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
    const std::size_t worker_count =
        std::min(thread_count, leaves.size());
    ThreadBarrier barrier(worker_count);
    const auto run_jobs = [&](std::size_t worker_index) {
        try {
            // Parallel topology proceeds bottom-up:
            //
            //   independent leaves -> barrier -> height-1 merges -> ...
            //
            // An abortable barrier retains the single worker-team run used by
            // the fast path while releasing peers if allocation throws.
            while (!cancelled.load(std::memory_order_relaxed)) {
                const std::size_t job =
                    next_job.fetch_add(1, std::memory_order_relaxed);
                if (job >= leaves.size()) {
                    break;
                }
                EdgeCursor& cursor = leaf_cursors[job];
                ParallelNode& node = nodes[leaves[job]];
                node.hull =
                    build_morton_range<WidePredicates, true>(
                        node.first,
                        node.last,
                        &cursor);
                detail::EdgeArena::finish_cursor(cursor);
            }
            if (!barrier.wait()) {
                return;
            }

            for (std::size_t level = 1;
                 level < merge_levels.size();
                 ++level) {
                if (worker_index == 0) {
                    next_job.store(0, std::memory_order_relaxed);
                }
                if (!barrier.wait()) {
                    return;
                }
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
                    const DirectionalHulls& left =
                        nodes[node.left].hull;
                    const DirectionalHulls& right =
                        nodes[node.right].hull;
                    const bool horizontal =
                        (node.split_bit & 1U) != 0;
                    node.hull =
                        merge_directional_hulls<WidePredicates, true>(
                            left, right, horizontal, &cursor);
                    detail::EdgeArena::finish_cursor(cursor);
                }
                if (!barrier.wait()) {
                    return;
                }
            }
        } catch (...) {
            cancelled.store(true, std::memory_order_relaxed);
            barrier.abort();
            throw;
        }
    };

    workers.run(worker_count, run_jobs);
    arena_.count = 0;
    const auto collect_ranges = [&](const EdgeCursor& cursor) {
        for (const EdgeRange range : cursor.ranges) {
            arena_.ranges.push_back(range);
            arena_.count += range.last - range.first;
        }
    };
    for (const EdgeCursor& cursor : leaf_cursors) {
        collect_ranges(cursor);
    }
    for (const EdgeCursor& cursor : node_cursors) {
        collect_ranges(cursor);
    }
    std::sort(
        arena_.ranges.begin(),
        arena_.ranges.end(),
        [](const EdgeRange& a, const EdgeRange& b) {
            return a.first < b.first;
        });
    return nodes.front().hull;
}

Triangulator::Impl::DirectionalHulls
Triangulator::Impl::scan_directional_hulls(
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

void Triangulator::Impl::splice(std::uint32_t a, std::uint32_t b) {
    const std::uint32_t a_next = arena_.next[a];
    const std::uint32_t b_next = arena_.next[b];
    arena_.next[a] = b_next;
    arena_.prev[b_next] = a;
    arena_.next[b] = a_next;
    arena_.prev[a_next] = b;
}

template <bool ParallelAllocation>
std::uint32_t Triangulator::Impl::connect(
    std::uint32_t a,
    std::uint32_t b,
    EdgeCursor* cursor) {
    const std::uint32_t edge =
        arena_.make_edge<ParallelAllocation>(dest(a), org(b), cursor);
    splice(edge, lnext(a));
    splice(sym(edge), b);
    return edge;
}

void Triangulator::Impl::delete_edge(std::uint32_t edge) {
    splice(edge, oprev(edge));
    splice(sym(edge), oprev(sym(edge)));
    const std::uint32_t pair = edge & ~1U;
    arena_.origin[pair] = kDeletedEdge;
    arena_.origin[pair + 1] = kDeletedEdge;
}

bool Triangulator::Impl::is_live_edge(std::uint32_t edge) const {
    const std::uint32_t pair = edge & ~1U;
    return pair + 1U < arena_.origin.size() &&
           arena_.origin[pair] != kDeletedEdge;
}

bool Triangulator::Impl::is_constrained(std::uint32_t edge) const {
    return (constraints_.edge_flags[edge & ~1U] & kConstrainedBit) != 0;
}

void Triangulator::Impl::mark_constrained(std::uint32_t edge) {
    const std::uint32_t pair = edge & ~1U;
    constraints_.edge_flags[pair] |= kConstrainedBit;
    constraints_.edge_flags[pair + 1U] |= kConstrainedBit;
}

bool Triangulator::Impl::left_face_opposite(
    std::uint32_t edge,
    std::uint32_t& opposite) const {
    if (!is_live_edge(edge)) {
        return false;
    }
    const std::uint32_t second = lnext(edge);
    const std::uint32_t third = lnext(second);
    if (lnext(third) != edge) {
        return false;
    }
    opposite = dest(second);
    return orient(org(edge), dest(edge), opposite) > 0;
}

bool Triangulator::Impl::can_flip(std::uint32_t edge) const {
    if (!is_live_edge(edge) || is_constrained(edge)) {
        return false;
    }
    std::uint32_t left = 0;
    std::uint32_t right = 0;
    if (!left_face_opposite(edge, left) ||
        !left_face_opposite(sym(edge), right)) {
        return false;
    }
    const std::int64_t origin_side =
        orient(left, right, org(edge));
    const std::int64_t destination_side =
        orient(left, right, dest(edge));
    return (origin_side > 0 && destination_side < 0) ||
           (origin_side < 0 && destination_side > 0);
}

void Triangulator::Impl::flip_edge(std::uint32_t edge) {
    if (!can_flip(edge)) {
        throw std::logic_error("attempted to flip a nonflippable edge");
    }
    const std::uint32_t old_origin = org(edge);
    const std::uint32_t old_destination = dest(edge);
    const std::uint32_t origin_previous = oprev(edge);
    const std::uint32_t destination_previous = oprev(sym(edge));

    // Guibas-Stolfi diagonal swap. Reusing the dart pair avoids growing the
    // append-only construction arena during constraint recovery.
    splice(edge, origin_previous);
    splice(sym(edge), destination_previous);
    splice(edge, lnext(origin_previous));
    splice(sym(edge), lnext(destination_previous));
    arena_.origin[edge] = dest(origin_previous);
    arena_.origin[sym(edge)] = dest(destination_previous);

    if (!constraints_.site_edge.empty()) {
        constraints_.site_edge[old_origin] = origin_previous;
        constraints_.site_edge[old_destination] = destination_previous;
        constraints_.site_edge[org(edge)] = edge;
        constraints_.site_edge[dest(edge)] = sym(edge);
    }
}

std::int64_t Triangulator::Impl::orient(
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

bool Triangulator::Impl::left_of(
    std::uint32_t point,
    std::uint32_t edge) const {
    return orient(org(edge), dest(edge), point) > 0;
}

bool Triangulator::Impl::right_of(
    std::uint32_t point,
    std::uint32_t edge) const {
    return orient(org(edge), dest(edge), point) < 0;
}

bool Triangulator::Impl::active_in_circle(
    std::uint32_t a,
    std::uint32_t b,
    std::uint32_t c,
    std::uint32_t d) const {
    return wide_lifts_.empty()
               ? in_circle<false>(a, b, c, d)
               : in_circle<true>(a, b, c, d);
}

template <bool WidePredicates>
bool Triangulator::Impl::in_circle(
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
template Triangulator::Impl::DirectionalHulls
Triangulator::Impl::build_morton_range<false, false>(
    std::size_t,
    std::size_t,
    EdgeCursor*);
template Triangulator::Impl::DirectionalHulls
Triangulator::Impl::build_morton_range<true, false>(
    std::size_t,
    std::size_t,
    EdgeCursor*);
template Triangulator::Impl::DirectionalHulls
Triangulator::Impl::build_parallel<false>(
    std::size_t,
    detail::WorkerTeam&);
template Triangulator::Impl::DirectionalHulls
Triangulator::Impl::build_parallel<true>(
    std::size_t,
    detail::WorkerTeam&);

}  // namespace delaunay32

#undef DELAUNAY32_ALWAYS_INLINE
