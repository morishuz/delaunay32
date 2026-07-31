// SPDX-License-Identifier: MIT

#include "delaunay32/delaunay.hpp"
#include "internal.hpp"

#include <algorithm>
#include <limits>
#include <memory>
#include <stdexcept>
#include <thread>
#include <utility>
#include <vector>

namespace delaunay32 {
using detail::ThreadBarrier;

// Public entry points, input normalization, and Morton ordering live here.

Triangulator::Triangulator(std::size_t thread_count)
    : thread_count_(thread_count) {}

Triangulator::~Triangulator() = default;

Triangulator::Triangulator(
    Triangulator&&) noexcept = default;

Triangulator& Triangulator::operator=(
    Triangulator&&) noexcept = default;

namespace {

// Shared predicate-span checks for Int64 vs Int128 certification.
template <typename UnsignedWide>
bool product_fits(
    UnsignedWide lift_bound,
    UnsignedWide factor,
    UnsignedWide limit) noexcept {
    return factor == 0 || lift_bound <= limit / factor;
}

template <typename UnsignedWide>
bool predicates_fit_limit(
    UnsignedWide sx,
    UnsignedWide sy,
    UnsignedWide limit) noexcept {
    const UnsignedWide lift_bound = sx * sx + sy * sy;
    const UnsignedWide intermediate_factor = 2U * std::max(sx, sy);
    const UnsignedWide determinant_factor = 6U * sx * sy;
    return product_fits(lift_bound, intermediate_factor, limit) &&
           product_fits(lift_bound, determinant_factor, limit);
}

}  // namespace

void Triangulator::require_point_count(std::size_t point_count) {
    if (point_count < 3) {
        throw std::invalid_argument("need at least 3 points");
    }
    if (point_count > kIndexMask) {
        throw std::invalid_argument("point count exceeds internal index range");
    }
}

detail::WorkerTeam* Triangulator::ensure_worker_team(
    std::size_t thread_count) {
    if (thread_count <= 1) {
        return nullptr;
    }
    if (worker_team_size_ != thread_count || worker_team_ == nullptr) {
        worker_team_ = std::make_unique<detail::WorkerTeam>(thread_count);
        worker_team_size_ = thread_count;
        // Keep per-worker export buffers warm across multi-shot clients.
        export_scratch_.assign(thread_count, {});
    }
    return worker_team_.get();
}

PredicateWidth
Triangulator::predicate_width_for_spans(
    std::uint64_t x_span,
    std::uint64_t y_span) noexcept {
    if (x_span > std::numeric_limits<std::uint32_t>::max() ||
        y_span > std::numeric_limits<std::uint32_t>::max()) {
        return PredicateWidth::Unsupported;
    }

#if defined(__SIZEOF_INT128__)
    using UnsignedWide = __uint128_t;
    const UnsignedWide sx = x_span;
    const UnsignedWide sy = y_span;
    const UnsignedWide lift_bound = sx * sx + sy * sy;
    const UnsignedWide orientation_bound = 2U * sx * sy;
    // Inner incircle products and the completed determinant have different
    // worst cases for thin domains, so both bounds must fit.
    const UnsignedWide signed64 =
        static_cast<UnsignedWide>(std::numeric_limits<std::int64_t>::max());

    if (lift_bound <= std::numeric_limits<std::uint32_t>::max() &&
        orientation_bound <= signed64 &&
        predicates_fit_limit(sx, sy, signed64)) {
        return PredicateWidth::Int64;
    }

    constexpr UnsignedWide kSigned128Max =
        (UnsignedWide{1} << 127U) - 1U;
    if (lift_bound <= std::numeric_limits<std::uint64_t>::max() &&
        orientation_bound <= signed64 &&
        predicates_fit_limit(sx, sy, kSigned128Max)) {
        return PredicateWidth::Int128;
    }
#else
    if (x_span <= 65535U && y_span <= 65535U) {
        const std::uint64_t lift_bound = x_span * x_span + y_span * y_span;
        const std::uint64_t orientation_bound = 2U * x_span * y_span;
        const std::uint64_t signed_limit =
            static_cast<std::uint64_t>(
                std::numeric_limits<std::int64_t>::max());
        if (lift_bound <= std::numeric_limits<std::uint32_t>::max() &&
            orientation_bound <= signed_limit &&
            predicates_fit_limit(x_span, y_span, signed_limit)) {
            return PredicateWidth::Int64;
        }
    }
#endif
    return PredicateWidth::Unsupported;
}

bool Triangulator::int64_wide_intermediates_for_spans(
    std::uint64_t x_span,
    std::uint64_t y_span) noexcept {
#if defined(__SIZEOF_INT128__)
    if (x_span > std::numeric_limits<std::uint32_t>::max() ||
        y_span > std::numeric_limits<std::uint32_t>::max()) {
        return false;
    }
    using UnsignedWide = __uint128_t;
    const UnsignedWide sx = x_span;
    const UnsignedWide sy = y_span;
    const UnsignedWide lift_bound = sx * sx + sy * sy;
    const UnsignedWide factor = 2U * std::max(sx, sy);
    const UnsignedWide limit =
        static_cast<UnsignedWide>(
            std::numeric_limits<std::int64_t>::max());
    return product_fits(lift_bound, factor, limit);
#else
    (void)x_span;
    (void)y_span;
    return false;
#endif
}

std::vector<Triangle> Triangulator::triangulate(
    const std::vector<Point>& points) {
    return triangulate_points_impl(points);
}

std::vector<Triangle> Triangulator::triangulate_points_impl(
    const std::vector<Point>& points) {
    require_point_count(points.size());
    points_.resize(points.size());
    min_x_ = max_x_ = points[0].x;
    min_y_ = max_y_ = points[0].y;
    points_[0] = {points[0].x, points[0].y, 0, 0};
    for (std::size_t i = 1; i < points.size(); ++i) {
        min_x_ = std::min(min_x_, points[i].x);
        max_x_ = std::max(max_x_, points[i].x);
        min_y_ = std::min(min_y_, points[i].y);
        max_y_ = std::max(max_y_, points[i].y);
        points_[i] = {
            points[i].x,
            points[i].y,
            static_cast<std::uint32_t>(i),
            0,
        };
    }
    prepare_points(points.size());
    return std::move(triangles_out_);
}

std::vector<Triangle> Triangulator::triangulate_int(
    const std::int32_t* xs,
    const std::int32_t* ys,
    std::size_t point_count) {
    return triangulate_int_impl(xs, ys, point_count);
}

std::vector<Triangle>
Triangulator::triangulate_int_impl(
    const std::int32_t* xs,
    const std::int32_t* ys,
    std::size_t point_count) {
    require_point_count(point_count);
    if (xs == nullptr || ys == nullptr) {
        throw std::invalid_argument("coordinate arrays must not be null");
    }
    points_.resize(point_count);
    min_x_ = max_x_ = xs[0];
    min_y_ = max_y_ = ys[0];
    points_[0] = {xs[0], ys[0], 0, 0};
    for (std::size_t i = 1; i < point_count; ++i) {
        min_x_ = std::min(min_x_, xs[i]);
        max_x_ = std::max(max_x_, xs[i]);
        min_y_ = std::min(min_y_, ys[i]);
        max_y_ = std::max(max_y_, ys[i]);
        points_[i] = {xs[i], ys[i], static_cast<std::uint32_t>(i), 0};
    }
    prepare_points(point_count);
    return std::move(triangles_out_);
}

void Triangulator::prepare_points(std::size_t point_count) {
    const std::int64_t x_span = static_cast<std::int64_t>(max_x_) - min_x_;
    const std::int64_t y_span = static_cast<std::int64_t>(max_y_) - min_y_;
    const PredicateWidth predicate_width =
        predicate_width_for_spans(
            static_cast<std::uint64_t>(x_span),
            static_cast<std::uint64_t>(y_span));
    if (predicate_width == PredicateWidth::Unsupported) {
        throw std::invalid_argument(
            "coordinate spans exceed exact predicate range");
    }
    const bool wide_predicates =
        predicate_width == PredicateWidth::Int128;
    int64_wide_intermediates_ =
        wide_predicates &&
        int64_wide_intermediates_for_spans(
            static_cast<std::uint64_t>(x_span),
            static_cast<std::uint64_t>(y_span));
    const std::uint64_t maximum_span =
        static_cast<std::uint64_t>(std::max(x_span, y_span));
    std::uint32_t morton_shift = 0;
    while ((maximum_span >> morton_shift) > 65535U) {
        ++morton_shift;
    }
    morton_keys_.resize(points_.size());
    if (wide_predicates) {
        wide_lifts_.resize(points_.size());
    } else {
        wide_lifts_.clear();
    }
    std::uint32_t maximum_morton_key = 0;
    for (std::size_t i = 0; i < points_.size(); ++i) {
        Site& point = points_[i];
        const std::uint32_t x = static_cast<std::uint32_t>(
            static_cast<std::int64_t>(point.x) - min_x_);
        const std::uint32_t y = static_cast<std::uint32_t>(
            static_cast<std::int64_t>(point.y) - min_y_);
        if (!wide_predicates) {
            point.lift = x * x + y * y;
        }
        morton_keys_[i] = morton_code(x >> morton_shift, y >> morton_shift);
        maximum_morton_key = std::max(maximum_morton_key, morton_keys_[i]);
    }
    std::size_t effective_threads = thread_count_;
    if (effective_threads == 0) {
        effective_threads =
            std::max<std::size_t>(1, std::thread::hardware_concurrency());
    }
    effective_threads =
        std::min(effective_threads, kMaxParallelThreads);
    if (point_count < kParallelMinPoints) {
        effective_threads = 1;
    }
    detail::WorkerTeam* workers =
        effective_threads > 1
            ? ensure_worker_team(effective_threads)
            : nullptr;
    sort_points_morton(effective_threads, maximum_morton_key, workers);
    // Equal coordinates share a Morton key. Keep the unique-input scan
    // write-free; SiteLessXY puts the lowest original index first, and only
    // the suffix after the first duplicate needs in-place compaction.
    std::size_t duplicate_group = points_.size();
    for (std::size_t first = 0; first < points_.size();) {
        std::size_t last = first + 1;
        while (last < points_.size() &&
               morton_keys_[last] == morton_keys_[first]) {
            ++last;
        }
        if (last - first > 1) {
            std::sort(
                points_.begin() + static_cast<std::ptrdiff_t>(first),
                points_.begin() + static_cast<std::ptrdiff_t>(last),
                SiteLessXY{});
            for (std::size_t i = first + 1; i < last; ++i) {
                if (points_[i - 1].x == points_[i].x &&
                    points_[i - 1].y == points_[i].y) {
                    duplicate_group = first;
                    break;
                }
            }
            if (duplicate_group != points_.size()) {
                break;
            }
        }
        first = last;
    }

    std::size_t unique_count = points_.size();
    if (duplicate_group != points_.size()) {
        unique_count = duplicate_group;
        bool group_already_sorted = true;
        for (std::size_t first = duplicate_group;
             first < points_.size();) {
            std::size_t last = first + 1;
            while (last < points_.size() &&
                   morton_keys_[last] == morton_keys_[first]) {
                ++last;
            }
            if (!group_already_sorted && last - first > 1) {
                std::sort(
                    points_.begin() + static_cast<std::ptrdiff_t>(first),
                    points_.begin() + static_cast<std::ptrdiff_t>(last),
                    SiteLessXY{});
            }

            bool have_coordinate = false;
            std::int32_t previous_x = 0;
            std::int32_t previous_y = 0;
            const std::uint32_t key = morton_keys_[first];
            for (std::size_t i = first; i < last; ++i) {
                const Site point = points_[i];
                if (have_coordinate &&
                    point.x == previous_x &&
                    point.y == previous_y) {
                    continue;
                }
                have_coordinate = true;
                previous_x = point.x;
                previous_y = point.y;
                points_[unique_count] = point;
                morton_keys_[unique_count] = key;
                ++unique_count;
            }
            group_already_sorted = false;
            first = last;
        }
        points_.resize(unique_count);
        morton_keys_.resize(unique_count);
        if (wide_predicates) {
            wide_lifts_.resize(unique_count);
        }
    }
    if (unique_count < 3) {
        throw std::invalid_argument("need at least 3 unique points");
    }
    point_count = unique_count;
    if (point_count < kParallelMinPoints) {
        effective_threads = 1;
    }
    edge_count_ = 0;
    edge_ranges_.clear();
    if (point_count >
        static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max()) /
            kEdgeArenaDartsPerPoint) {
        throw std::invalid_argument(
            "point count exceeds edge arena index range");
    }
    const std::size_t edge_capacity =
        point_count * kEdgeArenaDartsPerPoint;
    edge_capacity_limit_ = edge_capacity;
    if (edge_origin_.size() < edge_capacity) {
        edge_origin_.resize(edge_capacity);
        edge_next_.resize(edge_capacity);
        edge_prev_.resize(edge_capacity);
    }
    active_thread_count_ = effective_threads;
    const DirectionalHulls hull =
        effective_threads > 1
            ? (wide_predicates
                   ? build_parallel<true>(
                         effective_threads, *workers)
                   : build_parallel<false>(
                         effective_threads, *workers))
            : (wide_predicates
                   ? build_morton_range<true>(0, point_count)
                   : build_morton_range<false>(0, point_count));
    if (effective_threads == 1) {
        edge_ranges_.push_back(
            {0, static_cast<std::uint32_t>(edge_count_)});
    }
    outer_seed_ = sym(hull.x.left_outer);
    if (active_thread_count_ > 1) {
        export_triangles_parallel(active_thread_count_, *workers);
    } else {
        export_triangles();
    }
}

void Triangulator::sort_points_morton(
    std::size_t thread_count,
    std::uint32_t maximum_key,
    detail::WorkerTeam* workers) {
    static constexpr unsigned kRadixBits = 10;
    static constexpr std::uint32_t kRadixSize = 1U << kRadixBits;
    static constexpr std::uint32_t kRadixMask = kRadixSize - 1U;

    sort_scratch_.resize(points_.size());
    morton_scratch_.resize(points_.size());

    std::size_t pass_count = 1;
    while (pass_count < 4 &&
           (maximum_key >> (pass_count * kRadixBits)) != 0) {
        ++pass_count;
    }

    if (thread_count > 1 &&
        points_.size() >= kParallelRadixMinPoints) {
        // One stable radix pass:
        //
        //   private histograms -> barrier -> worker 0 assigns bucket offsets
        //                      -> barrier -> disjoint stable scatter -> barrier
        //
        // Workers own contiguous input slices. Assigning offsets in worker
        // order therefore preserves the input order within every bucket.
        sort_counts_.assign(thread_count * kRadixSize, 0);
        ThreadBarrier barrier(thread_count);
        const auto run = [&](std::size_t worker_index) {
            const std::size_t first =
                points_.size() * worker_index / thread_count;
            const std::size_t last =
                points_.size() * (worker_index + 1) / thread_count;
            std::uint32_t* counts =
                sort_counts_.data() + worker_index * kRadixSize;

            for (std::size_t pass = 0; pass < pass_count; ++pass) {
                std::fill(counts, counts + kRadixSize, 0);
                const unsigned shift =
                    static_cast<unsigned>(pass * kRadixBits);
                const std::vector<Site>& input_points =
                    (pass & 1U) == 0 ? points_ : sort_scratch_;
                const std::vector<std::uint32_t>& input_keys =
                    (pass & 1U) == 0 ? morton_keys_ : morton_scratch_;
                std::vector<Site>& output_points =
                    (pass & 1U) == 0 ? sort_scratch_ : points_;
                std::vector<std::uint32_t>& output_keys =
                    (pass & 1U) == 0 ? morton_scratch_ : morton_keys_;

                for (std::size_t i = first; i < last; ++i) {
                    ++counts[(input_keys[i] >> shift) & kRadixMask];
                }
                barrier.wait();
                if (worker_index == 0) {
                    std::uint32_t offset = 0;
                    for (std::size_t bucket = 0;
                         bucket < kRadixSize;
                         ++bucket) {
                        for (std::size_t worker = 0;
                             worker < thread_count;
                             ++worker) {
                            std::uint32_t& count =
                                sort_counts_[worker * kRadixSize + bucket];
                            const std::uint32_t next = offset + count;
                            count = offset;
                            offset = next;
                        }
                    }
                }
                barrier.wait();
                for (std::size_t i = first; i < last; ++i) {
                    const std::uint32_t key = input_keys[i];
                    const std::uint32_t destination =
                        counts[(key >> shift) & kRadixMask]++;
                    output_points[destination] = input_points[i];
                    output_keys[destination] = key;
                }
                barrier.wait();
            }
        };

        if (workers == nullptr) {
            throw std::logic_error(
                "parallel radix sort requires a worker team");
        }
        workers->run(thread_count, run);
        if ((pass_count & 1U) != 0) {
            points_.swap(sort_scratch_);
            morton_keys_.swap(morton_scratch_);
        }
        return;
    }

    const auto scatter = [&](const std::vector<Site>& input_points,
                             const std::vector<std::uint32_t>& input_keys,
                             std::vector<Site>& output_points,
                             std::vector<std::uint32_t>& output_keys,
                             unsigned shift) {
        sort_counts_.assign(kRadixSize, 0);
        for (const std::uint32_t key : input_keys) {
            ++sort_counts_[(key >> shift) & kRadixMask];
        }
        std::uint32_t offset = 0;
        for (std::uint32_t& count : sort_counts_) {
            const std::uint32_t next = offset + count;
            count = offset;
            offset = next;
        }
        for (std::size_t i = 0; i < input_points.size(); ++i) {
            const std::uint32_t key = input_keys[i];
            const std::uint32_t destination =
                sort_counts_[(key >> shift) & kRadixMask]++;
            output_points[destination] = input_points[i];
            output_keys[destination] = key;
        }
    };

    sort_counts_.assign(kRadixSize, 0);
    for (const std::uint32_t key : morton_keys_) {
        ++sort_counts_[key & kRadixMask];
    }
    std::uint32_t offset = 0;
    for (std::uint32_t& count : sort_counts_) {
        const std::uint32_t next = offset + count;
        count = offset;
        offset = next;
    }
    for (std::size_t i = 0; i < points_.size(); ++i) {
        const std::uint32_t key = morton_keys_[i];
        const std::uint32_t destination =
            sort_counts_[key & kRadixMask]++;
        sort_scratch_[destination] = points_[i];
        morton_scratch_[destination] = key;
    }
    if (maximum_key <= kRadixMask) {
        points_.swap(sort_scratch_);
        morton_keys_.swap(morton_scratch_);
        return;
    }

    scatter(
        sort_scratch_,
        morton_scratch_,
        points_,
        morton_keys_,
        kRadixBits);
    if (maximum_key <= ((1U << (2U * kRadixBits)) - 1U)) {
        return;
    }

    scatter(
        points_,
        morton_keys_,
        sort_scratch_,
        morton_scratch_,
        2U * kRadixBits);
    if (maximum_key <= ((1U << (3U * kRadixBits)) - 1U)) {
        points_.swap(sort_scratch_);
        morton_keys_.swap(morton_scratch_);
        return;
    }

    scatter(
        sort_scratch_,
        morton_scratch_,
        points_,
        morton_keys_,
        3U * kRadixBits);
}

std::uint32_t Triangulator::morton_code(
    std::uint32_t x,
    std::uint32_t y) {
    const auto spread = [](std::uint32_t value) {
        value &= 0x0000ffffU;
        value = (value | (value << 8U)) & 0x00ff00ffU;
        value = (value | (value << 4U)) & 0x0f0f0f0fU;
        value = (value | (value << 2U)) & 0x33333333U;
        value = (value | (value << 1U)) & 0x55555555U;
        return value;
    };
    return spread(x) | (spread(y) << 1U);
}

}  // namespace delaunay32
