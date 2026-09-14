// SPDX-License-Identifier: MIT

#include "triangulator_impl.hpp"

#include <algorithm>
#include <limits>
#include <memory>
#include <stdexcept>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

namespace delaunay32 {
using detail::ThreadBarrier;

// Public entry points, input normalization, and Morton ordering live here.

Triangulator::Triangulator() = default;
Triangulator::~Triangulator() = default;
Triangulator::Triangulator(Triangulator&&) noexcept = default;
Triangulator& Triangulator::operator=(Triangulator&&) noexcept = default;

Triangulator::Impl& Triangulator::implementation() {
    // Lazy ownership also lets a moved-from instance start another problem.
    if (!impl_) {
        impl_ = std::make_unique<Impl>();
    }
    return *impl_;
}

void Triangulator::set_options(TriangulationOptions options) {
    implementation().set_options(options);
}

void Triangulator::set_points(const std::vector<Point>& points) {
    implementation().set_points(points);
}

void Triangulator::set_constraints(std::vector<Constraint> constraints) {
    implementation().set_constraints(std::move(constraints));
}

void Triangulator::set_polygons(std::vector<PolygonDomain> polygons) {
    implementation().set_polygons(std::move(polygons));
}

TriangulationResult Triangulator::triangulate() {
    return implementation().triangulate();
}

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

void Triangulator::Impl::require_point_count(std::size_t point_count) {
    if (point_count < 3) {
        throw std::invalid_argument("need at least 3 points");
    }
    if (point_count > kIndexMask) {
        throw std::invalid_argument("point count exceeds internal index range");
    }
}

detail::WorkerTeam* Triangulator::Impl::ensure_worker_team(
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

void Triangulator::Impl::load_int_points(const std::vector<Point>& points) {
    points_.resize(points.size());
    min_x_ = max_x_ = points[0].x;
    min_y_ = max_y_ = points[0].y;
    points_[0] = {min_x_, min_y_, 0, 0};
    for (std::size_t i = 1; i < points.size(); ++i) {
        const std::int32_t x = points[i].x;
        const std::int32_t y = points[i].y;
        min_x_ = std::min(min_x_, x);
        max_x_ = std::max(max_x_, x);
        min_y_ = std::min(min_y_, y);
        max_y_ = std::max(max_y_, y);
        points_[i] = {
            x,
            y,
            static_cast<std::uint32_t>(i),
            0,
        };
    }
}

void Triangulator::Impl::set_options(TriangulationOptions options) {
    if (options.result_detail != ResultDetail::Triangles &&
        options.result_detail != ResultDetail::Full) {
        throw std::invalid_argument("unknown result detail");
    }
    thread_count_ = options.thread_count;
    result_detail_ = options.result_detail;
}

void Triangulator::Impl::set_points(const std::vector<Point>& points) {
    problem_ready_ = false;
    constraints_.clear_problem();
    triangles_out_.clear();
    halfedges_out_.clear();
    hull_out_.clear();
    require_point_count(points.size());
    input_point_count_ = points.size();
    load_int_points(points);
    problem_ready_ = true;
}

void Triangulator::Impl::require_ready_problem() const {
    if (!problem_ready_) {
        throw std::logic_error(
            "set_points() is required before configuring or triangulating");
    }
}

void Triangulator::Impl::set_constraints(
    std::vector<Constraint> constraints) {
    require_ready_problem();
    for (const Constraint constraint : constraints) {
        if (constraint.i0 >= input_point_count_ ||
            constraint.i1 >= input_point_count_) {
            throw std::invalid_argument(
                "constraint endpoint is outside the point array");
        }
        if (constraint.i0 == constraint.i1) {
            throw std::invalid_argument(
                "constraint endpoints are coincident");
        }
    }
    constraints_.segments = std::move(constraints);
}

void Triangulator::Impl::set_polygons(
    std::vector<PolygonDomain> polygons) {
    require_ready_problem();
    const auto validate_ring = [&](const std::vector<std::uint32_t>& ring) {
        if (ring.size() < 3) {
            throw std::invalid_argument(
                "polygon rings need at least three indices");
        }
        for (const std::uint32_t index : ring) {
            if (index >= input_point_count_) {
                throw std::invalid_argument(
                    "polygon ring index is outside the point array");
            }
        }
    };
    for (const PolygonDomain& polygon : polygons) {
        validate_ring(polygon.outer_ring);
        for (const std::vector<std::uint32_t>& hole : polygon.holes) {
            validate_ring(hole);
        }
    }
    constraints_.polygons = std::move(polygons);
}

TriangulationResult Triangulator::Impl::triangulate() {
    require_ready_problem();
    // A failed run is consumed as well: topology construction and constraint
    // recovery mutate the loaded sites and edge arena in place.
    problem_ready_ = false;

    const bool need_representatives =
        result_detail_ == ResultDetail::Full ||
        !constraints_.segments.empty() || !constraints_.polygons.empty();
    std::vector<std::uint32_t> representatives;
    const PredicateWidth predicate_width = build_loaded_topology(
        need_representatives ? &representatives : nullptr);

    if (!constraints_.segments.empty() || !constraints_.polygons.empty()) {
        build_constraint_indices(representatives, input_point_count_);
        const auto domains = prepare_polygon_domains();
        std::vector<Constraint> combined = constraints_.segments;
        std::size_t boundary_count = 0;
        for (const auto& domain : domains) {
            for (const auto& ring : domain) {
                boundary_count += ring.size();
            }
        }
        combined.reserve(combined.size() + boundary_count);
        for (const auto& domain : domains) {
            for (const auto& ring : domain) {
                for (std::size_t i = 0; i < ring.size(); ++i) {
                    combined.push_back({
                        points_[ring[i]].original,
                        points_[ring[(i + 1) % ring.size()]].original,
                    });
                }
            }
        }
        recover_constraints(combined);
        if (!domains.empty()) {
            mark_polygon_excluded_faces(domains);
        }
    }

    if (result_detail_ == ResultDetail::Full) {
        finish_full_export();
    } else {
        finish_triangle_export();
        representatives.clear();
        halfedges_out_.clear();
        hull_out_.clear();
    }
    return make_result(predicate_width, std::move(representatives));
}

PredicateWidth Triangulator::Impl::build_loaded_topology(
    std::vector<std::uint32_t>* representatives) {
    std::size_t point_count = points_.size();
    if (representatives != nullptr) {
        representatives->resize(point_count);
        for (std::size_t i = 0; i < point_count; ++i) {
            (*representatives)[i] = static_cast<std::uint32_t>(i);
        }
    }
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
            std::uint32_t representative = 0;
            const std::uint32_t key = morton_keys_[first];
            for (std::size_t i = first; i < last; ++i) {
                const Site point = points_[i];
                if (have_coordinate &&
                    point.x == previous_x &&
                    point.y == previous_y) {
                    if (representatives != nullptr) {
                        (*representatives)[point.original] = representative;
                    }
                    continue;
                }
                have_coordinate = true;
                previous_x = point.x;
                previous_y = point.y;
                representative = point.original;
                if (representatives != nullptr) {
                    (*representatives)[point.original] = representative;
                }
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
    arena_.count = 0;
    arena_.ranges.clear();
    if (point_count >
        static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max()) /
            detail::EdgeArena::kDartsPerPoint) {
        throw std::invalid_argument(
            "point count exceeds edge arena index range");
    }
    std::size_t edge_capacity =
        point_count * detail::EdgeArena::kDartsPerPoint;
#if defined(DELAUNAY32_TEST_PARALLEL_EDGE_ARENA_DART_LIMIT)
    // Match the physical allocation to the private test limit so sanitizer
    // runs also detect a write past a partially filled final edge block.
    edge_capacity = std::min(
        edge_capacity,
        static_cast<std::size_t>(
            DELAUNAY32_TEST_PARALLEL_EDGE_ARENA_DART_LIMIT));
#endif
    arena_.capacity_limit = edge_capacity;
    if (arena_.origin.size() < edge_capacity) {
        arena_.resize(edge_capacity);
    }
    active_thread_count_ = effective_threads;
    DirectionalHulls hull;
    bool serial_build = effective_threads == 1;
    if (!serial_build) {
        try {
            hull = wide_predicates
                       ? build_parallel<true>(
                             effective_threads, *workers)
                       : build_parallel<false>(
                             effective_threads, *workers);
        } catch (const detail::ParallelEdgeArenaExhausted&) {
            // The parallel arena deliberately cannot grow while workers hold
            // indices into it. Discard the partial topology and rebuild using
            // the serial allocator, which grows its arrays as needed. Leaf
            // sorting and wide lifts are recomputed during the serial build.
            arena_.count = 0;
            arena_.ranges.clear();
            active_thread_count_ = 1;
            serial_build = true;
        }
    }
    if (serial_build) {
        hull = wide_predicates
                   ? build_morton_range<true>(0, point_count)
                   : build_morton_range<false>(0, point_count);
        arena_.ranges.push_back(
            {0, static_cast<std::uint32_t>(arena_.count)});
    }
    outer_seed_ = sym(hull.x.left_outer);
    return predicate_width;
}

void Triangulator::Impl::sort_points_morton(
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

    const auto sort_serial = [&](auto radix_bits) {
        static constexpr unsigned kSerialRadixBits =
            decltype(radix_bits)::value;
        static constexpr std::uint32_t kSerialRadixSize =
            1U << kSerialRadixBits;
        static constexpr std::uint32_t kSerialRadixMask =
            kSerialRadixSize - 1U;
        static constexpr std::uint32_t kSerialLastRadixSize =
            1U << (32U - 2U * kSerialRadixBits);

        const unsigned serial_pass_count =
            maximum_key <= kSerialRadixMask
                ? 1U
                : maximum_key < (1U << (2U * kSerialRadixBits)) ? 2U : 3U;
        const std::size_t histogram_size =
            serial_pass_count == 3
                ? 2U * kSerialRadixSize + kSerialLastRadixSize
                : serial_pass_count * kSerialRadixSize;
        sort_counts_.assign(histogram_size, 0);

        // Digit frequencies do not change when points are permuted. Count all
        // active digits in one scan. Each stable scatter moves the point and
        // its Morton key together.
        if (serial_pass_count == 1) {
            for (const std::uint32_t key : morton_keys_) {
                ++sort_counts_[key & kSerialRadixMask];
            }
        } else if (serial_pass_count == 2) {
            for (const std::uint32_t key : morton_keys_) {
                ++sort_counts_[key & kSerialRadixMask];
                ++sort_counts_[kSerialRadixSize +
                               ((key >> kSerialRadixBits) & kSerialRadixMask)];
            }
        } else {
            for (const std::uint32_t key : morton_keys_) {
                ++sort_counts_[key & kSerialRadixMask];
                ++sort_counts_[kSerialRadixSize +
                               ((key >> kSerialRadixBits) & kSerialRadixMask)];
                ++sort_counts_[2U * kSerialRadixSize +
                               (key >> (2U * kSerialRadixBits))];
            }
        }

        for (unsigned pass = 0; pass < serial_pass_count; ++pass) {
            std::uint32_t* counts =
                sort_counts_.data() + pass * kSerialRadixSize;
            const std::uint32_t bucket_count =
                pass == 2 ? kSerialLastRadixSize : kSerialRadixSize;
            std::uint32_t offset = 0;
            for (std::uint32_t bucket = 0; bucket < bucket_count; ++bucket) {
                const std::uint32_t next = offset + counts[bucket];
                counts[bucket] = offset;
                offset = next;
            }

            const std::vector<Site>& input_points =
                (pass & 1U) == 0 ? points_ : sort_scratch_;
            const std::vector<std::uint32_t>& input_keys =
                (pass & 1U) == 0 ? morton_keys_ : morton_scratch_;
            std::vector<Site>& output_points =
                (pass & 1U) == 0 ? sort_scratch_ : points_;
            std::vector<std::uint32_t>& output_keys =
                (pass & 1U) == 0 ? morton_scratch_ : morton_keys_;
            const unsigned shift = pass * kSerialRadixBits;
            for (std::size_t i = 0; i < input_points.size(); ++i) {
                const std::uint32_t key = input_keys[i];
                const std::uint32_t destination =
                    counts[(key >> shift) & (bucket_count - 1U)]++;
                output_points[destination] = input_points[i];
                output_keys[destination] = key;
            }
        }
        if ((serial_pass_count & 1U) != 0) {
            points_.swap(sort_scratch_);
            morton_keys_.swap(morton_scratch_);
        }
    };

    // Keep the smaller tables for keys that already fit in two 10-bit
    // passes. Wider keys use at most three 11/11/10-bit passes.
    if (maximum_key < (1U << (2U * kRadixBits))) {
        sort_serial(std::integral_constant<unsigned, 10>{});
    } else {
        sort_serial(std::integral_constant<unsigned, 11>{});
    }
}

std::uint32_t Triangulator::Impl::morton_code(
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
