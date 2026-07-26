// SPDX-License-Identifier: MIT

#include "delaunay32/delaunay.hpp"
#include "internal.hpp"

#include <algorithm>
#include <vector>

namespace delaunay32 {
using detail::elapsed_ms;
using detail::ProfileClock;
using detail::ThreadBarrier;

// Face discovery and materialization are isolated from the topology kernel.
void Triangulator::mark_outer_face() {
    std::uint32_t outer = outer_seed_;
    do {
        edge_origin_[outer] |= kVisitedBit;
        outer = lnext(outer);
    } while (outer != outer_seed_);
}

template <bool CountOperations>
void Triangulator::export_triangles(Profile* profile) {
    ProfileClock::time_point phase_start;
    if (profile != nullptr) {
        phase_start = ProfileClock::now();
    }
    triangles_out_.clear();
    triangles_out_.reserve(points_.size() * 2);
    if constexpr (CountOperations) {
        profile->export_darts_scanned += edge_count_;
    }

    mark_outer_face();
    if (profile != nullptr) {
        profile->export_prepare_ms = elapsed_ms(phase_start);
        phase_start = ProfileClock::now();
    }

    for (const EdgeRange range : edge_ranges_) {
        for (std::uint32_t start = range.first; start < range.last; ++start) {
            if ((edge_origin_[start] & kVisitedBit) != 0) {
                continue;
            }
            const std::uint32_t second = lnext(start);
            const std::uint32_t third = lnext(second);
            const std::uint32_t a = org(start);
            const std::uint32_t b = org(second);
            const std::uint32_t c = org(third);
            edge_origin_[start] |= kVisitedBit;
            edge_origin_[second] |= kVisitedBit;
            edge_origin_[third] |= kVisitedBit;
            triangles_out_.push_back({
                points_[a].original,
                points_[b].original,
                points_[c].original,
            });
        }
    }
    if (profile != nullptr) {
        profile->export_write_ms = elapsed_ms(phase_start);
    }
}

bool Triangulator::find_export_face(
    std::uint32_t start,
    std::uint32_t& second,
    std::uint32_t& third) const {
    if ((edge_origin_[start] & kVisitedBit) != 0) {
        return false;
    }
    second = lnext(start);
    third = lnext(second);
    return start <= second && start <= third;
}

void Triangulator::export_triangles_parallel(
    std::size_t thread_count,
    Profile* profile,
    detail::WorkerTeam& workers) {
    ProfileClock::time_point phase_start;
    if (profile != nullptr) {
        phase_start = ProfileClock::now();
    }
    const std::size_t worker_count =
        std::min(thread_count, edge_ranges_.size());
    if (worker_count <= 1) {
        export_triangles<false>(profile);
        return;
    }

    mark_outer_face();

    export_scratch_.resize(worker_count);
    const std::size_t expected_triangles =
        (points_.size() * 2 + worker_count - 1) / worker_count;
    for (std::vector<Triangle>& buffer : export_scratch_) {
        buffer.clear();
        buffer.reserve(expected_triangles);
    }
    std::vector<std::size_t> counts(worker_count, 0);
    std::vector<std::size_t> offsets(worker_count + 1, 0);
    ThreadBarrier barrier(worker_count);

    // Export has three synchronized phases in one worker lifetime:
    //
    //   one face walk into worker buffers -> prefix sum/resize
    //                                     -> parallel contiguous copies
    //
    // Contiguous range partitions preserve the previous deterministic face
    // order. A face is owned by its smallest dart index, so workers never
    // mutate visited state while scanning.
    ProfileClock::time_point scan_start;
    ProfileClock::time_point scan_end;
    ProfileClock::time_point write_start;
    if (profile != nullptr) {
        profile->export_prepare_ms = elapsed_ms(phase_start);
        scan_start = ProfileClock::now();
    }
    const auto run = [&](std::size_t worker_index) {
        std::vector<Triangle> buffer;
        buffer.swap(export_scratch_[worker_index]);
        const std::size_t first_range =
            edge_ranges_.size() * worker_index / worker_count;
        const std::size_t last_range =
            edge_ranges_.size() * (worker_index + 1) / worker_count;
        for (std::size_t index = first_range;
             index < last_range;
             ++index) {
            const EdgeRange range = edge_ranges_[index];
            std::uint32_t second = 0;
            std::uint32_t third = 0;
            for (std::uint32_t start = range.first;
                 start < range.last;
                 ++start) {
                if (find_export_face(start, second, third)) {
                    buffer.push_back({
                        points_[org(start)].original,
                        points_[org(second)].original,
                        points_[org(third)].original,
                    });
                }
            }
        }
        counts[worker_index] = buffer.size();

        barrier.wait();
        if (worker_index == 0) {
            if (profile != nullptr) {
                scan_end = ProfileClock::now();
            }
            for (std::size_t i = 0; i < worker_count; ++i) {
                offsets[i + 1] = offsets[i] + counts[i];
            }
            triangles_out_.resize(offsets.back());
            if (profile != nullptr) {
                write_start = ProfileClock::now();
            }
        }

        barrier.wait();
        std::copy(
            buffer.begin(),
            buffer.end(),
            triangles_out_.begin() +
                static_cast<std::ptrdiff_t>(offsets[worker_index]));
        buffer.clear();
        export_scratch_[worker_index].swap(buffer);
    };

    workers.run(worker_count, run);
    if (profile != nullptr) {
        const ProfileClock::time_point write_end =
            ProfileClock::now();
        profile->export_count_ms =
            std::chrono::duration<double, std::milli>(
                scan_end - scan_start)
                .count();
        profile->export_prepare_ms +=
            std::chrono::duration<double, std::milli>(
                write_start - scan_end)
                .count();
        profile->export_write_ms =
            std::chrono::duration<double, std::milli>(
                write_end - write_start)
                .count();
    }
}

// The API translation unit references both serial profiling specializations.
template void Triangulator::export_triangles<false>(
    Profile*);
template void Triangulator::export_triangles<true>(
    Profile*);

}  // namespace delaunay32
