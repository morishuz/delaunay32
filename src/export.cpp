// SPDX-License-Identifier: MIT

#include "delaunay32/delaunay.hpp"
#include "internal.hpp"

#include <algorithm>
#include <vector>

namespace delaunay32 {
using detail::ThreadBarrier;

// Face discovery and materialization are isolated from the topology kernel.
void Triangulator::mark_outer_face() {
    std::uint32_t outer = outer_seed_;
    do {
        edge_origin_[outer] |= kVisitedBit;
        outer = lnext(outer);
    } while (outer != outer_seed_);
}

void Triangulator::export_triangles() {
    triangles_out_.clear();
    triangles_out_.reserve(points_.size() * 2);

    mark_outer_face();

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
    detail::WorkerTeam& workers) {
    const std::size_t worker_count =
        std::min(thread_count, edge_ranges_.size());
    if (worker_count <= 1) {
        export_triangles();
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

    // Export keeps all phases in one worker lifetime. The barrier is abortable
    // so a failed buffer growth or output resize releases waiting peers.
    //
    //   face walks into worker buffers -> caller prefix sum/resize
    //                                  -> parallel contiguous copies
    //
    // Contiguous range partitions preserve the previous deterministic face
    // order. A face is owned by its smallest dart index, so workers never
    // mutate visited state while scanning.
    const auto run = [&](std::size_t worker_index) {
        try {
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

            if (!barrier.wait()) {
                return;
            }
            if (worker_index == 0) {
                for (std::size_t i = 0; i < worker_count; ++i) {
                    offsets[i + 1] = offsets[i] + counts[i];
                }
                triangles_out_.resize(offsets.back());
            }

            if (!barrier.wait()) {
                return;
            }
            std::copy(
                buffer.begin(),
                buffer.end(),
                triangles_out_.begin() +
                    static_cast<std::ptrdiff_t>(offsets[worker_index]));
            buffer.clear();
            export_scratch_[worker_index].swap(buffer);
        } catch (...) {
            barrier.abort();
            throw;
        }
    };
    workers.run(worker_count, run);
}

}  // namespace delaunay32
