// SPDX-License-Identifier: MIT

#include "delaunay32/delaunay.hpp"
#include "internal.hpp"

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace delaunay32 {
using detail::ThreadBarrier;

std::size_t Triangulator::checked_flat_edge_count(
    std::size_t triangle_count) {
    if (triangle_count >
        std::numeric_limits<std::size_t>::max() / 3) {
        throw std::length_error("halfedge index range exceeded");
    }
    const std::size_t edge_count = triangle_count * 3;
    if (edge_count >
        static_cast<std::size_t>(
            std::numeric_limits<std::uint32_t>::max())) {
        throw std::length_error("internal halfedge index range exceeded");
    }
    return edge_count;
}

void Triangulator::prepare_full_export() {
    halfedges_out_.clear();
    hull_out_.clear();

    // The outer darts have no opposite output edge. edge_next_ is no longer
    // needed by face traversal, so reuse it as a dense dart-to-output-edge
    // map instead of allocating and populating a hash table.
    std::uint32_t outer = outer_seed_;
    do {
        const std::uint32_t next = lnext(outer);
        edge_origin_[outer] |= kVisitedBit;
        outer = next;
    } while (outer != outer_seed_);

    // Every retained face will populate its map entries. Outer and clipped
    // faces already carry kVisitedBit, so opposite-edge resolution can reject
    // those darts directly without clearing the complete map first.
}

void Triangulator::export_hull() {
    if (triangles_out_.empty()) {
        const auto endpoints = std::minmax_element(
            points_.begin(), points_.end(), SiteLessXY{});
        hull_out_ = {
            endpoints.first->original,
            endpoints.second->original,
        };
        return;
    }

    std::uint32_t outer = outer_seed_;
    do {
        const std::uint32_t point = org(outer) & kIndexMask;
        hull_out_.push_back(points_[point].original);
        outer = lnext(outer);
    } while (outer != outer_seed_);

    // lnext traverses the outer face clockwise. Public hull indices are
    // counterclockwise and start at the lowest original input index.
    std::reverse(hull_out_.begin(), hull_out_.end());
    const auto start = std::min_element(
        hull_out_.begin(), hull_out_.end());
    std::rotate(hull_out_.begin(), start, hull_out_.end());
}

void Triangulator::export_full_result() {
    triangles_out_.clear();
    triangles_out_.reserve(points_.size() * 2);
    prepare_full_export();

    for (const EdgeRange range : edge_ranges_) {
        for (std::uint32_t start = range.first; start < range.last; ++start) {
            std::uint32_t second = 0;
            std::uint32_t third = 0;
            if (!find_export_face(start, second, third)) {
                continue;
            }
            const std::size_t flat_edge = triangles_out_.size() * 3;
            if (flat_edge >
                static_cast<std::size_t>(
                    std::numeric_limits<std::uint32_t>::max() - 2U)) {
                throw std::length_error(
                    "internal halfedge index range exceeded");
            }
            // Keep the face's darts until every dart-to-output entry exists.
            // Reuse this buffer so adjacency can be written in output order.
            triangles_out_.push_back({start, second, third});
            edge_next_[start] = static_cast<std::uint32_t>(flat_edge);
            edge_next_[second] =
                static_cast<std::uint32_t>(flat_edge + 1);
            edge_next_[third] =
                static_cast<std::uint32_t>(flat_edge + 2);
        }
    }
    halfedges_out_.resize(
        checked_flat_edge_count(triangles_out_.size()));

    // Convert the temporary dart triples only after the complete map is ready.
    // Visit only retained faces and write halfedges in contiguous output order.
    std::size_t flat_edge = 0;
    for (Triangle& triangle : triangles_out_) {
        for (const std::uint32_t dart :
             {triangle.i0, triangle.i1, triangle.i2}) {
            const std::uint32_t opposite = sym(dart);
            halfedges_out_[flat_edge++] =
                (edge_origin_[opposite] & kVisitedBit) != 0
                    ? -1
                    : static_cast<std::int64_t>(edge_next_[opposite]);
        }
        triangle = {
            points_[org(triangle.i0)].original,
            points_[org(triangle.i1)].original,
            points_[org(triangle.i2)].original,
        };
    }
    export_hull();
}

void Triangulator::export_full_result_parallel(
    std::size_t thread_count,
    detail::WorkerTeam& workers) {
    const std::size_t worker_count =
        std::min(thread_count, edge_ranges_.size());
    if (worker_count <= 1) {
        export_full_result();
        return;
    }

    triangles_out_.clear();
    prepare_full_export();
    std::vector<std::size_t> counts(worker_count, 0);
    std::vector<std::size_t> offsets(worker_count + 1, 0);
    ThreadBarrier barrier(worker_count);

    // Count faces, allocate exact output sizes, and store face darts in the
    // triangle buffer while building the dart-to-output map. Once every worker
    // has published its map entries, convert those triples to original vertex
    // indices and write adjacent halfedges in contiguous output slices.
    // An abortable barrier makes allocation failure phase-safe.
    const auto run = [&](std::size_t worker_index) {
        try {
            const std::size_t first_range =
                edge_ranges_.size() * worker_index / worker_count;
            const std::size_t last_range =
                edge_ranges_.size() * (worker_index + 1) / worker_count;

            std::size_t count = 0;
            for (std::size_t index = first_range;
                 index < last_range;
                 ++index) {
                const EdgeRange range = edge_ranges_[index];
                for (std::uint32_t start = range.first;
                     start < range.last;
                     ++start) {
                    std::uint32_t second = 0;
                    std::uint32_t third = 0;
                    if (find_export_face(start, second, third)) {
                        ++count;
                    }
                }
            }
            counts[worker_index] = count;

            if (!barrier.wait()) {
                return;
            }
            if (worker_index == 0) {
                for (std::size_t i = 0; i < worker_count; ++i) {
                    offsets[i + 1] = offsets[i] + counts[i];
                }
                triangles_out_.resize(offsets.back());
                halfedges_out_.resize(
                    checked_flat_edge_count(offsets.back()));
            }

            if (!barrier.wait()) {
                return;
            }
            std::size_t triangle_index = offsets[worker_index];
            for (std::size_t index = first_range;
                 index < last_range;
                 ++index) {
                const EdgeRange range = edge_ranges_[index];
                for (std::uint32_t start = range.first;
                     start < range.last;
                     ++start) {
                    std::uint32_t second = 0;
                    std::uint32_t third = 0;
                    if (!find_export_face(start, second, third)) {
                        continue;
                    }
                    const std::size_t flat_edge = triangle_index * 3;
                    triangles_out_[triangle_index] = {start, second, third};
                    edge_next_[start] =
                        static_cast<std::uint32_t>(flat_edge);
                    edge_next_[second] =
                        static_cast<std::uint32_t>(flat_edge + 1);
                    edge_next_[third] =
                        static_cast<std::uint32_t>(flat_edge + 2);
                    ++triangle_index;
                }
            }

            if (!barrier.wait()) {
                return;
            }
            for (std::size_t i = offsets[worker_index];
                 i < offsets[worker_index + 1];
                 ++i) {
                Triangle& triangle = triangles_out_[i];
                std::size_t flat_edge = i * 3;
                for (const std::uint32_t dart :
                     {triangle.i0, triangle.i1, triangle.i2}) {
                    const std::uint32_t opposite = sym(dart);
                    halfedges_out_[flat_edge++] =
                        (edge_origin_[opposite] & kVisitedBit) != 0
                            ? -1
                            : static_cast<std::int64_t>(edge_next_[opposite]);
                }
                triangle = {
                    points_[org(triangle.i0)].original,
                    points_[org(triangle.i1)].original,
                    points_[org(triangle.i2)].original,
                };
            }
        } catch (...) {
            barrier.abort();
            throw;
        }
    };
    workers.run(worker_count, run);
    export_hull();
}

void Triangulator::finish_full_export() {
    if (active_thread_count_ > 1) {
        export_full_result_parallel(
            active_thread_count_, *worker_team_);
    } else {
        export_full_result();
    }
}

TriangulationResult Triangulator::make_result(
    PredicateWidth predicate_width,
    std::vector<std::uint32_t>&& representatives) {
    TriangulationResult result;
    result.triangles = std::move(triangles_out_);
    result.halfedges = std::move(halfedges_out_);
    result.hull = std::move(hull_out_);
    result.representatives = std::move(representatives);
    result.report.predicate_width = predicate_width;
    result.report.actual_thread_count = active_thread_count_;
    result.report.input_points = input_point_count_;
    result.report.unique_points = points_.size();
    result.report.collapsed_points = input_point_count_ - points_.size();
    return result;
}

}  // namespace delaunay32
