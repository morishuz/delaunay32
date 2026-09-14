// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

namespace delaunay32 {

// A signed integer site. Absolute offsets do not affect predicate selection.
struct Point {
    std::int32_t x = 0;
    std::int32_t y = 0;
};

// Indices into the caller's original point array, counterclockwise in the
// integer coordinates used for triangulation.
struct Triangle {
    std::uint32_t i0 = 0;
    std::uint32_t i1 = 0;
    std::uint32_t i2 = 0;
};

// An undirected edge that must be present in a constrained triangulation.
// Endpoints index the caller's original point array.
struct Constraint {
    std::uint32_t i0 = 0;
    std::uint32_t i1 = 0;
};

// Rings contain indices into the configured point array. Their closing edge
// is implicit; a repeated first index at the end is optional.
struct PolygonDomain {
    std::vector<std::uint32_t> outer_ring;
    std::vector<std::vector<std::uint32_t>> holes;
};

enum class PredicateWidth {
    Int64,
    Int128,
    Unsupported,
};

enum class ResultDetail {
    Triangles,
    Full,
};

struct TriangulationOptions {
    // One preserves the serial path. Zero selects the hardware thread count
    // and still falls back to serial for small inputs.
    std::size_t thread_count = 1;
    ResultDetail result_detail = ResultDetail::Triangles;
};

struct TriangulationReport {
    PredicateWidth predicate_width = PredicateWidth::Unsupported;
    std::size_t actual_thread_count = 1;
    std::size_t input_points = 0;
    std::size_t unique_points = 0;
    std::size_t collapsed_points = 0;
};

struct TriangulationResult {
    // For each flattened edge e = 3 * triangle + local_edge,
    // halfedges[e] is the oppositely directed neighboring edge, or -1 on the
    // convex hull or a clipped domain boundary. Local edges are i0->i1,
    // i1->i2, and i2->i0.
    std::vector<Triangle> triangles;
    std::vector<std::int64_t> halfedges;

    // Original input indices around the convex hull, counterclockwise on the
    // triangulation grid and rotated to start at the lowest input index. A
    // collinear result contains its two geometric endpoints.
    std::vector<std::uint32_t> hull;

    // One entry per input point. Coincident inputs map to the lowest original
    // index retained at that coordinate.
    std::vector<std::uint32_t> representatives;

    TriangulationReport report;
};

// Integer divide-and-conquer triangulator using a compact two-dart primal
// edge ring. Each configured problem is consumed by triangulate(), including
// failed runs. Call set_points() to begin another problem. Allocations and
// worker threads are retained across problems. An instance is not concurrent.
class Triangulator {
public:
    // Largest equal x/y spans certified by the conservative runtime bounds.
    // Asymmetric inputs may support a larger span on one axis.
    static constexpr std::int64_t kFastCoordinateSpan = 29609;
#if defined(__SIZEOF_INT128__)
    static constexpr std::int64_t kMaxCoordinateSpan = 1940470527;
#else
    static constexpr std::int64_t kMaxCoordinateSpan = kFastCoordinateSpan;
#endif

    static PredicateWidth predicate_width_for_spans(
        std::uint64_t x_span,
        std::uint64_t y_span) noexcept;
    static bool int64_wide_intermediates_for_spans(
        std::uint64_t x_span,
        std::uint64_t y_span) noexcept;

    // Special members are out-of-line so implementation details stay private.
    Triangulator();
    ~Triangulator();
    Triangulator(const Triangulator&) = delete;
    Triangulator& operator=(const Triangulator&) = delete;
    // Moving transfers the problem and retained storage. A moved-from instance
    // can be configured again, starting with default options.
    Triangulator(Triangulator&&) noexcept;
    Triangulator& operator=(Triangulator&&) noexcept;

    void set_options(TriangulationOptions options);

    // Starts a new problem and clears its previous geometry and run state.
    void set_points(const std::vector<Point>& points);
    void set_constraints(std::vector<Constraint> constraints);
    void set_polygons(std::vector<PolygonDomain> polygons);

    // Coincident sites are collapsed deterministically. Triangles reference
    // the lowest original input index for each retained site.
    TriangulationResult triangulate();

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
    Impl& implementation();
};

}  // namespace delaunay32
