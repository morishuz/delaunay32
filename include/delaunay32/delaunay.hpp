// SPDX-License-Identifier: MIT

#pragma once

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

namespace delaunay32 {

namespace detail {
class WorkerTeam;
}

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

    // Special members are out-of-line: WorkerTeam is incomplete in this header.
    Triangulator();
    ~Triangulator();
    Triangulator(const Triangulator&) = delete;
    Triangulator& operator=(const Triangulator&) = delete;
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
    static constexpr std::size_t kMortonLeafSize = 16;
    static constexpr std::size_t kParallelMinPoints = 50000;
    static constexpr std::size_t kParallelRadixMinPoints = 100000;
    static constexpr std::size_t kParallelMinJobPoints = 16384;
    static constexpr std::size_t kParallelJobsPerThread = 4;
    static constexpr std::size_t kMaxParallelThreads = 256;
    static constexpr std::uint32_t kEdgeBlockDarts = 4096;
    // Deleted darts remain in the append-only arena. Production benchmarks
    // use about 8.1 darts/point, so 9 provides measured headroom but is not a
    // geometric upper bound. Serial allocation can grow; parallel allocation
    // is fixed because workers hold indices into these arrays. Exhaustion
    // transparently discards the partial topology and retries serially.
    static constexpr std::size_t kEdgeArenaDartsPerPoint = 9;
    static constexpr std::uint32_t kVisitedBit = 0x80000000U;
    static constexpr std::uint32_t kIndexMask = 0x7fffffffU;
    static constexpr std::uint32_t kDeletedEdge = UINT32_MAX;
    static constexpr std::uint8_t kConstrainedBit = 0x01U;
    static constexpr std::uint8_t kLegalizationQueuedBit = 0x02U;

    struct Site {
        std::int32_t x = 0;
        std::int32_t y = 0;
        std::uint32_t original = 0;
        std::uint32_t lift = 0;
    };

    struct SiteLessXY {
        bool operator()(const Site& a, const Site& b) const noexcept {
            if (a.x != b.x) {
                return a.x < b.x;
            }
            if (a.y != b.y) {
                return a.y < b.y;
            }
            return a.original < b.original;
        }
    };

    struct HullEdges {
        std::uint32_t left_outer = 0;
        std::uint32_t right_outer = 0;
    };

    struct DirectionalHulls {
        HullEdges x;
        HullEdges y;
    };

    struct EdgeRange {
        std::uint32_t first = 0;
        std::uint32_t last = 0;
    };

    struct EdgeCursor {
        std::uint32_t next = 0;
        std::uint32_t end = 0;
        std::uint32_t range_first = 0;
        std::atomic<std::size_t>* block_counter = nullptr;
        std::vector<EdgeRange> ranges;
    };

    struct MortonSplit {
        std::size_t middle = 0;
        unsigned split_bit = 0;
        bool valid = false;
    };

    struct ParallelNode {
        std::size_t first = 0;
        std::size_t last = 0;
        std::size_t left = 0;
        std::size_t right = 0;
        unsigned split_bit = 0;
        bool leaf = true;
        DirectionalHulls hull;
    };

    std::vector<Site> points_;
    std::vector<Site> sort_scratch_;
    std::vector<std::uint32_t> morton_keys_;
    std::vector<std::uint32_t> morton_scratch_;
    std::vector<std::uint64_t> wide_lifts_;
    std::vector<std::uint32_t> sort_counts_;
    std::vector<std::uint32_t> edge_origin_;
    std::vector<std::uint32_t> edge_next_;
    std::vector<std::uint32_t> edge_prev_;
    std::vector<std::uint8_t> edge_constrained_;
    std::vector<EdgeRange> edge_ranges_;
    std::vector<std::uint32_t> original_to_site_;
    std::vector<std::uint32_t> site_edge_;
    std::vector<Triangle> triangles_out_;
    std::vector<std::int64_t> halfedges_out_;
    std::vector<std::uint32_t> hull_out_;
    std::vector<Constraint> constraints_;
    std::vector<PolygonDomain> polygons_;
    std::vector<std::vector<Triangle>> export_scratch_;
    // Retained across calls so multi-shot clients do not rebuild OS threads.
    std::unique_ptr<detail::WorkerTeam> worker_team_;
    std::size_t worker_team_size_ = 0;
    std::size_t edge_count_ = 0;
    std::size_t thread_count_ = 1;
    std::size_t active_thread_count_ = 1;
    std::size_t input_point_count_ = 0;
    std::size_t edge_capacity_limit_ = 0;
    std::uint32_t outer_seed_ = 0;
    std::int32_t min_x_ = 0;
    std::int32_t min_y_ = 0;
    std::int32_t max_x_ = 0;
    std::int32_t max_y_ = 0;
    bool int64_wide_intermediates_ = false;
    bool problem_ready_ = false;
    ResultDetail result_detail_ = ResultDetail::Triangles;

    static void require_point_count(std::size_t point_count);
    void require_ready_problem() const;
    detail::WorkerTeam* ensure_worker_team(std::size_t thread_count);
    void load_int_points(const std::vector<Point>& points);
    PredicateWidth build_loaded_topology(
        std::vector<std::uint32_t>* representatives = nullptr);
    TriangulationResult make_result(
        PredicateWidth predicate_width,
        std::vector<std::uint32_t>&& representatives);
    void sort_points_morton(
        std::size_t thread_count,
        std::uint32_t maximum_key,
        detail::WorkerTeam* workers);
    void mark_outer_face();
    template <
        bool WidePredicates,
        bool ParallelAllocation = false>
    HullEdges build_range(
        std::size_t first,
        std::size_t last,
        EdgeCursor* cursor);
    template <
        bool WidePredicates,
        bool ParallelAllocation = false>
    DirectionalHulls build_morton_range(
        std::size_t first,
        std::size_t last,
        EdgeCursor* cursor = nullptr);
    template <
        bool WidePredicates,
        bool ParallelAllocation = false>
    inline HullEdges merge_hulls_inline(
        HullEdges left,
        HullEdges right,
        EdgeCursor* cursor = nullptr);
    template <
        bool WidePredicates,
        bool ParallelAllocation = false>
    HullEdges merge_hulls(
        HullEdges left,
        HullEdges right,
        EdgeCursor* cursor = nullptr);
    DirectionalHulls scan_directional_hulls(
        std::uint32_t outer_seed) const;
    template <bool Horizontal>
    DirectionalHulls scan_merged_hulls(
        std::uint32_t outer_seed,
        HullEdges split_hull) const;
    static std::uint32_t morton_code(std::uint32_t x, std::uint32_t y);
    MortonSplit find_morton_split(std::size_t first, std::size_t last) const;
    std::size_t add_parallel_node(
        std::size_t first,
        std::size_t last,
        std::size_t target_size,
        std::vector<ParallelNode>& nodes,
        std::vector<std::size_t>& leaves) const;
    template <bool WidePredicates>
    DirectionalHulls build_parallel(
        std::size_t thread_count,
        detail::WorkerTeam& workers);
    void acquire_edge_block(EdgeCursor& cursor);
    static void finish_edge_cursor(EdgeCursor& cursor);

    template <bool ParallelAllocation = false>
    std::uint32_t make_edge(
        std::uint32_t origin,
        std::uint32_t destination,
        EdgeCursor* cursor = nullptr);
    template <bool ParallelAllocation = false>
    std::uint32_t connect(
        std::uint32_t a,
        std::uint32_t b,
        EdgeCursor* cursor = nullptr);
    void splice(std::uint32_t a, std::uint32_t b);
    void delete_edge(std::uint32_t edge);
    void flip_edge(std::uint32_t edge);
    bool is_live_edge(std::uint32_t edge) const;
    bool is_constrained(std::uint32_t edge) const;
    void mark_constrained(std::uint32_t edge);
    bool left_face_opposite(
        std::uint32_t edge,
        std::uint32_t& opposite) const;
    bool can_flip(std::uint32_t edge) const;
    bool active_in_circle(
        std::uint32_t a,
        std::uint32_t b,
        std::uint32_t c,
        std::uint32_t d) const;
    void build_constraint_indices(
        const std::vector<std::uint32_t>& representatives,
        std::size_t original_point_count);
    void recover_constraints(
        const std::vector<Constraint>& constraints);
    void recover_constraint(
        std::uint32_t a,
        std::uint32_t b,
        std::vector<std::uint32_t>& legalization_queue);
    void queue_constraint_legalization(
        std::uint32_t edge,
        std::vector<std::uint32_t>& legalization_queue);
    void seed_constraint_legalization(
        std::uint32_t edge,
        std::vector<std::uint32_t>& legalization_queue);
    std::uint32_t find_edge(
        std::uint32_t origin,
        std::uint32_t destination) const;
    std::uint32_t first_collinear_edge(
        std::uint32_t origin,
        std::uint32_t destination) const;
    std::vector<std::uint32_t> crossed_edges(
        std::uint32_t a,
        std::uint32_t b,
        std::uint32_t& reached) const;
    bool properly_intersects(
        std::uint32_t edge,
        std::uint32_t a,
        std::uint32_t b) const;
    void legalize_unconstrained_edges(
        std::vector<std::uint32_t>& legalization_queue);
    std::vector<std::vector<std::uint32_t>> prepare_polygon_rings(
        const std::vector<std::uint32_t>& outer_ring,
        const std::vector<std::vector<std::uint32_t>>& holes) const;
    std::vector<std::vector<std::vector<std::uint32_t>>>
    prepare_polygon_domains() const;
    std::uint32_t first_boundary_edge(
        std::uint32_t origin,
        std::uint32_t destination) const;
    void mark_polygon_excluded_faces(
        const std::vector<
            std::vector<std::vector<std::uint32_t>>>& domains);
    void finish_triangle_export();
    void finish_full_export();
    void export_triangles();
    void export_triangles_parallel(
        std::size_t thread_count,
        detail::WorkerTeam& workers);
    void export_full_result();
    void export_full_result_parallel(
        std::size_t thread_count,
        detail::WorkerTeam& workers);
    void prepare_full_export();
    void export_hull();
    static std::size_t checked_flat_edge_count(
        std::size_t triangle_count);
    bool find_export_face(
        std::uint32_t start,
        std::uint32_t& second,
        std::uint32_t& third) const;

    static std::uint32_t sym(std::uint32_t edge) { return edge ^ 1U; }
    std::uint32_t org(std::uint32_t edge) const { return edge_origin_[edge]; }
    std::uint32_t dest(std::uint32_t edge) const { return edge_origin_[sym(edge)]; }
    std::uint32_t onext(std::uint32_t edge) const { return edge_next_[edge]; }
    std::uint32_t oprev(std::uint32_t edge) const { return edge_prev_[edge]; }
    std::uint32_t lnext(std::uint32_t edge) const { return oprev(sym(edge)); }

    std::int64_t orient(std::uint32_t a, std::uint32_t b, std::uint32_t c) const;
    bool left_of(
        std::uint32_t point,
        std::uint32_t edge) const;
    bool right_of(
        std::uint32_t point,
        std::uint32_t edge) const;
    template <bool WidePredicates>
    bool in_circle(
        std::uint32_t a,
        std::uint32_t b,
        std::uint32_t c,
        std::uint32_t d) const;
};

}  // namespace delaunay32
