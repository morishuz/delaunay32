// SPDX-License-Identifier: MIT

#pragma once

#include "delaunay32/delaunay.hpp"
#include "edge_arena.hpp"
#include "internal.hpp"

namespace delaunay32 {

// The public facade owns this state once per instance. Kernel methods operate
// directly on it, keeping storage access and predicate inlining unchanged.
struct Triangulator::Impl {
    void set_options(TriangulationOptions options);
    void set_points(const std::vector<Point>& points);
    void set_constraints(std::vector<Constraint> constraints);
    void set_polygons(std::vector<PolygonDomain> polygons);
    TriangulationResult triangulate();

    static constexpr std::size_t kMortonLeafSize = 16;
    static constexpr std::size_t kParallelMinPoints = 50000;
    static constexpr std::size_t kParallelRadixMinPoints = 100000;
    static constexpr std::size_t kParallelMinJobPoints = 16384;
    static constexpr std::size_t kParallelJobsPerThread = 4;
    static constexpr std::size_t kMaxParallelThreads = 256;
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

    struct OuterBridges {
        std::uint32_t first_outer = 0;
        std::uint32_t last_outer = 0;
    };

    using EdgeRange = detail::EdgeArena::Range;
    using EdgeCursor = detail::EdgeArena::Cursor;

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

    // Configured boundaries and their per-run lookup/flag storage are kept
    // together. Ordinary triangulation never initializes the workspace maps.
    struct ConstraintWorkspace {
        std::vector<Constraint> segments;
        std::vector<PolygonDomain> polygons;
        std::vector<std::uint32_t> original_to_site;
        std::vector<std::uint32_t> site_edge;
        std::vector<std::uint8_t> edge_flags;

        void clear_problem() {
            segments.clear();
            polygons.clear();
        }
    };

    detail::EdgeArena arena_;
    ConstraintWorkspace constraints_;
    std::vector<Site> points_;
    std::vector<Site> sort_scratch_;
    std::vector<std::uint32_t> morton_keys_;
    std::vector<std::uint32_t> morton_scratch_;
    std::vector<std::uint64_t> wide_lifts_;
    std::vector<std::uint32_t> sort_counts_;
    std::vector<Triangle> triangles_out_;
    std::vector<std::int64_t> halfedges_out_;
    std::vector<std::uint32_t> hull_out_;
    std::vector<std::vector<Triangle>> export_scratch_;
    // Retained across calls so multi-shot clients do not rebuild OS threads.
    std::unique_ptr<detail::WorkerTeam> worker_team_;
    std::size_t worker_team_size_ = 0;
    std::size_t thread_count_ = 1;
    std::size_t active_thread_count_ = 1;
    std::size_t input_point_count_ = 0;
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
        EdgeCursor* cursor = nullptr,
        OuterBridges* bridges = nullptr);
    template <
        bool WidePredicates,
        bool ParallelAllocation = false>
    HullEdges merge_hulls(
        HullEdges left,
        HullEdges right,
        EdgeCursor* cursor = nullptr,
        OuterBridges* bridges = nullptr);
    template <
        bool WidePredicates,
        bool ParallelAllocation = false>
    DirectionalHulls merge_directional_hulls(
        const DirectionalHulls& left,
        const DirectionalHulls& right,
        bool horizontal,
        EdgeCursor* cursor);
    DirectionalHulls scan_directional_hulls(
        std::uint32_t outer_seed) const;
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
    // Keep discovery inline in both triangle-only and full-result exporters.
    bool find_export_face(
        std::uint32_t start,
        std::uint32_t& second,
        std::uint32_t& third) const {
        if ((arena_.origin[start] & kVisitedBit) != 0) {
            return false;
        }
        second = lnext(start);
        third = lnext(second);
        return start <= second && start <= third;
    }

    static std::uint32_t sym(std::uint32_t edge) { return edge ^ 1U; }
    std::uint32_t org(std::uint32_t edge) const { return arena_.origin[edge]; }
    std::uint32_t dest(std::uint32_t edge) const { return arena_.origin[sym(edge)]; }
    std::uint32_t onext(std::uint32_t edge) const { return arena_.next[edge]; }
    std::uint32_t oprev(std::uint32_t edge) const { return arena_.prev[edge]; }
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
