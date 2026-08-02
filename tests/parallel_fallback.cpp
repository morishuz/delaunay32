// SPDX-License-Identifier: MIT

#include "support.hpp"

#include "delaunay32/delaunay.hpp"

#include <algorithm>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void require_fallback_matches(
    const std::vector<delaunay32::Point>& points,
    const std::string& label) {
    using delaunay32::Triangle;
    namespace support = delaunay32::benchmark_support;

    delaunay32::Triangulator serial(1);
    const std::vector<Triangle> expected = serial.triangulate_int(points);

    // This executable links a private library variant with an arena large
    // enough to build substantial partial topology but too small to finish
    // this input. A successful result therefore proves that the partial
    // topology was discarded and the input was rebuilt through the growable
    // serial allocator.
    delaunay32::Triangulator parallel_with_fallback(2);
    const std::vector<Triangle> candidate =
        parallel_with_fallback.triangulate_int(points);
    require(
        support::meshes_equal(expected, candidate),
        label + ": serial fallback produced a different triangle set");

    std::string validation_error;
    require(
        support::validate_mesh(points, candidate, validation_error),
        label + ": serial fallback produced an invalid mesh: " +
            validation_error);

    const std::vector<Triangle> repeated =
        parallel_with_fallback.triangulate_int(points);
    require(
        support::meshes_equal(expected, repeated),
        label + ": reused triangulator produced a different fallback mesh");

    const delaunay32::TriangulationResult full =
        parallel_with_fallback.triangulate_int_full(points);
    require(
        support::meshes_equal(expected, full.triangles),
        label + ": full result changed the fallback mesh");
    require(
        full.halfedges.size() == full.triangles.size() * 3 &&
            !full.hull.empty() &&
            full.representatives.size() == points.size(),
        label + ": full fallback result omitted auxiliary data");
    require(
        full.actual_thread_count == 1,
        label + ": full result did not report serial fallback execution");
}

}  // namespace

int main() {
    try {
        using delaunay32::Point;
        namespace support = delaunay32::benchmark_support;

        const std::vector<Point> points = support::generate_points(
            support::Dataset::Uniform,
            60000,
            0xfa11bacULL,
            20000);
        require_fallback_matches(points, "int64 predicates");

        std::vector<Point> grid_points;
        grid_points.reserve(60000);
        for (std::int32_t y = 0; y < 200; ++y) {
            for (std::int32_t x = 0; x < 300; ++x) {
                grid_points.push_back({x, y});
            }
        }
        require_fallback_matches(grid_points, "cocircular grid");

#if defined(__SIZEOF_INT128__)
        const std::vector<Point> wide_points = support::generate_points(
            support::Dataset::Diagonal,
            60000,
            0xfa11badULL,
            100000);
        const auto [minimum_x, maximum_x] = std::minmax_element(
            wide_points.begin(),
            wide_points.end(),
            [](const Point& a, const Point& b) { return a.x < b.x; });
        const auto [minimum_y, maximum_y] = std::minmax_element(
            wide_points.begin(),
            wide_points.end(),
            [](const Point& a, const Point& b) { return a.y < b.y; });
        const std::uint64_t x_span = static_cast<std::uint64_t>(
            static_cast<std::int64_t>(maximum_x->x) - minimum_x->x);
        const std::uint64_t y_span = static_cast<std::uint64_t>(
            static_cast<std::int64_t>(maximum_y->y) - minimum_y->y);
        require(
            delaunay32::Triangulator::predicate_width_for_spans(
                x_span, y_span) == delaunay32::PredicateWidth::Int128,
            "wide fallback input did not select int128 predicates");
        require_fallback_matches(wide_points, "int128 predicates");
#endif

        std::cout << "parallel edge-arena fallback passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "parallel edge-arena fallback failed: "
                  << error.what() << '\n';
        return 1;
    }
}
