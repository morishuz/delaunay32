// SPDX-License-Identifier: MIT

#include "support.hpp"
#include "delaunay32/extras/sampling.hpp"
#include "delaunay32/extras/svg.hpp"

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using namespace delaunay32;
using namespace delaunay32::extras;

// Input generation is outside the timer. Each operation includes destruction
// of its output, and triangulation reuses its instance and includes set_points.
template <typename Operation>
void measure(const char* name, std::size_t size, std::size_t repetitions,
             Operation operation) {
    std::uint64_t checksum = operation(); // Warm up and surface invalid inputs.
    std::vector<double> samples;
    for (std::size_t i = 0; i < repetitions; ++i) {
        const auto start = std::chrono::steady_clock::now();
        checksum += operation();
        const auto end = std::chrono::steady_clock::now();
        samples.push_back(
            std::chrono::duration<double, std::milli>(end - start).count());
    }
    std::sort(samples.begin(), samples.end());
    std::cout << name << ',' << size << ',' << std::fixed << std::setprecision(4)
              << samples[samples.size() / 2] << ',' << repetitions << ','
              << checksum << '\n';
}

void run(bool quick) {
    const std::size_t point_count = quick ? 10000 : 100000;
    const std::size_t repetitions = quick ? 3 : 7;
    const auto points = benchmark_support::generate_points(
        benchmark_support::Dataset::Uniform, point_count, 1741, 20000);
    std::vector<Constraint> constraints;
    for (std::uint32_t i = 1; i <= 128; ++i) {
        constraints.push_back({0, i}); // Noncrossing rays with shared origin.
    }
    std::cout << "workload,input_size,median_ms,samples,checksum\n";
    Triangulator triangulator;
    for (const auto detail : {ResultDetail::Triangles, ResultDetail::Full}) {
        triangulator.set_options({1, detail});
        measure(detail == ResultDetail::Full ? "triangulate_full" : "triangulate",
                point_count, repetitions, [&] {
            triangulator.set_points(points);
            const auto result = triangulator.triangulate();
            return result.triangles.size() + result.halfedges.size();
        });
        measure(detail == ResultDetail::Full ? "constraints_full" : "constraints",
                point_count, repetitions, [&] {
            triangulator.set_points(points);
            triangulator.set_constraints(constraints);
            const auto result = triangulator.triangulate();
            return result.triangles.size() + result.halfedges.size();
        });
    }
    triangulator.set_points(points);
    const auto mesh = triangulator.triangulate();
    measure("svg_record", mesh.triangles.size(), repetitions, [&] {
        Svg svg(800, 600);
        svg.draw_triangles(points, mesh.triangles);
        return mesh.triangles.size();
    });
    Svg svg(800, 600);
    svg.draw_triangles(points, mesh.triangles);
    measure("svg_serialize", mesh.triangles.size(), repetitions, [&] {
        return svg.to_svg().size();
    });

    // Two equal rectangles with a 100x-wide empty gap. Rejection from the
    // combined bounds still succeeds with the default per-point attempt cap.
    PointSampler sampler;
    sampler.set_polygon_interiors(
        std::vector<Point>{{0, 0}, {100, 0}, {100, 100}, {0, 100},
                           {10000, 0}, {10100, 0}, {10100, 100}, {10000, 100}},
        {{{0, 1, 2, 3}, {}}, {{4, 5, 6, 7}, {}}});
    UniformSamplingOptions uniform;
    uniform.point_count = quick ? 1000 : 10000;
    measure("uniform_disconnected", uniform.point_count, repetitions, [&] {
        return sampler.generate_uniform(uniform).size();
    });
    BlueNoiseSamplingOptions blue;
    blue.point_count = quick ? 100 : 1000;
    measure("blue_noise_disconnected", blue.point_count, repetitions, [&] {
        return sampler.generate_blue_noise(blue).size();
    });
    sampler.set_polygon_interiors(
        std::vector<Point>{{0, 0}, {100, 0}, {100, 100}, {0, 100},
                           {110, 0}, {210, 0}, {210, 100}, {110, 100}},
        {{{0, 1, 2, 3}, {}}, {{4, 5, 6, 7}, {}}});
    measure("uniform_compact_domains", uniform.point_count, repetitions, [&] {
        return sampler.generate_uniform(uniform).size();
    });
    measure("blue_noise_compact_domains", blue.point_count, repetitions, [&] {
        return sampler.generate_blue_noise(blue).size();
    });
    sampler.set_bounds({0, 100, 0, 100});
    measure("uniform_bounds", uniform.point_count, repetitions, [&] {
        return sampler.generate_uniform(uniform).size();
    });
    measure("blue_noise_bounds", blue.point_count, repetitions, [&] {
        return sampler.generate_blue_noise(blue).size();
    });
}
} // namespace

int main(int argc, char** argv) {
    try {
        const bool quick = argc == 2 && std::string(argv[1]) == "--quick";
        if (argc != 1 && !quick) {
            throw std::invalid_argument("usage: delaunay_workload_benchmark [--quick]");
        }
        run(quick);
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
