// SPDX-License-Identifier: MIT

#include "support.hpp"

#include "delaunay32/delaunay.hpp"

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace delaunay32 {
namespace {

using benchmark_support::Dataset;

#if defined(__SIZEOF_INT128__)
constexpr std::int32_t kBenchmarkDomain = 100000;
#else
constexpr std::int32_t kBenchmarkDomain = 20000;
#endif
constexpr std::size_t kParallelMinPoints = 50000;
constexpr std::size_t kMaxParallelThreads = 256;

void configure_triangulator(
    Triangulator& triangulator,
    std::size_t thread_count) {
    TriangulationOptions options;
    options.thread_count = thread_count;
    triangulator.set_options(options);
}

std::vector<Triangle> triangulate_points(
    Triangulator& triangulator,
    const std::vector<Point>& points) {
    triangulator.set_points(points);
    return triangulator.triangulate().triangles;
}

enum class OutputFormat {
    Text,
    Markdown,
    Csv,
};

struct Options {
    bool quick = false;
    bool reuse = false;
    OutputFormat format = OutputFormat::Text;
    std::vector<std::size_t> sizes = {1000, 10000, 100000, 1000000};
};

struct Row {
    Dataset dataset = Dataset::Uniform;
    std::size_t points = 0;
    std::size_t triangles = 0;
    std::size_t auto_threads = 1;
    std::size_t samples = 0;
    bool reuse = true;
    double serial_ms = 0.0;
    double parallel_ms = 0.0;
};

std::size_t automatic_thread_count(std::size_t point_count) {
    // Mirrors Triangulator's initial automatic scheduling policy. Benchmark
    // generators produce unique points, so compaction cannot change the count.
    if (point_count < kParallelMinPoints) {
        return 1;
    }
    return std::min<std::size_t>(
        kMaxParallelThreads,
        std::max<std::size_t>(1, std::thread::hardware_concurrency()));
}

void print_usage(const char* program) {
    std::cout
        << "Usage: " << program << " [options]\n"
        << "  --quick              Benchmark through 100,000 points\n"
        << "  --fresh              New Delaunay32 instance (default)\n"
        << "  --reuse              Retain working storage between samples\n"
        << "  --format FORMAT      text, markdown, or csv\n"
        << "  --markdown           Alias for --format markdown\n"
        << "  --sizes N[,N...]     Override point counts\n"
        << "  --help                Show this message\n";
}

std::vector<std::size_t> parse_sizes(const std::string& value) {
    std::vector<std::size_t> sizes;
    std::size_t first = 0;
    while (first < value.size()) {
        const std::size_t comma = value.find(',', first);
        const std::string token =
            value.substr(first, comma == std::string::npos
                                    ? std::string::npos
                                    : comma - first);
        if (token.empty()) {
            throw std::invalid_argument("empty point count in --sizes");
        }
        std::size_t consumed = 0;
        const unsigned long long parsed = std::stoull(token, &consumed);
        if (consumed != token.size()) {
            throw std::invalid_argument(
                "invalid point count in --sizes: " + token);
        }
        if (parsed < 3 ||
            parsed > std::numeric_limits<std::uint32_t>::max()) {
            throw std::invalid_argument(
                "benchmark point counts must be between 3 and uint32 max");
        }
        sizes.push_back(static_cast<std::size_t>(parsed));
        if (comma == std::string::npos) {
            break;
        }
        first = comma + 1;
    }
    if (sizes.empty()) {
        throw std::invalid_argument("--sizes needs at least one point count");
    }
    return sizes;
}

Options parse_options(int argc, char** argv) {
    Options options;
    bool custom_sizes = false;
    for (int i = 1; i < argc; ++i) {
        const std::string argument = argv[i];
        if (argument == "--help" || argument == "-h") {
            print_usage(argv[0]);
            std::exit(0);
        }
        if (argument == "--quick") {
            options.quick = true;
            continue;
        }
        if (argument == "--fresh") {
            options.reuse = false;
            continue;
        }
        if (argument == "--reuse") {
            options.reuse = true;
            continue;
        }
        if (argument == "--markdown") {
            options.format = OutputFormat::Markdown;
            continue;
        }
        if (argument == "--format") {
            if (++i == argc) {
                throw std::invalid_argument("--format needs a value");
            }
            const std::string format = argv[i];
            if (format == "text") {
                options.format = OutputFormat::Text;
            } else if (format == "markdown" || format == "md") {
                options.format = OutputFormat::Markdown;
            } else if (format == "csv") {
                options.format = OutputFormat::Csv;
            } else {
                throw std::invalid_argument(
                    "--format must be text, markdown, or csv");
            }
            continue;
        }
        if (argument == "--sizes") {
            if (++i == argc) {
                throw std::invalid_argument("--sizes needs a value");
            }
            options.sizes = parse_sizes(argv[i]);
            custom_sizes = true;
            continue;
        }
        throw std::invalid_argument("unknown option: " + argument);
    }
    if (options.quick && !custom_sizes) {
        options.sizes = {1000, 10000, 100000};
    }
    return options;
}

double median(std::vector<double> samples) {
    std::sort(samples.begin(), samples.end());
    const std::size_t middle = samples.size() / 2;
    if ((samples.size() & 1U) != 0) {
        return samples[middle];
    }
    return 0.5 * (samples[middle - 1] + samples[middle]);
}

template <typename Operation>
double measure_once(
    Operation&& operation,
    std::vector<Triangle>& output,
    std::uint64_t& checksum) {
    const auto start = std::chrono::steady_clock::now();
    output = operation();
    const auto end = std::chrono::steady_clock::now();
    checksum += output.size();
    if (!output.empty()) {
        checksum += output[output.size() / 2].i0;
    }
    return std::chrono::duration<double, std::milli>(end - start).count();
}

std::pair<std::size_t, std::size_t> repetitions(
    std::size_t point_count,
    bool quick) {
    if (quick) {
        return {1, 3};
    }
    if (point_count <= 10000) {
        return {3, 11};
    }
    if (point_count <= 100000) {
        return {2, 7};
    }
    return {1, 5};
}

void require_valid(
    const char* implementation,
    const std::vector<Point>& points,
    const std::vector<Triangle>& candidate) {
    std::string error;
    if (!benchmark_support::validate_mesh(points, candidate, error)) {
        throw std::runtime_error(
            std::string(implementation) + " validation failed: " + error);
    }
}

Row benchmark_case(
    Dataset dataset,
    std::size_t point_count,
    std::uint64_t seed,
    bool quick,
    bool reuse,
    Triangulator& serial,
    Triangulator& parallel,
    std::uint64_t& checksum) {
    std::vector<Point> points =
        benchmark_support::generate_points(
            dataset, point_count, seed, kBenchmarkDomain);

    std::vector<Triangle> serial_output =
        triangulate_points(serial, points);
    std::vector<Triangle> parallel_output =
        triangulate_points(parallel, points);
    require_valid("serial Delaunay32", points, serial_output);
    require_valid("automatic Delaunay32", points, parallel_output);
    if (!benchmark_support::meshes_equal(
            serial_output, parallel_output)) {
        throw std::runtime_error(
            "serial and automatic Delaunay32 meshes differ");
    }

    const auto [warmups, iterations] = repetitions(point_count, quick);
    const auto run_serial = [&] {
        if (reuse) {
            return triangulate_points(serial, points);
        }
        Triangulator one_shot;
        configure_triangulator(one_shot, 1);
        return triangulate_points(one_shot, points);
    };
    const auto run_parallel = [&] {
        if (reuse) {
            return triangulate_points(parallel, points);
        }
        Triangulator one_shot;
        configure_triangulator(one_shot, 0);
        return triangulate_points(one_shot, points);
    };
    for (std::size_t i = 0; i < warmups; ++i) {
        if ((i & 1U) == 0) {
            measure_once(run_serial, serial_output, checksum);
            measure_once(run_parallel, parallel_output, checksum);
        } else {
            measure_once(run_parallel, parallel_output, checksum);
            measure_once(run_serial, serial_output, checksum);
        }
    }

    std::vector<double> serial_samples;
    std::vector<double> parallel_samples;
    serial_samples.reserve(iterations);
    parallel_samples.reserve(iterations);
    for (std::size_t i = 0; i < iterations; ++i) {
        if ((i & 1U) == 0) {
            serial_samples.push_back(
                measure_once(run_serial, serial_output, checksum));
            parallel_samples.push_back(
                measure_once(run_parallel, parallel_output, checksum));
        } else {
            parallel_samples.push_back(
                measure_once(run_parallel, parallel_output, checksum));
            serial_samples.push_back(
                measure_once(run_serial, serial_output, checksum));
        }
    }

    return {
        dataset,
        point_count,
        serial_output.size(),
        automatic_thread_count(point_count),
        iterations,
        reuse,
        median(std::move(serial_samples)),
        median(std::move(parallel_samples)),
    };
}

double automatic_speedup(const Row& row) {
    return row.serial_ms / row.parallel_ms;
}

void print_preamble(const Options& options) {
    if (options.format == OutputFormat::Csv) {
        return;
    }
    const unsigned hardware_threads = std::thread::hardware_concurrency();
    const char* lifecycle =
        options.reuse ? "reused working storage" : "fresh instances";
    if (options.format == OutputFormat::Markdown) {
        std::cout
            << "Median end-to-end triangulation time in milliseconds. "
               "Generation and validation are excluded; Delaunay32 input "
               "ingestion and triangle materialization are included. "
               "Delaunay32 lifecycle: "
            << lifecycle
            << ". Automatic mode may "
               "fall back to one thread for small inputs. Detected hardware "
               "threads: "
            << hardware_threads << ". Coordinate domain: [0, "
            << kBenchmarkDomain << ") squared.\n\n";
        return;
    }
    std::cout
        << "Delaunay32 benchmark\n"
        << "median end-to-end time; generation and validation excluded\n"
        << "input ingestion and triangle materialization included\n"
        << "Delaunay32 lifecycle: " << lifecycle << '\n'
        << "coordinate domain: [0," << kBenchmarkDomain << ")^2\n"
        << "hardware threads: " << hardware_threads
        << " (auto falls back to serial for small inputs)\n\n";
#ifndef NDEBUG
    std::cout << "WARNING: assertions are enabled; use a Release build.\n\n";
#endif
}

void print_rows(const std::vector<Row>& rows, OutputFormat format) {
    if (format == OutputFormat::Csv) {
        std::cout
            << "distribution,points,lifecycle,delaunay32_1t_ms,"
               "delaunay32_auto_ms,auto_threads,auto_speedup,triangles,"
               "samples\n";
        for (const Row& row : rows) {
            std::cout << benchmark_support::dataset_name(row.dataset) << ','
                      << row.points << ','
                      << (row.reuse ? "reuse" : "fresh") << ','
                      << std::fixed
                      << std::setprecision(4) << row.serial_ms << ','
                      << row.parallel_ms << ',' << row.auto_threads << ','
                      << automatic_speedup(row) << ','
                      << row.triangles << ',' << row.samples << '\n';
        }
        return;
    }

    if (format == OutputFormat::Markdown) {
        std::cout
            << "| distribution | points | Delaunay32 1T ms | "
               "Delaunay32 auto ms | "
               "threads | auto speedup | triangles |\n"
            << "|:--|--:|--:|--:|--:|--:|--:|\n";
        for (const Row& row : rows) {
            std::cout << "| "
                      << benchmark_support::dataset_name(row.dataset)
                      << " | " << row.points << " | " << std::fixed
                      << std::setprecision(3) << row.serial_ms << " | "
                      << row.parallel_ms << " | " << row.auto_threads
                      << " | " << automatic_speedup(row)
                      << "x | " << row.triangles << " |\n";
        }
        return;
    }

    std::cout
        << std::left << std::setw(12) << "distribution"
        << std::right << std::setw(11) << "points"
        << std::setw(13) << "D32 1T"
        << std::setw(13) << "D32 auto"
        << std::setw(9) << "threads"
        << std::setw(13) << "speedup"
        << std::setw(13) << "triangles" << '\n'
        << std::string(84, '-') << '\n';
    for (const Row& row : rows) {
        std::cout
            << std::left << std::setw(12)
            << benchmark_support::dataset_name(row.dataset)
            << std::right << std::setw(11) << row.points
            << std::setw(13) << std::fixed << std::setprecision(3)
            << row.serial_ms
            << std::setw(13) << row.parallel_ms
            << std::setw(9) << row.auto_threads
            << std::setw(12) << automatic_speedup(row) << 'x'
            << std::setw(13) << row.triangles << '\n';
    }
    std::cout << std::string(84, '-') << '\n';
}

}  // namespace
}  // namespace delaunay32

int main(int argc, char** argv) {
    try {
        const delaunay32::Options options =
            delaunay32::parse_options(argc, argv);
        delaunay32::print_preamble(options);

        delaunay32::Triangulator serial;
        delaunay32::Triangulator parallel;
        delaunay32::configure_triangulator(serial, 1);
        delaunay32::configure_triangulator(parallel, 0);
        std::vector<delaunay32::Row> rows;
        rows.reserve(options.sizes.size() * 3);
        std::uint64_t checksum = 0;
        const std::vector<delaunay32::benchmark_support::Dataset> datasets = {
            delaunay32::benchmark_support::Dataset::Uniform,
            delaunay32::benchmark_support::Dataset::Clustered,
            delaunay32::benchmark_support::Dataset::Diagonal,
        };
        for (std::size_t dataset_index = 0;
             dataset_index < datasets.size();
             ++dataset_index) {
            for (std::size_t size_index = 0;
                 size_index < options.sizes.size();
                 ++size_index) {
                const std::uint64_t seed =
                    0x5eed1234ULL +
                    dataset_index * 0x10000ULL +
                    size_index * 97ULL;
                rows.push_back(delaunay32::benchmark_case(
                    datasets[dataset_index],
                    options.sizes[size_index],
                    seed,
                    options.quick,
                    options.reuse,
                    serial,
                    parallel,
                    checksum));
            }
        }
        delaunay32::print_rows(rows, options.format);
        if (checksum == 0) {
            std::cerr << "unexpected empty benchmark checksum\n";
            return 2;
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "benchmark error: " << error.what() << '\n';
        return 1;
    }
}
