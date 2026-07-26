// SPDX-License-Identifier: MIT

#include "delaunator_adapter.hpp"
#include "support.hpp"

#include "delaunay32/delaunay.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
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
    double delaunator_ms = 0.0;
};

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
    const std::vector<Triangle>& candidate,
    const std::vector<Triangle>& reference) {
    std::string error;
    if (!benchmark_support::validate_against_reference(
            points, candidate, reference, error)) {
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
    DelaunatorBaseline& delaunator,
    std::uint64_t& checksum) {
    std::vector<Point> points =
        benchmark_support::generate_points(
            dataset, point_count, seed, kBenchmarkDomain);

    delaunator.prepare(points);
    std::vector<Triangle> reference =
        delaunator.triangulate_prepared();
    std::vector<Triangle> serial_output = serial.triangulate(points);
    Profile parallel_profile;
    std::vector<Triangle> parallel_output =
        parallel.triangulate_profiled(points, parallel_profile, false);
    require_valid("serial Delaunay32", points, serial_output, reference);
    require_valid("parallel Delaunay32", points, parallel_output, reference);

    const auto [warmups, iterations] = repetitions(point_count, quick);
    const auto run_serial = [&] {
        if (reuse) {
            return serial.triangulate(points);
        }
        Triangulator one_shot(1);
        return one_shot.triangulate(points);
    };
    const auto run_parallel = [&] {
        if (reuse) {
            return parallel.triangulate(points);
        }
        Triangulator one_shot(0);
        return one_shot.triangulate(points);
    };
    const auto run_delaunator = [&] {
        return delaunator.triangulate_prepared();
    };

    for (std::size_t i = 0; i < warmups; ++i) {
        if (i % 3 == 0) {
            measure_once(run_serial, serial_output, checksum);
            measure_once(run_parallel, parallel_output, checksum);
            measure_once(run_delaunator, reference, checksum);
        } else if (i % 3 == 1) {
            measure_once(run_parallel, parallel_output, checksum);
            measure_once(run_delaunator, reference, checksum);
            measure_once(run_serial, serial_output, checksum);
        } else {
            measure_once(run_delaunator, reference, checksum);
            measure_once(run_serial, serial_output, checksum);
            measure_once(run_parallel, parallel_output, checksum);
        }
    }

    std::vector<double> serial_samples;
    std::vector<double> parallel_samples;
    std::vector<double> delaunator_samples;
    serial_samples.reserve(iterations);
    parallel_samples.reserve(iterations);
    delaunator_samples.reserve(iterations);
    for (std::size_t i = 0; i < iterations; ++i) {
        if (i % 3 == 0) {
            serial_samples.push_back(
                measure_once(run_serial, serial_output, checksum));
            parallel_samples.push_back(
                measure_once(run_parallel, parallel_output, checksum));
            delaunator_samples.push_back(
                measure_once(run_delaunator, reference, checksum));
        } else if (i % 3 == 1) {
            parallel_samples.push_back(
                measure_once(run_parallel, parallel_output, checksum));
            delaunator_samples.push_back(
                measure_once(run_delaunator, reference, checksum));
            serial_samples.push_back(
                measure_once(run_serial, serial_output, checksum));
        } else {
            delaunator_samples.push_back(
                measure_once(run_delaunator, reference, checksum));
            serial_samples.push_back(
                measure_once(run_serial, serial_output, checksum));
            parallel_samples.push_back(
                measure_once(run_parallel, parallel_output, checksum));
        }
    }

    return {
        dataset,
        point_count,
        serial_output.size(),
        parallel_profile.thread_count,
        iterations,
        reuse,
        median(std::move(serial_samples)),
        median(std::move(parallel_samples)),
        median(std::move(delaunator_samples)),
    };
}

double ratio(double numerator, double denominator) {
    return numerator / denominator;
}

double geometric_mean(
    const std::vector<Row>& rows,
    bool parallel) {
    double logarithm_sum = 0.0;
    for (const Row& row : rows) {
        logarithm_sum += std::log(
            ratio(
                parallel ? row.parallel_ms : row.serial_ms,
                row.delaunator_ms));
    }
    return std::exp(logarithm_sum / static_cast<double>(rows.size()));
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
               "ingestion and both implementations' triangle materialization "
               "are included. Delaunator's double coordinates are prepared "
               "outside timing. Delaunay32 lifecycle: "
            << lifecycle
            << ". Lower ratios are better. Automatic mode may "
               "fall back to one thread for small inputs. Detected hardware "
               "threads: "
            << hardware_threads << ". Coordinate domain: [0, "
            << kBenchmarkDomain << ") squared.\n\n";
        return;
    }
    std::cout
        << "Delaunay32 benchmark against Delaunator\n"
        << "median end-to-end time; generation and validation excluded\n"
        << "Delaunay32 ingestion and both output materializations included\n"
        << "Delaunator input conversion excluded\n"
        << "Delaunay32 lifecycle: " << lifecycle << '\n'
        << "coordinate domain: [0," << kBenchmarkDomain << ")^2\n"
        << "hardware threads: " << hardware_threads
        << " (auto falls back to serial for small inputs)\n\n";
#ifndef NDEBUG
    std::cout << "WARNING: assertions are enabled; use a Release build.\n\n";
#endif
}

void print_rows(const std::vector<Row>& rows, OutputFormat format) {
    const double serial_geomean = geometric_mean(rows, false);
    const double parallel_geomean = geometric_mean(rows, true);
    if (format == OutputFormat::Csv) {
        std::cout
            << "distribution,points,lifecycle,fast_1t_ms,fast_auto_ms,"
               "auto_threads,delaunator_ms,fast_1t_over_delaunator,"
               "fast_auto_over_delaunator,triangles,samples\n";
        for (const Row& row : rows) {
            std::cout << benchmark_support::dataset_name(row.dataset) << ','
                      << row.points << ','
                      << (row.reuse ? "reuse" : "fresh") << ','
                      << std::fixed
                      << std::setprecision(4) << row.serial_ms << ','
                      << row.parallel_ms << ',' << row.auto_threads << ','
                      << row.delaunator_ms << ','
                      << ratio(row.serial_ms, row.delaunator_ms) << ','
                      << ratio(row.parallel_ms, row.delaunator_ms) << ','
                      << row.triangles << ',' << row.samples << '\n';
        }
        return;
    }

    if (format == OutputFormat::Markdown) {
        std::cout
            << "| distribution | points | Delaunay32 1T ms | "
               "Delaunay32 auto ms | "
               "threads | Delaunator ms | 1T / del | auto / del | "
               "triangles |\n"
            << "|:--|--:|--:|--:|--:|--:|--:|--:|--:|\n";
        for (const Row& row : rows) {
            std::cout << "| "
                      << benchmark_support::dataset_name(row.dataset)
                      << " | " << row.points << " | " << std::fixed
                      << std::setprecision(3) << row.serial_ms << " | "
                      << row.parallel_ms << " | " << row.auto_threads
                      << " | " << row.delaunator_ms << " | "
                      << ratio(row.serial_ms, row.delaunator_ms) << "x | "
                      << ratio(row.parallel_ms, row.delaunator_ms)
                      << "x | " << row.triangles << " |\n";
        }
        std::cout << "\nGeometric mean: **"
                  << std::fixed << std::setprecision(3)
                  << serial_geomean << "x** (one thread), **"
                  << parallel_geomean
                  << "x** (automatic threads) versus Delaunator.\n";
        return;
    }

    std::cout
        << std::left << std::setw(12) << "distribution"
        << std::right << std::setw(11) << "points"
        << std::setw(13) << "D32 1T"
        << std::setw(13) << "D32 auto"
        << std::setw(9) << "threads"
        << std::setw(15) << "Delaunator"
        << std::setw(11) << "1T/del"
        << std::setw(11) << "auto/del"
        << std::setw(13) << "triangles" << '\n'
        << std::string(108, '-') << '\n';
    for (const Row& row : rows) {
        std::cout
            << std::left << std::setw(12)
            << benchmark_support::dataset_name(row.dataset)
            << std::right << std::setw(11) << row.points
            << std::setw(13) << std::fixed << std::setprecision(3)
            << row.serial_ms
            << std::setw(13) << row.parallel_ms
            << std::setw(9) << row.auto_threads
            << std::setw(15) << row.delaunator_ms
            << std::setw(10) << ratio(row.serial_ms, row.delaunator_ms)
            << 'x'
            << std::setw(10) << ratio(row.parallel_ms, row.delaunator_ms)
            << 'x'
            << std::setw(13) << row.triangles << '\n';
    }
    std::cout << std::string(108, '-') << '\n'
              << "geometric mean vs Delaunator: 1T "
              << std::fixed << std::setprecision(3)
              << serial_geomean << "x, auto " << parallel_geomean
              << "x (lower is faster)\n";
}

}  // namespace
}  // namespace delaunay32

int main(int argc, char** argv) {
    try {
        const delaunay32::Options options =
            delaunay32::parse_options(argc, argv);
        delaunay32::print_preamble(options);

        delaunay32::Triangulator serial(1);
        delaunay32::Triangulator parallel(0);
        delaunay32::DelaunatorBaseline delaunator;
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
                    delaunator,
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
