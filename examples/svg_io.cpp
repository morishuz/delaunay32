// SPDX-License-Identifier: MIT

#include "svg_io.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

namespace delaunay32_example {
namespace {

using delaunay32::Point;
using delaunay32::Triangle;

constexpr std::int32_t kDomainMaximum = 999;
constexpr double kCanvasSize = 1200.0;
constexpr double kHeaderHeight = 56.0;
constexpr double kPlotMargin = 30.0;
constexpr double kPlotSize =
    kCanvasSize - kHeaderHeight - 2.0 * kPlotMargin;
constexpr std::size_t kInteriorCapacity =
    static_cast<std::size_t>(kDomainMaximum - 1) *
    static_cast<std::size_t>(kDomainMaximum - 1);
constexpr std::size_t kMaximumPoints = kInteriorCapacity + 4;

std::uint64_t point_key(const Point& point) {
    return (static_cast<std::uint64_t>(
                static_cast<std::uint32_t>(point.x))
            << 32U) |
           static_cast<std::uint32_t>(point.y);
}

std::vector<Point> generate_points(
    std::size_t point_count,
    std::uint64_t seed) {
    if (point_count < 4 || point_count > kMaximumPoints) {
        throw std::invalid_argument(
            "point count must be between 4 and " +
            std::to_string(kMaximumPoints));
    }

    std::vector<Point> points = {
        {0, 0},
        {kDomainMaximum, 0},
        {kDomainMaximum, kDomainMaximum},
        {0, kDomainMaximum},
    };
    points.reserve(point_count);

    std::unordered_set<std::uint64_t> occupied;
    occupied.reserve(point_count * 2);
    for (const Point& point : points) {
        occupied.insert(point_key(point));
    }

    std::mt19937_64 random(seed);
    std::uniform_int_distribution<std::int32_t> coordinate(
        1, kDomainMaximum - 1);
    while (points.size() < point_count) {
        const Point point = {coordinate(random), coordinate(random)};
        if (occupied.insert(point_key(point)).second) {
            points.push_back(point);
        }
    }
    return points;
}

std::string trim(const std::string& text) {
    const std::size_t first = text.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) {
        return {};
    }
    const std::size_t last = text.find_last_not_of(" \t\r\n");
    return text.substr(first, last - first + 1);
}

std::int32_t parse_coordinate(
    const std::string& value,
    const std::string& input_path,
    std::size_t line_number) {
    const std::string text = trim(value);
    std::size_t consumed = 0;
    long long parsed = 0;
    try {
        parsed = std::stoll(text, &consumed);
    } catch (const std::exception&) {
        throw std::invalid_argument(
            input_path + ":" + std::to_string(line_number) +
            ": coordinate is not a signed 32-bit integer");
    }
    if (consumed != text.size() ||
        parsed < std::numeric_limits<std::int32_t>::min() ||
        parsed > std::numeric_limits<std::int32_t>::max()) {
        throw std::invalid_argument(
            input_path + ":" + std::to_string(line_number) +
            ": coordinate is not a signed 32-bit integer");
    }
    return static_cast<std::int32_t>(parsed);
}

std::vector<Point> read_csv(const std::string& input_path) {
    std::ifstream input(input_path);
    if (!input) {
        throw std::runtime_error("could not open CSV: " + input_path);
    }

    std::vector<Point> points;
    std::string line;
    std::size_t line_number = 0;
    bool first_record = true;
    while (std::getline(input, line)) {
        ++line_number;
        const std::string record = trim(line);
        if (record.empty() || record[0] == '#') {
            continue;
        }

        const std::size_t comma = record.find(',');
        if (comma == std::string::npos ||
            record.find(',', comma + 1) != std::string::npos) {
            throw std::invalid_argument(
                input_path + ":" + std::to_string(line_number) +
                ": expected exactly two comma-separated fields");
        }

        const std::string x_text = trim(record.substr(0, comma));
        const std::string y_text = trim(record.substr(comma + 1));
        if (first_record && x_text == "x" && y_text == "y") {
            first_record = false;
            continue;
        }
        first_record = false;
        points.push_back({
            parse_coordinate(x_text, input_path, line_number),
            parse_coordinate(y_text, input_path, line_number),
        });
    }
    if (!input.eof()) {
        throw std::runtime_error("failed while reading CSV: " + input_path);
    }
    if (points.size() < 3) {
        throw std::invalid_argument(
            input_path + ": expected at least three point records");
    }
    if (points.size() >
        static_cast<std::size_t>(
            std::numeric_limits<std::uint32_t>::max() >> 1U)) {
        throw std::invalid_argument(input_path + ": too many point records");
    }
    return points;
}

std::size_t parse_point_count(const char* value) {
    std::size_t consumed = 0;
    const std::string text = value;
    const unsigned long long parsed = std::stoull(text, &consumed);
    if (consumed != text.size()) {
        throw std::invalid_argument("point count is not an integer");
    }
    if (parsed > std::numeric_limits<std::size_t>::max()) {
        throw std::invalid_argument("point count is too large");
    }
    return static_cast<std::size_t>(parsed);
}

std::uint64_t parse_seed(const char* value) {
    std::size_t consumed = 0;
    const std::string text = value;
    const std::uint64_t seed = std::stoull(text, &consumed);
    if (consumed != text.size()) {
        throw std::invalid_argument("seed is not an integer");
    }
    return seed;
}

struct SvgTransform {
    std::int32_t min_x = 0;
    std::int32_t max_x = 0;
    std::int32_t min_y = 0;
    std::int32_t max_y = 0;
    double scale = 1.0;
    double left = kPlotMargin;
    double top = kHeaderHeight + kPlotMargin;

    explicit SvgTransform(const std::vector<Point>& points) {
        min_x = max_x = points[0].x;
        min_y = max_y = points[0].y;
        for (const Point& point : points) {
            min_x = std::min(min_x, point.x);
            max_x = std::max(max_x, point.x);
            min_y = std::min(min_y, point.y);
            max_y = std::max(max_y, point.y);
        }
        const double span_x = static_cast<double>(
            static_cast<std::int64_t>(max_x) - min_x);
        const double span_y = static_cast<double>(
            static_cast<std::int64_t>(max_y) - min_y);
        const double longest_span = std::max(span_x, span_y);
        if (longest_span > 0.0) {
            scale = kPlotSize / longest_span;
        }
        const double width = span_x * scale;
        const double height = span_y * scale;
        left = (kCanvasSize - width) / 2.0;
        top = kHeaderHeight + kPlotMargin +
              (kPlotSize - height) / 2.0;
    }

    double x(std::int32_t value) const {
        return left +
               static_cast<double>(
                   static_cast<std::int64_t>(value) - min_x) *
                   scale;
    }

    double y(std::int32_t value) const {
        return top +
               static_cast<double>(
                   static_cast<std::int64_t>(max_y) - value) *
                   scale;
    }
};

void write_point(
    std::ostream& output,
    const Point& point,
    const SvgTransform& transform) {
    output << std::fixed << std::setprecision(2)
           << transform.x(point.x) << ',' << transform.y(point.y);
}

}  // namespace

Options parse_options(int argc, char** argv) {
    Options options;
    if (argc >= 2 && std::string(argv[1]) == "--input") {
        options.csv_mode = true;
        for (int i = 1; i < argc; ++i) {
            const std::string argument = argv[i];
            if (argument == "--input" && i + 1 < argc) {
                options.input_path = argv[++i];
            } else if (argument == "--output" && i + 1 < argc) {
                options.output_path = argv[++i];
            } else {
                throw std::invalid_argument(
                    "unknown or incomplete option: " + argument);
            }
        }
        if (options.input_path.empty()) {
            throw std::invalid_argument("--input requires a CSV path");
        }
        return options;
    }

    if (argc > 4) {
        throw std::invalid_argument("too many positional arguments");
    }
    options.point_count =
        argc >= 2 ? parse_point_count(argv[1]) : 1000;
    options.output_path = argc >= 3 ? argv[2] : "mesh.svg";
    options.seed = argc >= 4 ? parse_seed(argv[3]) : 1;
    return options;
}

void print_usage(const char* executable) {
    std::cerr
        << "Usage:\n"
        << "  " << executable
        << " [point-count] [output.svg] [seed]\n"
        << "  " << executable
        << " --input points.csv [--output output.svg]\n";
}

std::vector<Point> load_points(const Options& options) {
    return options.csv_mode
               ? read_csv(options.input_path)
               : generate_points(options.point_count, options.seed);
}

void write_svg(
    const std::string& output_path,
    const std::vector<Point>& points,
    const std::vector<Triangle>& triangles) {
    std::ofstream output(output_path);
    if (!output) {
        throw std::runtime_error("could not create SVG: " + output_path);
    }
    const SvgTransform transform(points);

    static constexpr std::array<const char*, 7> colors = {
        "#dcefe8",
        "#dbe8f4",
        "#f3e3a6",
        "#e8dff0",
        "#cde8ee",
        "#e7ebc5",
        "#f1d9cf",
    };

    output
        << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
        << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"1200\" "
           "height=\"1200\" viewBox=\"0 0 1200 1200\">\n"
        << "<rect width=\"1200\" height=\"1200\" fill=\"#f7f7f5\"/>\n"
        << "<text x=\"30\" y=\"35\" font-family=\"system-ui, sans-serif\" "
           "font-size=\"18\" fill=\"#202426\">Delaunay32 | "
        << points.size() << " points | " << triangles.size()
        << " triangles</text>\n";

    for (const Triangle& triangle : triangles) {
        const std::size_t color =
            (static_cast<std::size_t>(triangle.i0) * 17 +
             static_cast<std::size_t>(triangle.i1) * 31 +
             static_cast<std::size_t>(triangle.i2) * 43) %
            colors.size();
        output << "<polygon fill=\"" << colors[color]
               << "\" stroke=\"#344044\" stroke-width=\"0.72\" "
                  "stroke-linejoin=\"round\" points=\"";
        write_point(output, points[triangle.i0], transform);
        output << ' ';
        write_point(output, points[triangle.i1], transform);
        output << ' ';
        write_point(output, points[triangle.i2], transform);
        output << "\"/>\n";
    }

    const double radius =
        points.size() <= 2000 ? 2.2 : (points.size() <= 20000 ? 1.1 : 0.55);
    for (const Point& point : points) {
        output << "<circle cx=\"" << std::fixed << std::setprecision(2)
               << transform.x(point.x) << "\" cy=\""
               << transform.y(point.y)
               << "\" r=\"" << radius
               << "\" fill=\"#111719\" stroke=\"#ffffff\" "
                  "stroke-width=\"0.35\"/>\n";
    }
    output << "</svg>\n";
    if (!output) {
        throw std::runtime_error(
            "failed while writing SVG: " + output_path);
    }
}

}  // namespace delaunay32_example
