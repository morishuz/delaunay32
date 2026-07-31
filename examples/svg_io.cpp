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
constexpr double kMaximumCanvasSize = 1200.0;
constexpr double kMinimumCanvasWidth = 480.0;
constexpr double kHeaderHeight = 56.0;
constexpr double kPlotMargin = 30.0;
constexpr double kMaximumPlotWidth =
    kMaximumCanvasSize - 2.0 * kPlotMargin;
constexpr double kMaximumPlotHeight =
    kMaximumCanvasSize - kHeaderHeight - 2.0 * kPlotMargin;
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

template <typename PointType>
struct SvgTransform {
    double min_x = 0.0;
    double max_x = 0.0;
    double min_y = 0.0;
    double max_y = 0.0;
    double scale = 1.0;
    double left = kPlotMargin;
    double top = kHeaderHeight + kPlotMargin;
    double canvas_width = kMinimumCanvasWidth;
    double canvas_height = kHeaderHeight + 2.0 * kPlotMargin;

    explicit SvgTransform(const std::vector<PointType>& points) {
        min_x = max_x = points[0].x;
        min_y = max_y = points[0].y;
        for (const PointType& point : points) {
            min_x = std::min(min_x, static_cast<double>(point.x));
            max_x = std::max(max_x, static_cast<double>(point.x));
            min_y = std::min(min_y, static_cast<double>(point.y));
            max_y = std::max(max_y, static_cast<double>(point.y));
        }
        const double span_x = max_x - min_x;
        const double span_y = max_y - min_y;
        if (span_x > 0.0 || span_y > 0.0) {
            const double x_scale =
                span_x > 0.0
                    ? kMaximumPlotWidth / span_x
                    : std::numeric_limits<double>::max();
            const double y_scale =
                span_y > 0.0
                    ? kMaximumPlotHeight / span_y
                    : std::numeric_limits<double>::max();
            scale = std::min(x_scale, y_scale);
        }
        const double width = span_x * scale;
        const double height = span_y * scale;
        canvas_width =
            std::max(kMinimumCanvasWidth, width + 2.0 * kPlotMargin);
        canvas_height = height + kHeaderHeight + 2.0 * kPlotMargin;
        left = (canvas_width - width) / 2.0;
    }

    double x(double value) const {
        return left + (value - min_x) * scale;
    }

    double y(double value) const {
        return top + (max_y - value) * scale;
    }
};

template <typename PointType>
void write_point(
    std::ostream& output,
    const PointType& point,
    const SvgTransform<PointType>& transform) {
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

template <typename PointType>
void write_svg_impl(
    const std::string& output_path,
    const std::vector<PointType>& points,
    const std::vector<Triangle>& triangles,
    const delaunay32::QuantizationReport* report) {
    std::ofstream output(output_path);
    if (!output) {
        throw std::runtime_error("could not create SVG: " + output_path);
    }
    const SvgTransform<PointType> transform(points);

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
        << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\""
        << std::fixed << std::setprecision(2) << transform.canvas_width
        << "\" height=\"" << transform.canvas_height
        << "\" viewBox=\"0 0 " << transform.canvas_width << ' '
        << transform.canvas_height << "\">\n"
        << "<rect width=\"100%\" height=\"100%\" fill=\"#f7f7f5\"/>\n"
        << "<text x=\"30\" y=\"35\" font-family=\"system-ui, sans-serif\" "
           "font-size=\"18\" fill=\"#202426\">Delaunay32"
        << (report != nullptr ? " float | " : " | ")
        << points.size() << " points | " << triangles.size()
        << " triangles";
    if (report != nullptr) {
        output << " | grid step " << std::scientific << std::setprecision(3)
               << report->grid_step << " | collapsed "
               << report->collapsed_points;
    }
    output << "</text>\n";

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
    for (const PointType& point : points) {
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

void write_svg(
    const std::string& output_path,
    const std::vector<Point>& points,
    const std::vector<Triangle>& triangles) {
    write_svg_impl(output_path, points, triangles, nullptr);
}

void write_svg(
    const std::string& output_path,
    const std::vector<delaunay32::FloatPoint>& points,
    const std::vector<Triangle>& triangles,
    const delaunay32::QuantizationReport& report) {
    write_svg_impl(output_path, points, triangles, &report);
}

}  // namespace delaunay32_example
