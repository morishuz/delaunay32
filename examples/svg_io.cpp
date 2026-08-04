// SPDX-License-Identifier: MIT

#include "svg_io.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace delaunay32_example {
namespace {

using delaunay32::Point;
using delaunay32::Triangle;

constexpr double kMaximumCanvasSize = 1200.0;
constexpr double kMinimumCanvasWidth = 480.0;
constexpr double kHeaderHeight = 56.0;
constexpr double kPlotMargin = 30.0;
constexpr double kMaximumPlotWidth =
    kMaximumCanvasSize - 2.0 * kPlotMargin;
constexpr double kMaximumPlotHeight =
    kMaximumCanvasSize - kHeaderHeight - 2.0 * kPlotMargin;
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

template <typename PointType, typename Transform>
void write_point(
    std::ostream& output,
    const PointType& point,
    const Transform& transform) {
    output << std::fixed << std::setprecision(2)
           << transform.x(point.x) << ',' << transform.y(point.y);
}

}  // namespace

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

namespace {

struct ComparisonTransform {
    double min_x = 0.0;
    double max_x = 0.0;
    double min_y = 0.0;
    double max_y = 0.0;
    double scale = 1.0;
    double left = 0.0;
    double top = 0.0;

    ComparisonTransform(
        const std::vector<Point>& points,
        double panel_left,
        double panel_top,
        double panel_width,
        double panel_height) {
        min_x = max_x = points[0].x;
        min_y = max_y = points[0].y;
        for (const Point& point : points) {
            min_x = std::min(min_x, static_cast<double>(point.x));
            max_x = std::max(max_x, static_cast<double>(point.x));
            min_y = std::min(min_y, static_cast<double>(point.y));
            max_y = std::max(max_y, static_cast<double>(point.y));
        }
        const double span_x = max_x - min_x;
        const double span_y = max_y - min_y;
        const double x_scale =
            span_x == 0.0
                ? std::numeric_limits<double>::max()
                : panel_width / span_x;
        const double y_scale =
            span_y == 0.0
                ? std::numeric_limits<double>::max()
                : panel_height / span_y;
        scale = std::min(x_scale, y_scale);
        const double plot_width = span_x * scale;
        const double plot_height = span_y * scale;
        left = panel_left + (panel_width - plot_width) / 2.0;
        top = panel_top + (panel_height - plot_height) / 2.0;
    }

    double x(double value) const {
        return left + (value - min_x) * scale;
    }

    double y(double value) const {
        return top + (max_y - value) * scale;
    }
};

void write_comparison_panel(
    std::ostream& output,
    const std::vector<Point>& points,
    const std::vector<delaunay32::Constraint>& constraints,
    const std::vector<Triangle>& triangles,
    const ComparisonTransform& transform,
    double panel_left,
    const char* title,
    const char* subtitle,
    bool constrained) {
    static constexpr std::array<const char*, 6> colors = {
        "#dcefe8",
        "#dbe8f4",
        "#f3e3a6",
        "#e8dff0",
        "#cde8ee",
        "#f1d9cf",
    };

    output
        << "<rect x=\"" << panel_left
        << "\" y=\"92\" width=\"596\" height=\"636\" rx=\"14\" "
           "fill=\"#ffffff\" stroke=\"#d7dddf\"/>\n"
        << "<text x=\"" << panel_left + 24.0
        << "\" y=\"126\" font-family=\"system-ui, sans-serif\" "
           "font-size=\"20\" font-weight=\"650\" fill=\"#202426\">"
        << title << "</text>\n"
        << "<text x=\"" << panel_left + 24.0
        << "\" y=\"149\" font-family=\"system-ui, sans-serif\" "
           "font-size=\"13\" fill=\"#647075\">"
        << subtitle << "</text>\n";

    for (const Triangle& triangle : triangles) {
        const std::size_t color =
            (static_cast<std::size_t>(triangle.i0) * 17 +
             static_cast<std::size_t>(triangle.i1) * 31 +
             static_cast<std::size_t>(triangle.i2) * 43) %
            colors.size();
        output << "<polygon fill=\""
               << (constrained ? colors[color] : "#edf1f2")
               << "\" stroke=\"#536166\" stroke-width=\"0.8\" "
                  "stroke-linejoin=\"round\" points=\"";
        write_point(output, points[triangle.i0], transform);
        output << ' ';
        write_point(output, points[triangle.i1], transform);
        output << ' ';
        write_point(output, points[triangle.i2], transform);
        output << "\"/>\n";
    }

    for (const delaunay32::Constraint constraint : constraints) {
        const Point& a = points[constraint.i0];
        const Point& b = points[constraint.i1];
        output << "<line x1=\"" << std::fixed << std::setprecision(2)
               << transform.x(a.x) << "\" y1=\"" << transform.y(a.y)
               << "\" x2=\"" << transform.x(b.x) << "\" y2=\""
               << transform.y(b.y) << "\" stroke=\"#e0442e\" "
               << "stroke-width=\"" << (constrained ? 4.6 : 3.0)
               << "\" stroke-linecap=\"round\" opacity=\""
               << (constrained ? 0.96 : 0.58) << "\"";
        if (!constrained) {
            output << " stroke-dasharray=\"9 7\"";
        }
        output << "/>\n";
    }

    for (std::size_t i = 0; i < points.size(); ++i) {
        const Point& point = points[i];
        const double x = transform.x(point.x);
        const double y = transform.y(point.y);
        output << "<circle cx=\"" << x << "\" cy=\"" << y
               << "\" r=\"3.4\" fill=\"#111719\" stroke=\"#ffffff\" "
                  "stroke-width=\"1.2\"/>\n"
               << "<text x=\"" << x + 5.0 << "\" y=\"" << y - 5.0
               << "\" font-family=\"ui-monospace, monospace\" "
                  "font-size=\"10\" fill=\"#30383b\">"
               << i << "</text>\n";
    }
}

}  // namespace

void write_constrained_comparison_svg(
    const std::string& output_path,
    const std::vector<Point>& points,
    const std::vector<delaunay32::Constraint>& constraints,
    const std::vector<Triangle>& ordinary_triangles,
    const std::vector<Triangle>& constrained_triangles) {
    constexpr double canvas_width = 1280.0;
    constexpr double canvas_height = 760.0;
    constexpr double left_panel = 24.0;
    constexpr double right_panel = 660.0;
    constexpr double plot_top = 172.0;
    constexpr double plot_width = 540.0;
    constexpr double plot_height = 520.0;

    std::ofstream output(output_path);
    if (!output) {
        throw std::runtime_error("could not create SVG: " + output_path);
    }
    const ComparisonTransform ordinary_transform(
        points, left_panel + 28.0, plot_top, plot_width, plot_height);
    const ComparisonTransform constrained_transform(
        points, right_panel + 28.0, plot_top, plot_width, plot_height);

    output
        << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
        << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\""
        << canvas_width << "\" height=\"" << canvas_height
        << "\" viewBox=\"0 0 " << canvas_width << ' ' << canvas_height
        << "\">\n"
        << "<rect width=\"100%\" height=\"100%\" fill=\"#f4f6f5\"/>\n"
        << "<text x=\"24\" y=\"38\" "
           "font-family=\"system-ui, sans-serif\" font-size=\"24\" "
           "font-weight=\"700\" fill=\"#172023\">"
           "Delaunay32 constrained triangulation</text>\n"
        << "<text x=\"24\" y=\"65\" "
           "font-family=\"system-ui, sans-serif\" font-size=\"14\" "
           "fill=\"#566267\">"
        << points.size() << " points · " << constraints.size()
        << " requested segments · red lines show the same constraints in "
           "both panels</text>\n";

    write_comparison_panel(
        output,
        points,
        constraints,
        ordinary_triangles,
        ordinary_transform,
        left_panel,
        "Ordinary Delaunay",
        "Dashed requests may cut across existing triangle edges",
        false);
    write_comparison_panel(
        output,
        points,
        constraints,
        constrained_triangles,
        constrained_transform,
        right_panel,
        "Constrained Delaunay",
        "Every solid request is a mesh edge or edge chain",
        true);
    output << "</svg>\n";
    if (!output) {
        throw std::runtime_error(
            "failed while writing SVG: " + output_path);
    }
}

void write_polygon_svg(
    const std::string& output_path,
    const std::vector<Point>& points,
    const std::vector<std::uint32_t>& outer_ring,
    const std::vector<std::vector<std::uint32_t>>& holes,
    const std::vector<Triangle>& triangles) {
    constexpr double canvas_width = 1040.0;
    constexpr double canvas_height = 820.0;
    constexpr double plot_left = 54.0;
    constexpr double plot_top = 126.0;
    constexpr double plot_width = 932.0;
    constexpr double plot_height = 640.0;

    std::ofstream output(output_path);
    if (!output) {
        throw std::runtime_error("could not create SVG: " + output_path);
    }
    const ComparisonTransform transform(
        points, plot_left, plot_top, plot_width, plot_height);
    const auto write_ring_path = [&](const std::vector<std::uint32_t>& ring) {
        output << 'M';
        for (std::size_t i = 0; i < ring.size(); ++i) {
            output << (i == 0 ? " " : " L ");
            const std::uint32_t index = ring[i];
            write_point(output, points[index], transform);
        }
        output << " Z ";
    };

    output
        << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
        << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\""
        << canvas_width << "\" height=\"" << canvas_height
        << "\" viewBox=\"0 0 " << canvas_width << ' ' << canvas_height
        << "\">\n"
        << "<defs><filter id=\"shadow\" x=\"-10%\" y=\"-10%\" "
           "width=\"120%\" height=\"125%\"><feDropShadow dx=\"0\" "
           "dy=\"7\" stdDeviation=\"9\" flood-color=\"#142226\" "
           "flood-opacity=\"0.16\"/></filter></defs>\n"
        << "<rect width=\"100%\" height=\"100%\" fill=\"#f2f5f3\"/>\n"
        << "<text x=\"54\" y=\"45\" font-family=\"system-ui, sans-serif\" "
           "font-size=\"27\" font-weight=\"720\" fill=\"#172326\">"
           "Delaunay32 polygon triangulation</text>\n"
        << "<text x=\"54\" y=\"76\" font-family=\"system-ui, sans-serif\" "
           "font-size=\"15\" fill=\"#5b686c\">"
        << points.size() << " input points · " << holes.size()
        << " holes · " << triangles.size()
        << " domain triangles · hollow points are outside the domain</text>\n"
        << "<path d=\"";
    write_ring_path(outer_ring);
    for (const std::vector<std::uint32_t>& hole : holes) {
        write_ring_path(hole);
    }
    output
        << "\" fill=\"#ffffff\" fill-rule=\"evenodd\" "
           "filter=\"url(#shadow)\"/>\n";

    static constexpr std::array<const char*, 8> colors = {
        "#d8eee7",
        "#d9e7f3",
        "#f1e2a8",
        "#e6def0",
        "#cce7ec",
        "#e5ebc8",
        "#f0d9cf",
        "#d8e5d5",
    };
    std::vector<bool> used(points.size(), false);
    for (const Triangle& triangle : triangles) {
        used[triangle.i0] = true;
        used[triangle.i1] = true;
        used[triangle.i2] = true;
        const std::size_t color =
            (static_cast<std::size_t>(triangle.i0) * 17 +
             static_cast<std::size_t>(triangle.i1) * 31 +
             static_cast<std::size_t>(triangle.i2) * 43) %
            colors.size();
        output << "<polygon fill=\"" << colors[color]
               << "\" stroke=\"#58676b\" stroke-width=\"0.85\" "
                  "stroke-linejoin=\"round\" points=\"";
        write_point(output, points[triangle.i0], transform);
        output << ' ';
        write_point(output, points[triangle.i1], transform);
        output << ' ';
        write_point(output, points[triangle.i2], transform);
        output << "\"/>\n";
    }

    const auto write_ring_outline = [&](const std::vector<std::uint32_t>& ring,
                                        const char* color,
                                        double width) {
        output << "<polyline fill=\"none\" stroke=\"" << color
               << "\" stroke-width=\"" << width
               << "\" stroke-linejoin=\"round\" stroke-linecap=\"round\" "
                  "points=\"";
        for (const std::uint32_t index : ring) {
            write_point(output, points[index], transform);
            output << ' ';
        }
        write_point(output, points[ring.front()], transform);
        output << "\"/>\n";
    };
    write_ring_outline(outer_ring, "#172f35", 5.2);
    for (const std::vector<std::uint32_t>& hole : holes) {
        write_ring_outline(hole, "#b54b42", 4.3);
    }

    for (std::size_t i = 0; i < points.size(); ++i) {
        const Point& point = points[i];
        output << "<circle cx=\"" << std::fixed << std::setprecision(2)
               << transform.x(point.x) << "\" cy=\""
               << transform.y(point.y) << "\" r=\""
               << (used[i] ? 3.3 : 4.0) << "\" fill=\""
               << (used[i] ? "#132125" : "#f2f5f3")
               << "\" stroke=\""
               << (used[i] ? "#ffffff" : "#b54b42")
               << "\" stroke-width=\"" << (used[i] ? 1.0 : 2.0)
               << "\"/>\n";
    }
    output << "</svg>\n";
    if (!output) {
        throw std::runtime_error(
            "failed while writing SVG: " + output_path);
    }
}

void write_logo_polygon_svg(
    const std::string& output_path,
    const std::vector<Point>& points,
    std::size_t interior_point_count,
    const std::vector<PolygonDomain>& domains,
    const std::vector<Triangle>& triangles) {
    std::ofstream output(output_path);
    if (!output) {
        throw std::runtime_error("could not create SVG: " + output_path);
    }
    const SvgTransform<Point> transform(points);

    std::size_t hole_count = 0;
    for (const PolygonDomain& domain : domains) {
        hole_count += domain.holes.size();
    }

    output
        << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
        << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\""
        << std::fixed << std::setprecision(2) << transform.canvas_width
        << "\" height=\"" << transform.canvas_height
        << "\" viewBox=\"0 0 " << transform.canvas_width << ' '
        << transform.canvas_height << "\">\n"
        << "<rect width=\"100%\" height=\"100%\" fill=\"#f7f7f5\"/>\n"
        << "<text x=\"30\" y=\"35\" font-family=\"system-ui, sans-serif\" "
           "font-size=\"18\" fill=\"#202426\">"
           "Delaunay32 constrained polygon mesh | "
        << interior_point_count << " blue-noise interior points | "
        << hole_count
        << " holes | " << triangles.size() << " triangles</text>\n";

    static constexpr std::array<const char*, 8> colors = {
        "#d8eee7",
        "#d9e7f3",
        "#f1e2a8",
        "#e6def0",
        "#cce7ec",
        "#e5ebc8",
        "#f0d9cf",
        "#d8e5d5",
    };
    for (const Triangle& triangle : triangles) {
        const std::size_t color =
            (static_cast<std::size_t>(triangle.i0) * 17 +
             static_cast<std::size_t>(triangle.i1) * 31 +
             static_cast<std::size_t>(triangle.i2) * 43) %
            colors.size();
        output << "<polygon fill=\"" << colors[color]
               << "\" stroke=\"#4d5b60\" stroke-width=\"0.36\" "
                  "stroke-linejoin=\"round\" points=\"";
        write_point(output, points[triangle.i0], transform);
        output << ' ';
        write_point(output, points[triangle.i1], transform);
        output << ' ';
        write_point(output, points[triangle.i2], transform);
        output << "\"/>\n";
    }

    const auto write_ring_outline = [&](
        const std::vector<std::uint32_t>& ring,
        const char* color,
        double width) {
        output << "<polyline fill=\"none\" stroke=\"" << color
               << "\" stroke-width=\"" << width
               << "\" stroke-linejoin=\"round\" stroke-linecap=\"round\" "
                  "points=\"";
        for (const std::uint32_t index : ring) {
            write_point(output, points[index], transform);
            output << ' ';
        }
        write_point(output, points[ring.front()], transform);
        output << "\"/>\n";
    };
    for (const PolygonDomain& domain : domains) {
        write_ring_outline(domain.outer_ring, "#172f35", 0.36);
        for (const std::vector<std::uint32_t>& hole : domain.holes) {
            write_ring_outline(hole, "#172f35", 0.36);
        }
    }

    output << "</svg>\n";
    if (!output) {
        throw std::runtime_error(
            "failed while writing SVG: " + output_path);
    }
}

}  // namespace delaunay32_example
