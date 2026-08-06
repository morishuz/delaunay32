// SPDX-License-Identifier: MIT

#include "delaunay32/extras/svg.hpp"
#include "internal.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <optional>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace delaunay32::extras {
namespace {

struct SvgCoordinate {
    double x = 0.0;
    double y = 0.0;
};

enum class SvgElementKind {
    Point,
    Line,
    Path,
    Text,
};

struct SvgElement {
    SvgElementKind kind = SvgElementKind::Point;
    std::vector<SvgCoordinate> points;
    std::vector<std::vector<SvgCoordinate>> rings;
    SvgCoordinate first;
    SvgCoordinate second;
    SvgPointStyle point_style;
    SvgLineStyle line_style;
    SvgShapeStyle shape_style;
    SvgTextStyle text_style;
    std::string text;
};

struct SvgTransform {
    double x_scale = 1.0;
    double y_scale = 1.0;
    double x_offset = 0.0;
    double y_offset = 0.0;

    SvgCoordinate map(const SvgCoordinate& point) const {
        const SvgCoordinate result = {
            point.x * x_scale + x_offset,
            point.y * y_scale + y_offset,
        };
        if (!std::isfinite(result.x) || !std::isfinite(result.y)) {
            throw std::invalid_argument(
                "SVG transform produced a non-finite coordinate");
        }
        return result;
    }
};

struct SvgBounds {
    double min_x = 0.0;
    double max_x = 0.0;
    double min_y = 0.0;
    double max_y = 0.0;
    bool empty = true;

    void include(const SvgCoordinate& point) {
        if (empty) {
            min_x = max_x = point.x;
            min_y = max_y = point.y;
            empty = false;
            return;
        }
        min_x = std::min(min_x, point.x);
        max_x = std::max(max_x, point.x);
        min_y = std::min(min_y, point.y);
        max_y = std::max(max_y, point.y);
    }
};

void require_finite(double value, const char* label) {
    if (!std::isfinite(value)) {
        throw std::invalid_argument(std::string(label) + " must be finite");
    }
}

void require_nonnegative(double value, const char* label) {
    require_finite(value, label);
    if (value < 0.0) {
        throw std::invalid_argument(
            std::string(label) + " must not be negative");
    }
}

void validate_point_style(const SvgPointStyle& style) {
    require_nonnegative(style.radius, "SVG point radius");
    require_nonnegative(style.stroke_width, "SVG point stroke width");
}

void validate_line_style(const SvgLineStyle& style) {
    require_nonnegative(style.stroke_width, "SVG line stroke width");
}

void validate_shape_style(const SvgShapeStyle& style) {
    require_nonnegative(style.stroke_width, "SVG shape stroke width");
}

void validate_triangle_color_style(const SvgTriangleColorStyle& style) {
    require_nonnegative(
        style.stroke_width, "SVG colored triangle stroke width");
    if (style.palette.empty()) {
        throw std::invalid_argument(
            "SVG triangle color palette must not be empty");
    }
    std::unordered_set<std::string> colors;
    colors.reserve(style.palette.size());
    for (const std::string& color : style.palette) {
        if (!colors.insert(color).second) {
            throw std::invalid_argument(
                "SVG triangle color palette must contain unique colors");
        }
    }
}

void validate_text_style(const SvgTextStyle& style) {
    require_nonnegative(style.font_size, "SVG font size");
}

SvgCoordinate coordinate(double x, double y) {
    require_finite(x, "SVG x coordinate");
    require_finite(y, "SVG y coordinate");
    return {x, y};
}

template <typename PointType>
SvgCoordinate coordinate(const PointType& point) {
    return coordinate(
        static_cast<double>(point.x), static_cast<double>(point.y));
}

template <typename PointType>
void validate_triangles(
    const std::vector<PointType>& points,
    const std::vector<Triangle>& triangles) {
    for (const Triangle& triangle : triangles) {
        if (triangle.i0 >= points.size() || triangle.i1 >= points.size() ||
            triangle.i2 >= points.size()) {
            throw std::invalid_argument(
                "triangle index is outside the points array");
        }
    }
}

std::uint64_t edge_key(std::uint32_t first, std::uint32_t second) {
    const std::uint32_t low = std::min(first, second);
    const std::uint32_t high = std::max(first, second);
    return (static_cast<std::uint64_t>(low) << 32U) |
           static_cast<std::uint64_t>(high);
}

struct EdgeOwner {
    std::size_t triangle = 0;
    bool paired = false;
};

struct TriangleNeighbors {
    std::array<std::size_t, 3> triangles{};
    std::size_t count = 0;

    void add(std::size_t triangle) {
        triangles[count] = triangle;
        ++count;
    }
};

template <typename PointType>
void draw_adjacency_colored_triangles(
    Svg& svg,
    const std::vector<PointType>& points,
    const std::vector<Triangle>& triangles,
    const SvgTriangleColorStyle& style) {
    validate_triangle_color_style(style);
    validate_triangles(points, triangles);
    if (triangles.empty()) {
        return;
    }

    std::vector<TriangleNeighbors> neighbors(triangles.size());
    std::unordered_map<std::uint64_t, EdgeOwner> edge_owners;
    edge_owners.reserve(triangles.size() * 3U);
    for (std::size_t triangle_index = 0;
         triangle_index < triangles.size();
         ++triangle_index) {
        const Triangle& triangle = triangles[triangle_index];
        if (triangle.i0 == triangle.i1 || triangle.i1 == triangle.i2 ||
            triangle.i2 == triangle.i0) {
            throw std::invalid_argument(
                "SVG adjacency coloring requires distinct triangle indices");
        }
        coordinate(points[triangle.i0]);
        coordinate(points[triangle.i1]);
        coordinate(points[triangle.i2]);
        const std::array<std::uint32_t, 3> vertices = {
            triangle.i0,
            triangle.i1,
            triangle.i2,
        };
        for (std::size_t edge = 0; edge < vertices.size(); ++edge) {
            const std::uint64_t key = edge_key(
                vertices[edge], vertices[(edge + 1U) % vertices.size()]);
            const auto [iterator, inserted] =
                edge_owners.emplace(key, EdgeOwner{triangle_index, false});
            if (inserted) {
                continue;
            }
            if (iterator->second.paired) {
                throw std::invalid_argument(
                    "SVG adjacency coloring requires a manifold triangle "
                    "mesh");
            }
            const std::size_t other = iterator->second.triangle;
            neighbors[triangle_index].add(other);
            neighbors[other].add(triangle_index);
            iterator->second.paired = true;
        }
    }

    const std::size_t unassigned = std::numeric_limits<std::size_t>::max();
    std::vector<std::size_t> assigned(triangles.size(), unassigned);
    std::vector<bool> unavailable(style.palette.size(), false);
    for (std::size_t triangle_index = 0;
         triangle_index < triangles.size();
         ++triangle_index) {
        std::fill(unavailable.begin(), unavailable.end(), false);
        const TriangleNeighbors& adjacent = neighbors[triangle_index];
        for (std::size_t i = 0; i < adjacent.count; ++i) {
            const std::size_t color = assigned[adjacent.triangles[i]];
            if (color != unassigned) {
                unavailable[color] = true;
            }
        }

        const Triangle& triangle = triangles[triangle_index];
        const std::size_t preferred =
            (static_cast<std::size_t>(triangle.i0) * 17U +
             static_cast<std::size_t>(triangle.i1) * 31U +
             static_cast<std::size_t>(triangle.i2) * 43U) %
            style.palette.size();
        for (std::size_t offset = 0; offset < style.palette.size(); ++offset) {
            const std::size_t candidate =
                (preferred + offset) % style.palette.size();
            if (!unavailable[candidate]) {
                assigned[triangle_index] = candidate;
                break;
            }
        }
        if (assigned[triangle_index] == unassigned) {
            throw std::invalid_argument(
                "SVG triangle color palette contains " +
                std::to_string(style.palette.size()) +
                " colors and cannot distinguish all edge-adjacent "
                "triangles; provide at least 4 unique colors to guarantee "
                "adjacency coloring for a manifold triangle mesh");
        }
    }

    std::vector<std::vector<Triangle>> groups(style.palette.size());
    for (std::size_t triangle_index = 0;
         triangle_index < triangles.size();
         ++triangle_index) {
        groups[assigned[triangle_index]].push_back(triangles[triangle_index]);
    }
    for (std::size_t color = 0; color < groups.size(); ++color) {
        SvgShapeStyle shape_style;
        shape_style.fill = style.palette[color];
        shape_style.stroke = style.stroke;
        shape_style.stroke_width = style.stroke_width;
        svg.draw_triangles(points, groups[color], shape_style);
    }
}

template <typename PointType>
std::vector<SvgCoordinate> make_ring(
    const std::vector<PointType>& points,
    const std::vector<std::uint32_t>& ring) {
    if (ring.size() < 3) {
        throw std::invalid_argument(
            "SVG polygon ring requires at least three indices");
    }
    std::vector<SvgCoordinate> result;
    result.reserve(ring.size());
    for (const std::uint32_t index : ring) {
        if (index >= points.size()) {
            throw std::invalid_argument(
                "polygon index is outside the points array");
        }
        result.push_back(coordinate(points[index]));
    }
    return result;
}

template <typename PointType>
std::vector<std::vector<SvgCoordinate>> make_domain_rings(
    const std::vector<PointType>& points,
    const PolygonDomain& domain) {
    detail::validate_domain(domain, points.size(), "SVG polygon");
    std::vector<std::vector<SvgCoordinate>> result;
    result.reserve(1 + domain.holes.size());
    result.push_back(make_ring(points, domain.outer_ring));
    for (const std::vector<std::uint32_t>& hole : domain.holes) {
        result.push_back(make_ring(points, hole));
    }
    return result;
}

std::string escape_xml(const std::string& value) {
    std::string escaped;
    escaped.reserve(value.size());
    for (const char character : value) {
        switch (character) {
        case '&':
            escaped += "&amp;";
            break;
        case '<':
            escaped += "&lt;";
            break;
        case '>':
            escaped += "&gt;";
            break;
        case '\"':
            escaped += "&quot;";
            break;
        case '\'':
            escaped += "&apos;";
            break;
        default:
            escaped += character;
            break;
        }
    }
    return escaped;
}

void write_coordinate(
    std::ostream& output,
    const SvgCoordinate& coordinate_value) {
    output << coordinate_value.x << ' ' << coordinate_value.y;
}

SvgBounds element_bounds(const std::vector<SvgElement>& elements) {
    SvgBounds bounds;
    for (const SvgElement& element : elements) {
        switch (element.kind) {
        case SvgElementKind::Point:
            for (const SvgCoordinate& point : element.points) {
                bounds.include(point);
            }
            break;
        case SvgElementKind::Line:
            bounds.include(element.first);
            bounds.include(element.second);
            break;
        case SvgElementKind::Path:
            for (const std::vector<SvgCoordinate>& ring : element.rings) {
                for (const SvgCoordinate& point : ring) {
                    bounds.include(point);
                }
            }
            break;
        case SvgElementKind::Text:
            break;
        }
    }
    return bounds;
}

SvgTransform auto_fit_transform(
    const SvgBounds& bounds,
    double width,
    double height,
    const SvgMargins& margins) {
    if (bounds.empty) {
        return {};
    }

    const double available_width = width - margins.left - margins.right;
    const double available_height = height - margins.top - margins.bottom;
    const double span_x = bounds.max_x - bounds.min_x;
    const double span_y = bounds.max_y - bounds.min_y;
    if (!std::isfinite(span_x) || !std::isfinite(span_y)) {
        throw std::invalid_argument(
            "SVG geometry range is too large to auto-fit");
    }

    double scale = 1.0;
    if (span_x > 0.0 || span_y > 0.0) {
        const double x_scale =
            span_x > 0.0
                ? available_width / span_x
                : std::numeric_limits<double>::max();
        const double y_scale =
            span_y > 0.0
                ? available_height / span_y
                : std::numeric_limits<double>::max();
        scale = std::min(x_scale, y_scale);
    }

    const double plot_width = span_x * scale;
    const double plot_height = span_y * scale;
    const double left =
        margins.left + (available_width - plot_width) / 2.0;
    const double top =
        margins.top + (available_height - plot_height) / 2.0;
    return {
        scale,
        -scale,
        left - bounds.min_x * scale,
        top + bounds.max_y * scale,
    };
}

template <typename PointType>
std::pair<double, double> diagnostic_canvas_size(
    const std::vector<PointType>& points) {
    constexpr double maximum_canvas_size = 1200.0;
    constexpr double minimum_canvas_width = 480.0;
    constexpr double header_height = 56.0;
    constexpr double plot_margin = 30.0;
    constexpr double maximum_plot_width =
        maximum_canvas_size - 2.0 * plot_margin;
    constexpr double maximum_plot_height =
        maximum_canvas_size - header_height - 2.0 * plot_margin;

    double min_x = static_cast<double>(points.front().x);
    double max_x = min_x;
    double min_y = static_cast<double>(points.front().y);
    double max_y = min_y;
    for (const PointType& point : points) {
        min_x = std::min(min_x, static_cast<double>(point.x));
        max_x = std::max(max_x, static_cast<double>(point.x));
        min_y = std::min(min_y, static_cast<double>(point.y));
        max_y = std::max(max_y, static_cast<double>(point.y));
    }
    const double span_x = max_x - min_x;
    const double span_y = max_y - min_y;
    double scale = 1.0;
    if (span_x > 0.0 || span_y > 0.0) {
        const double x_scale =
            span_x > 0.0
                ? maximum_plot_width / span_x
                : std::numeric_limits<double>::max();
        const double y_scale =
            span_y > 0.0
                ? maximum_plot_height / span_y
                : std::numeric_limits<double>::max();
        scale = std::min(x_scale, y_scale);
    }
    return {
        std::max(minimum_canvas_width, span_x * scale + 2.0 * plot_margin),
        span_y * scale + header_height + 2.0 * plot_margin,
    };
}

template <typename PointType>
void write_mesh_svg_impl(
    const std::string& output_path,
    const std::vector<PointType>& points,
    const std::vector<Triangle>& triangles,
    const QuantizationReport* report) {
    if (points.empty()) {
        throw std::invalid_argument("SVG output requires at least one point");
    }
    validate_triangles(points, triangles);
    const auto [width, height] = diagnostic_canvas_size(points);
    Svg svg(width, height);
    svg.set_background("#f7f7f5");
    svg.set_auto_fit({30.0, 86.0, 30.0, 30.0});

    std::ostringstream title;
    title << "Delaunay32" << (report != nullptr ? " float | " : " | ")
          << points.size() << " points | " << triangles.size()
          << " triangles";
    if (report != nullptr) {
        title << " | grid step " << std::scientific
              << std::setprecision(3) << report->grid_step
              << " | collapsed " << report->collapsed_points;
    }
    SvgTextStyle title_style;
    title_style.font_size = 18.0;
    svg.draw_text(title.str(), 30.0, 35.0, title_style);
    SvgTriangleColorStyle triangle_style;
    triangle_style.stroke = "#344044";
    triangle_style.stroke_width = 0.72;
    svg.draw_colored_triangles(points, triangles, triangle_style);

    SvgPointStyle point_style;
    point_style.radius =
        points.size() <= 2000 ? 2.2 : (points.size() <= 20000 ? 1.1 : 0.55);
    point_style.fill = "#111719";
    point_style.stroke = "#ffffff";
    point_style.stroke_width = 0.35;
    svg.draw_points(points, point_style);
    svg.render_to_svg(output_path);
}

}  // namespace

struct Svg::Impl {
    double width = 0.0;
    double height = 0.0;
    std::optional<std::string> background;
    bool auto_fit = true;
    SvgMargins margins;
    SvgTransform transform;
    std::vector<SvgElement> elements;
};

Svg::Svg(double width, double height) : impl_(std::make_unique<Impl>()) {
    require_finite(width, "SVG width");
    require_finite(height, "SVG height");
    if (width <= 0.0 || height <= 0.0) {
        throw std::invalid_argument("SVG dimensions must be positive");
    }
    impl_->width = width;
    impl_->height = height;
    set_auto_fit();
}

Svg::~Svg() = default;

Svg::Svg(const Svg& other) : impl_(std::make_unique<Impl>(*other.impl_)) {}

Svg& Svg::operator=(const Svg& other) {
    if (this != &other) {
        impl_ = std::make_unique<Impl>(*other.impl_);
    }
    return *this;
}

Svg::Svg(Svg&& other) noexcept = default;

Svg& Svg::operator=(Svg&& other) noexcept = default;

double Svg::width() const noexcept {
    return impl_->width;
}

double Svg::height() const noexcept {
    return impl_->height;
}

Svg& Svg::set_background(std::string color) {
    impl_->background = std::move(color);
    return *this;
}

Svg& Svg::clear_background() {
    impl_->background.reset();
    return *this;
}

Svg& Svg::set_auto_fit(const SvgMargins& margins) {
    require_nonnegative(margins.left, "SVG left margin");
    require_nonnegative(margins.top, "SVG top margin");
    require_nonnegative(margins.right, "SVG right margin");
    require_nonnegative(margins.bottom, "SVG bottom margin");
    if (margins.left + margins.right >= impl_->width ||
        margins.top + margins.bottom >= impl_->height) {
        throw std::invalid_argument(
            "SVG margins must leave a positive drawing area");
    }
    impl_->auto_fit = true;
    impl_->margins = margins;
    return *this;
}

Svg& Svg::set_transform(
    double x_scale,
    double y_scale,
    double x_offset,
    double y_offset) {
    require_finite(x_scale, "SVG x scale");
    require_finite(y_scale, "SVG y scale");
    require_finite(x_offset, "SVG x offset");
    require_finite(y_offset, "SVG y offset");
    if (x_scale == 0.0 || y_scale == 0.0) {
        throw std::invalid_argument("SVG scales must be non-zero");
    }
    impl_->auto_fit = false;
    impl_->transform = {x_scale, y_scale, x_offset, y_offset};
    return *this;
}

Svg& Svg::draw_point(
    double x,
    double y,
    const SvgPointStyle& style) {
    validate_point_style(style);
    SvgElement element;
    element.kind = SvgElementKind::Point;
    element.points.push_back(coordinate(x, y));
    element.point_style = style;
    impl_->elements.push_back(std::move(element));
    return *this;
}

Svg& Svg::draw_point(const Point& point, const SvgPointStyle& style) {
    return draw_point(
        static_cast<double>(point.x), static_cast<double>(point.y), style);
}

Svg& Svg::draw_point(
    const FloatPoint& point,
    const SvgPointStyle& style) {
    return draw_point(
        static_cast<double>(point.x), static_cast<double>(point.y), style);
}

Svg& Svg::draw_line(
    double x0,
    double y0,
    double x1,
    double y1,
    const SvgLineStyle& style) {
    validate_line_style(style);
    SvgElement element;
    element.kind = SvgElementKind::Line;
    element.first = coordinate(x0, y0);
    element.second = coordinate(x1, y1);
    element.line_style = style;
    impl_->elements.push_back(std::move(element));
    return *this;
}

Svg& Svg::draw_line(
    const Point& first,
    const Point& second,
    const SvgLineStyle& style) {
    return draw_line(
        static_cast<double>(first.x),
        static_cast<double>(first.y),
        static_cast<double>(second.x),
        static_cast<double>(second.y),
        style);
}

Svg& Svg::draw_line(
    const FloatPoint& first,
    const FloatPoint& second,
    const SvgLineStyle& style) {
    return draw_line(
        static_cast<double>(first.x),
        static_cast<double>(first.y),
        static_cast<double>(second.x),
        static_cast<double>(second.y),
        style);
}

Svg& Svg::draw_points(
    const std::vector<Point>& points,
    const SvgPointStyle& style) {
    validate_point_style(style);
    if (points.empty()) {
        return *this;
    }
    SvgElement element;
    element.kind = SvgElementKind::Point;
    element.points.reserve(points.size());
    for (const Point& point : points) {
        element.points.push_back(coordinate(point));
    }
    element.point_style = style;
    impl_->elements.push_back(std::move(element));
    return *this;
}

Svg& Svg::draw_points(
    const std::vector<FloatPoint>& points,
    const SvgPointStyle& style) {
    validate_point_style(style);
    if (points.empty()) {
        return *this;
    }
    SvgElement element;
    element.kind = SvgElementKind::Point;
    element.points.reserve(points.size());
    for (const FloatPoint& point : points) {
        element.points.push_back(coordinate(point));
    }
    element.point_style = style;
    impl_->elements.push_back(std::move(element));
    return *this;
}

Svg& Svg::draw_polygon(
    const std::vector<Point>& points,
    const std::vector<std::uint32_t>& ring,
    const SvgShapeStyle& style) {
    validate_shape_style(style);
    SvgElement element;
    element.kind = SvgElementKind::Path;
    element.rings.push_back(make_ring(points, ring));
    element.shape_style = style;
    impl_->elements.push_back(std::move(element));
    return *this;
}

Svg& Svg::draw_polygon(
    const std::vector<FloatPoint>& points,
    const std::vector<std::uint32_t>& ring,
    const SvgShapeStyle& style) {
    validate_shape_style(style);
    SvgElement element;
    element.kind = SvgElementKind::Path;
    element.rings.push_back(make_ring(points, ring));
    element.shape_style = style;
    impl_->elements.push_back(std::move(element));
    return *this;
}

Svg& Svg::draw_polygon(
    const std::vector<Point>& points,
    const PolygonDomain& domain,
    const SvgShapeStyle& style) {
    validate_shape_style(style);
    SvgElement element;
    element.kind = SvgElementKind::Path;
    element.rings = make_domain_rings(points, domain);
    element.shape_style = style;
    impl_->elements.push_back(std::move(element));
    return *this;
}

Svg& Svg::draw_polygon(
    const std::vector<FloatPoint>& points,
    const PolygonDomain& domain,
    const SvgShapeStyle& style) {
    validate_shape_style(style);
    SvgElement element;
    element.kind = SvgElementKind::Path;
    element.rings = make_domain_rings(points, domain);
    element.shape_style = style;
    impl_->elements.push_back(std::move(element));
    return *this;
}

Svg& Svg::draw_triangles(
    const std::vector<Point>& points,
    const std::vector<Triangle>& triangles,
    const SvgShapeStyle& style) {
    validate_shape_style(style);
    validate_triangles(points, triangles);
    if (triangles.empty()) {
        return *this;
    }
    SvgElement element;
    element.kind = SvgElementKind::Path;
    element.rings.reserve(triangles.size());
    for (const Triangle& triangle : triangles) {
        element.rings.push_back({
            coordinate(points[triangle.i0]),
            coordinate(points[triangle.i1]),
            coordinate(points[triangle.i2]),
        });
    }
    element.shape_style = style;
    impl_->elements.push_back(std::move(element));
    return *this;
}

Svg& Svg::draw_triangles(
    const std::vector<FloatPoint>& points,
    const std::vector<Triangle>& triangles,
    const SvgShapeStyle& style) {
    validate_shape_style(style);
    validate_triangles(points, triangles);
    if (triangles.empty()) {
        return *this;
    }
    SvgElement element;
    element.kind = SvgElementKind::Path;
    element.rings.reserve(triangles.size());
    for (const Triangle& triangle : triangles) {
        element.rings.push_back({
            coordinate(points[triangle.i0]),
            coordinate(points[triangle.i1]),
            coordinate(points[triangle.i2]),
        });
    }
    element.shape_style = style;
    impl_->elements.push_back(std::move(element));
    return *this;
}

Svg& Svg::draw_colored_triangles(
    const std::vector<Point>& points,
    const std::vector<Triangle>& triangles,
    const SvgTriangleColorStyle& style) {
    draw_adjacency_colored_triangles(*this, points, triangles, style);
    return *this;
}

Svg& Svg::draw_colored_triangles(
    const std::vector<FloatPoint>& points,
    const std::vector<Triangle>& triangles,
    const SvgTriangleColorStyle& style) {
    draw_adjacency_colored_triangles(*this, points, triangles, style);
    return *this;
}

Svg& Svg::draw_text(
    std::string text,
    double x,
    double y,
    const SvgTextStyle& style) {
    validate_text_style(style);
    SvgElement element;
    element.kind = SvgElementKind::Text;
    element.first = coordinate(x, y);
    element.text_style = style;
    element.text = std::move(text);
    impl_->elements.push_back(std::move(element));
    return *this;
}

void Svg::write_svg(std::ostream& output) const {
    const SvgTransform transform =
        impl_->auto_fit
            ? auto_fit_transform(
                  element_bounds(impl_->elements),
                  impl_->width,
                  impl_->height,
                  impl_->margins)
            : impl_->transform;

    output << std::fixed << std::setprecision(3)
           << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
           << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\""
           << impl_->width << "\" height=\"" << impl_->height
           << "\" viewBox=\"0 0 " << impl_->width << ' ' << impl_->height
           << "\">\n";
    if (impl_->background.has_value()) {
        output << "<rect width=\"100%\" height=\"100%\" fill=\""
               << escape_xml(*impl_->background) << "\"/>\n";
    }

    for (const SvgElement& element : impl_->elements) {
        switch (element.kind) {
        case SvgElementKind::Point:
            for (const SvgCoordinate& world_point : element.points) {
                const SvgCoordinate point = transform.map(world_point);
                output << "<circle cx=\"" << point.x << "\" cy=\""
                       << point.y << "\" r=\""
                       << element.point_style.radius << "\" fill=\""
                       << escape_xml(element.point_style.fill)
                       << "\" stroke=\""
                       << escape_xml(element.point_style.stroke)
                       << "\" stroke-width=\""
                       << element.point_style.stroke_width << "\"/>\n";
            }
            break;
        case SvgElementKind::Line: {
            const SvgCoordinate first = transform.map(element.first);
            const SvgCoordinate second = transform.map(element.second);
            output << "<line x1=\"" << first.x << "\" y1=\"" << first.y
                   << "\" x2=\"" << second.x << "\" y2=\"" << second.y
                   << "\" stroke=\""
                   << escape_xml(element.line_style.stroke)
                   << "\" stroke-width=\""
                   << element.line_style.stroke_width
                   << "\" stroke-linecap=\"round\"/>\n";
            break;
        }
        case SvgElementKind::Path:
            output << "<path d=\"";
            for (const std::vector<SvgCoordinate>& ring : element.rings) {
                output << "M ";
                for (std::size_t i = 0; i < ring.size(); ++i) {
                    if (i != 0) {
                        output << " L ";
                    }
                    write_coordinate(output, transform.map(ring[i]));
                }
                output << " Z ";
            }
            output << "\" fill=\"" << escape_xml(element.shape_style.fill)
                   << "\" fill-rule=\"evenodd\" stroke=\""
                   << escape_xml(element.shape_style.stroke)
                   << "\" stroke-width=\""
                   << element.shape_style.stroke_width
                   << "\" stroke-linejoin=\"round\"/>\n";
            break;
        case SvgElementKind::Text:
            output << "<text x=\"" << element.first.x << "\" y=\""
                   << element.first.y << "\" font-family=\""
                   << escape_xml(element.text_style.font_family)
                   << "\" font-size=\"" << element.text_style.font_size
                   << "\" font-weight=\""
                   << escape_xml(element.text_style.font_weight)
                   << "\" fill=\"" << escape_xml(element.text_style.fill)
                   << "\">" << escape_xml(element.text) << "</text>\n";
            break;
        }
    }
    output << "</svg>\n";
}

std::string Svg::to_svg() const {
    std::ostringstream output;
    write_svg(output);
    return output.str();
}

void Svg::render_to_svg(const std::string& output_path) const {
    std::ofstream output(output_path, std::ios::binary);
    if (!output) {
        throw std::runtime_error("could not create SVG: " + output_path);
    }
    write_svg(output);
    if (!output) {
        throw std::runtime_error(
            "failed while writing SVG: " + output_path);
    }
}

void write_mesh_svg(
    const std::string& output_path,
    const std::vector<Point>& points,
    const std::vector<Triangle>& triangles) {
    write_mesh_svg_impl(output_path, points, triangles, nullptr);
}

void write_mesh_svg(
    const std::string& output_path,
    const std::vector<FloatPoint>& points,
    const std::vector<Triangle>& triangles,
    const QuantizationReport& report) {
    write_mesh_svg_impl(output_path, points, triangles, &report);
}

void write_polygon_mesh_svg(
    const std::string& output_path,
    const std::vector<Point>& points,
    const PolygonDomain& domain,
    const std::vector<Triangle>& triangles) {
    if (points.empty()) {
        throw std::invalid_argument("SVG output requires at least one point");
    }
    validate_triangles(points, triangles);
    detail::validate_domain(domain, points.size(), "polygon");

    Svg svg(1040.0, 820.0);
    svg.set_background("#f2f5f3");
    svg.set_auto_fit({54.0, 126.0, 54.0, 54.0});

    SvgTextStyle title_style;
    title_style.fill = "#172326";
    title_style.font_size = 27.0;
    title_style.font_weight = "720";
    svg.draw_text(
        "Delaunay32 polygon triangulation", 54.0, 45.0, title_style);

    SvgTextStyle subtitle_style;
    subtitle_style.fill = "#5b686c";
    subtitle_style.font_size = 15.0;
    std::ostringstream subtitle;
    subtitle << points.size() << " input points | " << domain.holes.size()
             << " holes | " << triangles.size()
             << " domain triangles | hollow points are outside the domain";
    svg.draw_text(subtitle.str(), 54.0, 76.0, subtitle_style);

    SvgShapeStyle domain_fill;
    domain_fill.fill = "#ffffff";
    domain_fill.stroke = "none";
    domain_fill.stroke_width = 0.0;
    svg.draw_polygon(points, domain, domain_fill);
    SvgTriangleColorStyle triangle_style;
    triangle_style.stroke = "#58676b";
    triangle_style.stroke_width = 0.85;
    svg.draw_colored_triangles(points, triangles, triangle_style);

    SvgShapeStyle outer_style;
    outer_style.fill = "none";
    outer_style.stroke = "#172f35";
    outer_style.stroke_width = 5.2;
    svg.draw_polygon(points, domain.outer_ring, outer_style);
    SvgShapeStyle hole_style = outer_style;
    hole_style.stroke_width = 4.3;
    for (const std::vector<std::uint32_t>& hole : domain.holes) {
        svg.draw_polygon(points, hole, hole_style);
    }

    std::vector<bool> used(points.size(), false);
    for (const Triangle& triangle : triangles) {
        used[triangle.i0] = true;
        used[triangle.i1] = true;
        used[triangle.i2] = true;
    }
    for (std::size_t i = 0; i < points.size(); ++i) {
        SvgPointStyle style;
        style.radius = used[i] ? 3.3 : 4.0;
        style.fill = used[i] ? "#132125" : "#f2f5f3";
        style.stroke = used[i] ? "#ffffff" : "#b54b42";
        style.stroke_width = used[i] ? 1.0 : 2.0;
        svg.draw_point(points[i], style);
    }
    svg.render_to_svg(output_path);
}

}  // namespace delaunay32::extras
