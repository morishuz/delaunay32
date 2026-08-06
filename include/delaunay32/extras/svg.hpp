// SPDX-License-Identifier: MIT

#pragma once

#include "delaunay32/delaunay.hpp"
#include "delaunay32/extras/geometry.hpp"
#include "delaunay32/quantization.hpp"

#include <cstdint>
#include <iosfwd>
#include <memory>
#include <string>
#include <vector>

namespace delaunay32::extras {

// Margins are ordered left, top, right, bottom and use the same canvas units
// as the Svg width and height.
struct SvgMargins {
    double left = 20.0;
    double top = 20.0;
    double right = 20.0;
    double bottom = 20.0;
};

struct SvgPointStyle {
    std::string fill = "#111719";
    std::string stroke = "none";
    double radius = 2.0;
    double stroke_width = 0.0;
};

struct SvgLineStyle {
    std::string stroke = "#344044";
    double stroke_width = 1.0;
};

struct SvgShapeStyle {
    std::string fill = "none";
    std::string stroke = "#344044";
    double stroke_width = 1.0;
};

struct SvgTriangleColorStyle {
    // Smaller palettes are accepted when the mesh permits them. Four unique
    // colors guarantee a result for a manifold triangle mesh.
    std::vector<std::string> palette = {
        "#d8eee7",
        "#d9e7f3",
        "#f1e2a8",
        "#e6def0",
        "#cce7ec",
        "#e5ebc8",
        "#f0d9cf",
        "#d8e5d5",
    };
    std::string stroke = "#344044";
    double stroke_width = 1.0;
};

struct SvgTextStyle {
    std::string fill = "#202426";
    std::string font_family = "system-ui, sans-serif";
    std::string font_weight = "normal";
    double font_size = 16.0;
};

// A retained SVG document. Geometry is supplied in world coordinates and is
// auto-fitted to the canvas by default. Style sizes and text coordinates are
// expressed in SVG canvas units, normally pixels.
class Svg {
public:
    Svg(double width, double height);
    ~Svg();

    Svg(const Svg& other);
    Svg& operator=(const Svg& other);
    Svg(Svg&& other) noexcept;
    Svg& operator=(Svg&& other) noexcept;

    double width() const noexcept;
    double height() const noexcept;

    // SVG is transparent unless a background is set.
    Svg& set_background(std::string color);
    Svg& clear_background();

    // Fits all recorded geometry uniformly inside the margined canvas and
    // centers it on any axis with unused space. Cartesian y coordinates are
    // flipped so that larger y values appear higher. Text, point radii, and
    // stroke widths do not affect the fitted geometry bounds.
    Svg& set_auto_fit(const SvgMargins& margins = {});

    // Maps world coordinates directly to SVG canvas coordinates:
    // (x * x_scale + x_offset, y * y_scale + y_offset).
    Svg& set_transform(
        double x_scale,
        double y_scale,
        double x_offset = 0.0,
        double y_offset = 0.0);

    Svg& draw_point(
        double x,
        double y,
        const SvgPointStyle& style = {});
    Svg& draw_point(
        const Point& point,
        const SvgPointStyle& style = {});
    Svg& draw_point(
        const FloatPoint& point,
        const SvgPointStyle& style = {});

    Svg& draw_line(
        double x0,
        double y0,
        double x1,
        double y1,
        const SvgLineStyle& style = {});
    Svg& draw_line(
        const Point& first,
        const Point& second,
        const SvgLineStyle& style = {});
    Svg& draw_line(
        const FloatPoint& first,
        const FloatPoint& second,
        const SvgLineStyle& style = {});

    Svg& draw_points(
        const std::vector<Point>& points,
        const SvgPointStyle& style = {});
    Svg& draw_points(
        const std::vector<FloatPoint>& points,
        const SvgPointStyle& style = {});

    // A ring's closing edge is implicit. The PolygonDomain overload draws its
    // outer ring and holes as one even-odd path.
    Svg& draw_polygon(
        const std::vector<Point>& points,
        const std::vector<std::uint32_t>& ring,
        const SvgShapeStyle& style = {});
    Svg& draw_polygon(
        const std::vector<FloatPoint>& points,
        const std::vector<std::uint32_t>& ring,
        const SvgShapeStyle& style = {});
    Svg& draw_polygon(
        const std::vector<Point>& points,
        const PolygonDomain& domain,
        const SvgShapeStyle& style = {});
    Svg& draw_polygon(
        const std::vector<FloatPoint>& points,
        const PolygonDomain& domain,
        const SvgShapeStyle& style = {});

    Svg& draw_triangles(
        const std::vector<Point>& points,
        const std::vector<Triangle>& triangles,
        const SvgShapeStyle& style = {});
    Svg& draw_triangles(
        const std::vector<FloatPoint>& points,
        const std::vector<Triangle>& triangles,
        const SvgShapeStyle& style = {});

    // Assigns different palette entries to triangles that share an indexed
    // edge. Triangles that meet only at one point are not adjacent.
    Svg& draw_colored_triangles(
        const std::vector<Point>& points,
        const std::vector<Triangle>& triangles,
        const SvgTriangleColorStyle& style = {});
    Svg& draw_colored_triangles(
        const std::vector<FloatPoint>& points,
        const std::vector<Triangle>& triangles,
        const SvgTriangleColorStyle& style = {});

    // Text positions are canvas coordinates and do not affect auto-fitting.
    Svg& draw_text(
        std::string text,
        double x,
        double y,
        const SvgTextStyle& style = {});

    std::string to_svg() const;
    void render_to_svg(const std::string& output_path) const;

private:
    void write_svg(std::ostream& output) const;

    struct Impl;
    std::unique_ptr<Impl> impl_;
};

// Opinionated one-call visualizations retained for simple diagnostics. Use
// Svg directly when layout or styling needs to be customized.
void write_mesh_svg(
    const std::string& output_path,
    const std::vector<Point>& points,
    const std::vector<Triangle>& triangles);

void write_mesh_svg(
    const std::string& output_path,
    const std::vector<FloatPoint>& points,
    const std::vector<Triangle>& triangles,
    const QuantizationReport& report);

void write_polygon_mesh_svg(
    const std::string& output_path,
    const std::vector<Point>& points,
    const PolygonDomain& domain,
    const std::vector<Triangle>& triangles);

}  // namespace delaunay32::extras
