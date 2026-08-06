// SPDX-License-Identifier: MIT

#include "delaunay32/quantization.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <unordered_set>

namespace delaunay32 {
namespace {

struct FloatBounds {
    double min_x = 0.0;
    double min_y = 0.0;
    double max_x = 0.0;
    double max_y = 0.0;
};

struct FloatQuantizer {
    double origin_x = 0.0;
    double origin_y = 0.0;
    double scale = 0.0;
    double target_span = 0.0;
    bool clamp_to_target = false;

    double map(double value, double origin) const {
        return (value - origin) * scale;
    }

    std::int32_t quantize_coordinate(double value, double origin) const {
        if (scale == 0.0) {
            return 0;
        }
        double mapped = map(value, origin);
        if (clamp_to_target) {
            mapped = std::clamp(mapped, 0.0, target_span);
            return static_cast<std::int32_t>(std::llround(mapped));
        }
        const double rounded = std::round(mapped);
        if (!std::isfinite(rounded) ||
            rounded < std::numeric_limits<std::int32_t>::min() ||
            rounded > std::numeric_limits<std::int32_t>::max()) {
            throw std::invalid_argument(
                "quantized coordinate is outside the int32 range");
        }
        return static_cast<std::int32_t>(rounded);
    }

    double error(
        double value,
        double origin,
        std::int32_t quantized) const {
        if (scale == 0.0) {
            return std::abs(value - origin);
        }
        return std::abs(
                   map(value, origin) -
                   static_cast<double>(quantized)) /
               scale;
    }
};

FloatBounds find_bounds(const std::vector<FloatPoint>& points) {
    if (points.empty()) {
        throw std::invalid_argument(
            "quantization requires at least one point");
    }
    const double first_x = points[0].x;
    const double first_y = points[0].y;
    if (!std::isfinite(first_x) || !std::isfinite(first_y)) {
        throw std::invalid_argument(
            "floating-point coordinates must be finite");
    }

    FloatBounds bounds{first_x, first_y, first_x, first_y};
    for (std::size_t i = 1; i < points.size(); ++i) {
        const double x = points[i].x;
        const double y = points[i].y;
        if (!std::isfinite(x) || !std::isfinite(y)) {
            throw std::invalid_argument(
                "floating-point coordinates must be finite");
        }
        bounds.min_x = std::min(bounds.min_x, x);
        bounds.min_y = std::min(bounds.min_y, y);
        bounds.max_x = std::max(bounds.max_x, x);
        bounds.max_y = std::max(bounds.max_y, y);
    }
    return bounds;
}

FloatQuantizer make_quantizer(
    const FloatBounds& bounds,
    const QuantizationOptions& options) {
    const double target_span =
        static_cast<double>(Triangulator::kMaxCoordinateSpan);
    if (!std::isfinite(options.max_coordinate_error) ||
        options.max_coordinate_error < 0.0) {
        throw std::invalid_argument(
            "maximum coordinate error must be finite and nonnegative");
    }
    if (options.collision_policy !=
            QuantizationCollisionPolicy::Allow &&
        options.collision_policy !=
            QuantizationCollisionPolicy::Reject) {
        throw std::invalid_argument(
            "unknown quantization collision policy");
    }

    switch (options.mode) {
        case QuantizationMode::Automatic: {
            const double maximum_span = std::max(
                bounds.max_x - bounds.min_x,
                bounds.max_y - bounds.min_y);
            if (!std::isfinite(maximum_span)) {
                throw std::invalid_argument(
                    "floating-point coordinate span is too large to "
                    "quantize");
            }
            return {
                bounds.min_x,
                bounds.min_y,
                maximum_span == 0.0
                    ? 0.0
                    : target_span / maximum_span,
                target_span,
                true,
            };
        }
        case QuantizationMode::GridStep: {
            if (!std::isfinite(options.grid_step) ||
                options.grid_step <= 0.0) {
                throw std::invalid_argument(
                    "quantization grid step must be finite and positive");
            }
            const double scale = 1.0 / options.grid_step;
            if (!std::isfinite(scale) || scale <= 0.0) {
                throw std::invalid_argument(
                    "quantization grid step produces an invalid scale");
            }
            return {
                bounds.min_x,
                bounds.min_y,
                scale,
                target_span,
                false,
            };
        }
        case QuantizationMode::FixedScale: {
            if (!std::isfinite(options.origin_x) ||
                !std::isfinite(options.origin_y)) {
                throw std::invalid_argument(
                    "fixed quantization origin must be finite");
            }
            if (!std::isfinite(options.scale) ||
                options.scale <= 0.0 ||
                !std::isfinite(1.0 / options.scale)) {
                throw std::invalid_argument(
                    "fixed quantization scale must be finite and positive");
            }
            return {
                options.origin_x,
                options.origin_y,
                options.scale,
                target_span,
                false,
            };
        }
    }
    throw std::invalid_argument("unknown quantization mode");
}

std::uint64_t point_key(const Point& point) {
    return (static_cast<std::uint64_t>(
                static_cast<std::uint32_t>(point.x))
            << 32U) |
           static_cast<std::uint32_t>(point.y);
}

}  // namespace

QuantizationResult quantize(
    const std::vector<FloatPoint>& points,
    const QuantizationOptions& options) {
    const FloatQuantizer quantizer =
        make_quantizer(find_bounds(points), options);

    QuantizationResult result;
    result.points.resize(points.size());
    result.report.origin_x = quantizer.origin_x;
    result.report.origin_y = quantizer.origin_y;
    result.report.scale = quantizer.scale;
    result.report.grid_step =
        quantizer.scale == 0.0 ? 0.0 : 1.0 / quantizer.scale;

    std::unordered_set<std::uint64_t> unique;
    unique.reserve(points.size());
    for (std::size_t i = 0; i < points.size(); ++i) {
        const double input_x = points[i].x;
        const double input_y = points[i].y;
        const std::int32_t x = quantizer.quantize_coordinate(
            input_x, quantizer.origin_x);
        const std::int32_t y = quantizer.quantize_coordinate(
            input_y, quantizer.origin_y);
        result.points[i] = {x, y};
        result.report.max_coordinate_error = std::max({
            result.report.max_coordinate_error,
            quantizer.error(input_x, quantizer.origin_x, x),
            quantizer.error(input_y, quantizer.origin_y, y),
        });
        unique.insert(point_key(result.points[i]));
    }

    result.report.unique_points = unique.size();
    result.report.collapsed_points = points.size() - unique.size();
    if (options.max_coordinate_error > 0.0 &&
        result.report.max_coordinate_error >
            options.max_coordinate_error) {
        throw std::invalid_argument(
            "quantization exceeds the requested maximum coordinate error");
    }
    if (options.collision_policy ==
            QuantizationCollisionPolicy::Reject &&
        result.report.collapsed_points != 0) {
        throw std::invalid_argument(
            "quantization produced coincident points");
    }
    return result;
}

}  // namespace delaunay32
