// SPDX-License-Identifier: MIT

#include "delaunay32/extras/geometry.hpp"
#include "internal.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace delaunay32::extras {
namespace {

struct SignedMagnitude {
    bool negative = false;
    std::uint64_t magnitude = 0;
};

SignedMagnitude multiply(std::int64_t a, std::int64_t b) {
    const bool negative = (a < 0) != (b < 0);
    const std::uint64_t magnitude_a = static_cast<std::uint64_t>(
        a < 0 ? -a : a);
    const std::uint64_t magnitude_b = static_cast<std::uint64_t>(
        b < 0 ? -b : b);
    return {negative, magnitude_a * magnitude_b};
}

int orient_sign(const Point& a, const Point& b, const Point& point) {
    const SignedMagnitude lhs = multiply(
        static_cast<std::int64_t>(b.x) - a.x,
        static_cast<std::int64_t>(point.y) - a.y);
    const SignedMagnitude rhs = multiply(
        static_cast<std::int64_t>(b.y) - a.y,
        static_cast<std::int64_t>(point.x) - a.x);
    if (lhs.negative != rhs.negative) {
        return lhs.negative ? -1 : 1;
    }
    if (lhs.magnitude == rhs.magnitude) {
        return 0;
    }
    const int magnitude_order = lhs.magnitude > rhs.magnitude ? 1 : -1;
    return lhs.negative ? -magnitude_order : magnitude_order;
}

bool point_on_segment(const Point& point, const Point& a, const Point& b) {
    if (orient_sign(a, b, point) != 0) {
        return false;
    }
    return point.x >= std::min(a.x, b.x) &&
           point.x <= std::max(a.x, b.x) &&
           point.y >= std::min(a.y, b.y) &&
           point.y <= std::max(a.y, b.y);
}

enum class RingLocation {
    outside,
    inside,
    boundary,
};

void validate_ring(
    const std::vector<std::uint32_t>& ring,
    std::size_t point_count,
    const std::string& label) {
    if (ring.size() < 3) {
        throw std::invalid_argument(
            label + " must contain at least three indices");
    }
    for (const std::uint32_t index : ring) {
        if (index >= point_count) {
            throw std::invalid_argument(
                label + " index is outside the points array");
        }
    }
}

RingLocation locate_in_ring(
    const Point& point,
    const std::vector<std::uint32_t>& ring,
    const std::vector<Point>& points) {
    bool inside = false;
    for (std::size_t i = 0; i < ring.size(); ++i) {
        const Point& a = points[ring[i]];
        const Point& b = points[ring[(i + 1) % ring.size()]];
        if (point_on_segment(point, a, b)) {
            return RingLocation::boundary;
        }
        if ((a.y > point.y) != (b.y > point.y)) {
            const int side = orient_sign(a, b, point);
            const bool crosses_to_right =
                (b.y > a.y && side > 0) ||
                (b.y < a.y && side < 0);
            if (crosses_to_right) {
                inside = !inside;
            }
        }
    }
    return inside ? RingLocation::inside : RingLocation::outside;
}

long double orient(
    const FloatPoint& a,
    const FloatPoint& b,
    const FloatPoint& point) {
    return
        (static_cast<long double>(b.x) - a.x) *
            (static_cast<long double>(point.y) - a.y) -
        (static_cast<long double>(b.y) - a.y) *
            (static_cast<long double>(point.x) - a.x);
}

bool point_on_segment(
    const FloatPoint& point,
    const FloatPoint& a,
    const FloatPoint& b) {
    if (orient(a, b, point) != 0.0L) {
        return false;
    }
    return point.x >= std::min(a.x, b.x) &&
           point.x <= std::max(a.x, b.x) &&
           point.y >= std::min(a.y, b.y) &&
           point.y <= std::max(a.y, b.y);
}

RingLocation locate_in_ring(
    const FloatPoint& point,
    const std::vector<std::uint32_t>& ring,
    const std::vector<FloatPoint>& points) {
    bool inside = false;
    for (std::size_t i = 0; i < ring.size(); ++i) {
        const FloatPoint& a = points[ring[i]];
        const FloatPoint& b = points[ring[(i + 1) % ring.size()]];
        if (point_on_segment(point, a, b)) {
            return RingLocation::boundary;
        }
        if ((a.y > point.y) != (b.y > point.y)) {
            const long double side = orient(a, b, point);
            const bool crosses_to_right =
                (b.y > a.y && side > 0.0L) ||
                (b.y < a.y && side < 0.0L);
            if (crosses_to_right) {
                inside = !inside;
            }
        }
    }
    return inside ? RingLocation::inside : RingLocation::outside;
}

void validate_finite(const std::vector<FloatPoint>& points) {
    for (const FloatPoint& point : points) {
        if (!std::isfinite(point.x) || !std::isfinite(point.y)) {
            throw std::invalid_argument(
                "floating-point polygon coordinates must be finite");
        }
    }
}

}  // namespace

namespace detail {

void validate_domain(
    const PolygonDomain& domain,
    std::size_t point_count,
    const std::string& label) {
    validate_ring(
        domain.outer_ring, point_count, label + " outer ring");
    for (std::size_t i = 0; i < domain.holes.size(); ++i) {
        validate_ring(
            domain.holes[i],
            point_count,
            label + " hole " + std::to_string(i));
    }
}

bool point_is_strictly_inside_domain_unchecked(
    const Point& point,
    const PolygonDomain& domain,
    const std::vector<Point>& points) {
    if (locate_in_ring(point, domain.outer_ring, points) !=
        RingLocation::inside) {
        return false;
    }
    for (const std::vector<std::uint32_t>& hole : domain.holes) {
        if (locate_in_ring(point, hole, points) != RingLocation::outside) {
            return false;
        }
    }
    return true;
}

bool point_is_strictly_inside_domain_unchecked(
    const FloatPoint& point,
    const PolygonDomain& domain,
    const std::vector<FloatPoint>& points) {
    if (locate_in_ring(point, domain.outer_ring, points) !=
        RingLocation::inside) {
        return false;
    }
    for (const std::vector<std::uint32_t>& hole : domain.holes) {
        if (locate_in_ring(point, hole, points) != RingLocation::outside) {
            return false;
        }
    }
    return true;
}

}  // namespace detail

bool point_is_strictly_inside_domain(
    const Point& point,
    const PolygonDomain& domain,
    const std::vector<Point>& points) {
    detail::validate_domain(domain, points.size(), "polygon");
    return detail::point_is_strictly_inside_domain_unchecked(
        point, domain, points);
}

bool point_is_strictly_inside_domain(
    const FloatPoint& point,
    const PolygonDomain& domain,
    const std::vector<FloatPoint>& points) {
    if (!std::isfinite(point.x) || !std::isfinite(point.y)) {
        throw std::invalid_argument(
            "floating-point query coordinate must be finite");
    }
    validate_finite(points);
    detail::validate_domain(domain, points.size(), "polygon");
    return detail::point_is_strictly_inside_domain_unchecked(
        point, domain, points);
}

}  // namespace delaunay32::extras
