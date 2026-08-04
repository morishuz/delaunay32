// SPDX-License-Identifier: MIT

#include "geometry_io.hpp"

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

namespace delaunay32_example {
namespace {

using delaunay32::Constraint;
using delaunay32::Point;

constexpr std::int32_t kDomainMaximum = 999;
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

std::int64_t orient(const Point& a, const Point& b, const Point& point) {
    return (static_cast<std::int64_t>(b.x) - a.x) *
               (static_cast<std::int64_t>(point.y) - a.y) -
           (static_cast<std::int64_t>(b.y) - a.y) *
               (static_cast<std::int64_t>(point.x) - a.x);
}

bool point_on_segment(const Point& point, const Point& a, const Point& b) {
    if (orient(a, b, point) != 0) {
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
        if ((a.y > point.y) == (b.y > point.y)) {
            continue;
        }
        const long double intersection_x =
            static_cast<long double>(a.x) +
            (static_cast<long double>(point.y) - a.y) *
                (static_cast<long double>(b.x) - a.x) /
                (static_cast<long double>(b.y) - a.y);
        if (static_cast<long double>(point.x) < intersection_x) {
            inside = !inside;
        }
    }
    return inside ? RingLocation::inside : RingLocation::outside;
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

class GeometryJsonParser {
public:
    GeometryJsonParser(
        const std::string& text,
        const std::string& input_path)
        : text_(text), input_path_(input_path) {}

    GeometryInput parse() {
        GeometryInput geometry;
        bool have_points = false;
        bool have_constraints = false;
        bool have_polygon = false;
        bool have_polygons = false;

        expect('{');
        if (!take('}')) {
            while (true) {
                const std::string key = parse_string();
                expect(':');
                if (key == "points") {
                    if (have_points) {
                        fail("duplicate points field");
                    }
                    geometry.points = parse_points();
                    have_points = true;
                } else if (key == "constraints") {
                    if (have_constraints) {
                        fail("duplicate constraints field");
                    }
                    geometry.constraints = parse_constraints();
                    have_constraints = true;
                } else if (key == "polygon") {
                    if (have_polygon) {
                        fail("duplicate polygon field");
                    }
                    parse_polygon(geometry);
                    have_polygon = true;
                } else if (key == "polygons") {
                    if (have_polygons) {
                        fail("duplicate polygons field");
                    }
                    geometry.polygons = parse_polygons();
                    have_polygons = true;
                } else {
                    fail("unknown top-level field: " + key);
                }
                if (take('}')) {
                    break;
                }
                expect(',');
            }
        }
        skip_whitespace();
        if (position_ != text_.size()) {
            fail("unexpected data after the top-level object");
        }
        if (!have_points) {
            fail("missing points field");
        }
        if (have_polygon && have_polygons) {
            fail("polygon and polygons cannot both be present");
        }
        validate(geometry, have_polygon, have_polygons);
        return geometry;
    }

private:
    const std::string& text_;
    const std::string& input_path_;
    std::size_t position_ = 0;

    [[noreturn]] void fail(const std::string& message) const {
        std::size_t line = 1;
        std::size_t column = 1;
        for (std::size_t i = 0; i < position_ && i < text_.size(); ++i) {
            if (text_[i] == '\n') {
                ++line;
                column = 1;
            } else {
                ++column;
            }
        }
        throw std::invalid_argument(
            input_path_ + ":" + std::to_string(line) + ":" +
            std::to_string(column) + ": " + message);
    }

    void skip_whitespace() {
        while (position_ < text_.size() &&
               std::isspace(
                   static_cast<unsigned char>(text_[position_])) != 0) {
            ++position_;
        }
    }

    bool take(char expected) {
        skip_whitespace();
        if (position_ < text_.size() && text_[position_] == expected) {
            ++position_;
            return true;
        }
        return false;
    }

    void expect(char expected) {
        if (!take(expected)) {
            fail(std::string("expected '") + expected + "'");
        }
    }

    static int hexadecimal_digit(char value) {
        if (value >= '0' && value <= '9') {
            return value - '0';
        }
        if (value >= 'a' && value <= 'f') {
            return 10 + value - 'a';
        }
        if (value >= 'A' && value <= 'F') {
            return 10 + value - 'A';
        }
        return -1;
    }

    std::string parse_string() {
        skip_whitespace();
        if (position_ >= text_.size() || text_[position_] != '"') {
            fail("expected a JSON string");
        }
        ++position_;
        std::string result;
        while (position_ < text_.size()) {
            const char value = text_[position_++];
            if (value == '"') {
                return result;
            }
            if (static_cast<unsigned char>(value) < 0x20U) {
                fail("control character in JSON string");
            }
            if (value != '\\') {
                result.push_back(value);
                continue;
            }
            if (position_ >= text_.size()) {
                fail("unterminated JSON escape");
            }
            const char escaped = text_[position_++];
            switch (escaped) {
                case '"':
                case '\\':
                case '/':
                    result.push_back(escaped);
                    break;
                case 'b':
                    result.push_back('\b');
                    break;
                case 'f':
                    result.push_back('\f');
                    break;
                case 'n':
                    result.push_back('\n');
                    break;
                case 'r':
                    result.push_back('\r');
                    break;
                case 't':
                    result.push_back('\t');
                    break;
                case 'u': {
                    unsigned code_point = 0;
                    for (int digit = 0; digit < 4; ++digit) {
                        if (position_ >= text_.size()) {
                            fail("incomplete Unicode escape");
                        }
                        const int value_digit =
                            hexadecimal_digit(text_[position_++]);
                        if (value_digit < 0) {
                            fail("invalid Unicode escape");
                        }
                        code_point =
                            code_point * 16U +
                            static_cast<unsigned>(value_digit);
                    }
                    if (code_point > 0x7fU) {
                        fail("non-ASCII field name is not supported");
                    }
                    result.push_back(static_cast<char>(code_point));
                    break;
                }
                default:
                    fail("invalid JSON escape");
            }
        }
        fail("unterminated JSON string");
    }

    std::int64_t parse_integer() {
        skip_whitespace();
        const std::size_t start = position_;
        if (position_ < text_.size() && text_[position_] == '-') {
            ++position_;
        }
        if (position_ >= text_.size() ||
            std::isdigit(
                static_cast<unsigned char>(text_[position_])) == 0) {
            fail("expected an integer");
        }
        if (text_[position_] == '0') {
            ++position_;
            if (position_ < text_.size() &&
                std::isdigit(
                    static_cast<unsigned char>(text_[position_])) != 0) {
                fail("leading zero in integer");
            }
        } else {
            while (position_ < text_.size() &&
                   std::isdigit(
                       static_cast<unsigned char>(text_[position_])) != 0) {
                ++position_;
            }
        }
        if (position_ < text_.size() &&
            (text_[position_] == '.' || text_[position_] == 'e' ||
             text_[position_] == 'E')) {
            fail("geometry values must be integers");
        }

        const std::string token = text_.substr(start, position_ - start);
        std::size_t consumed = 0;
        long long value = 0;
        try {
            value = std::stoll(token, &consumed);
        } catch (const std::exception&) {
            fail("integer is outside the signed 64-bit range");
        }
        if (consumed != token.size()) {
            fail("invalid integer");
        }
        return static_cast<std::int64_t>(value);
    }

    std::int32_t parse_coordinate() {
        const std::int64_t value = parse_integer();
        if (value < std::numeric_limits<std::int32_t>::min() ||
            value > std::numeric_limits<std::int32_t>::max()) {
            fail("coordinate is outside the signed 32-bit range");
        }
        return static_cast<std::int32_t>(value);
    }

    std::uint32_t parse_index() {
        const std::int64_t value = parse_integer();
        if (value < 0 ||
            static_cast<std::uint64_t>(value) >
                std::numeric_limits<std::uint32_t>::max()) {
            fail("index is outside the unsigned 32-bit range");
        }
        return static_cast<std::uint32_t>(value);
    }

    std::vector<Point> parse_points() {
        std::vector<Point> points;
        expect('[');
        if (take(']')) {
            return points;
        }
        while (true) {
            expect('[');
            const std::int32_t x = parse_coordinate();
            expect(',');
            const std::int32_t y = parse_coordinate();
            expect(']');
            points.push_back({x, y});
            if (take(']')) {
                return points;
            }
            expect(',');
        }
    }

    std::vector<Constraint> parse_constraints() {
        std::vector<Constraint> constraints;
        expect('[');
        if (take(']')) {
            return constraints;
        }
        while (true) {
            expect('[');
            const std::uint32_t i0 = parse_index();
            expect(',');
            const std::uint32_t i1 = parse_index();
            expect(']');
            constraints.push_back({i0, i1});
            if (take(']')) {
                return constraints;
            }
            expect(',');
        }
    }

    std::vector<std::uint32_t> parse_ring() {
        std::vector<std::uint32_t> ring;
        expect('[');
        if (take(']')) {
            return ring;
        }
        while (true) {
            ring.push_back(parse_index());
            if (take(']')) {
                return ring;
            }
            expect(',');
        }
    }

    std::vector<std::vector<std::uint32_t>> parse_holes() {
        std::vector<std::vector<std::uint32_t>> holes;
        expect('[');
        if (take(']')) {
            return holes;
        }
        while (true) {
            holes.push_back(parse_ring());
            if (take(']')) {
                return holes;
            }
            expect(',');
        }
    }

    PolygonDomain parse_polygon_domain() {
        PolygonDomain domain;
        bool have_outer = false;
        bool have_holes = false;
        expect('{');
        if (!take('}')) {
            while (true) {
                const std::string key = parse_string();
                expect(':');
                if (key == "outer") {
                    if (have_outer) {
                        fail("duplicate polygon outer field");
                    }
                    domain.outer_ring = parse_ring();
                    have_outer = true;
                } else if (key == "holes") {
                    if (have_holes) {
                        fail("duplicate polygon holes field");
                    }
                    domain.holes = parse_holes();
                    have_holes = true;
                } else {
                    fail("unknown polygon field: " + key);
                }
                if (take('}')) {
                    break;
                }
                expect(',');
            }
        }
        if (!have_outer) {
            fail("polygon object is missing its outer field");
        }
        return domain;
    }

    void parse_polygon(GeometryInput& geometry) {
        PolygonDomain domain = parse_polygon_domain();
        geometry.outer_ring = std::move(domain.outer_ring);
        geometry.holes = std::move(domain.holes);
    }

    std::vector<PolygonDomain> parse_polygons() {
        std::vector<PolygonDomain> polygons;
        expect('[');
        if (take(']')) {
            return polygons;
        }
        while (true) {
            polygons.push_back(parse_polygon_domain());
            if (take(']')) {
                return polygons;
            }
            expect(',');
        }
    }

    void validate_ring(
        const std::vector<std::uint32_t>& ring,
        std::size_t point_count,
        const std::string& label) const {
        if (ring.size() < 3) {
            fail(label + " must contain at least three indices");
        }
        for (const std::uint32_t index : ring) {
            if (index >= point_count) {
                fail(label + " index is outside the points array");
            }
        }
    }

    void validate(
        const GeometryInput& geometry,
        bool have_polygon,
        bool have_polygons) const {
        if (geometry.points.size() < 3) {
            fail("points must contain at least three entries");
        }
        if (geometry.points.size() >
            static_cast<std::size_t>(
                std::numeric_limits<std::uint32_t>::max() >> 1U)) {
            fail("points contains too many entries");
        }
        for (const Constraint constraint : geometry.constraints) {
            if (constraint.i0 >= geometry.points.size() ||
                constraint.i1 >= geometry.points.size()) {
                fail("constraint index is outside the points array");
            }
            if (constraint.i0 == constraint.i1) {
                fail("constraint endpoints must use different indices");
            }
        }
        if (have_polygon) {
            validate_ring(
                geometry.outer_ring,
                geometry.points.size(),
                "polygon outer ring");
            for (std::size_t i = 0; i < geometry.holes.size(); ++i) {
                validate_ring(
                    geometry.holes[i],
                    geometry.points.size(),
                    "polygon hole " + std::to_string(i));
            }
        }
        if (have_polygons) {
            if (geometry.polygons.empty()) {
                fail("polygons must contain at least one domain");
            }
            for (std::size_t polygon = 0;
                 polygon < geometry.polygons.size();
                 ++polygon) {
                const PolygonDomain& domain = geometry.polygons[polygon];
                validate_ring(
                    domain.outer_ring,
                    geometry.points.size(),
                    "polygons domain " + std::to_string(polygon) +
                        " outer ring");
                for (std::size_t hole = 0;
                     hole < domain.holes.size();
                     ++hole) {
                    validate_ring(
                        domain.holes[hole],
                        geometry.points.size(),
                        "polygons domain " + std::to_string(polygon) +
                            " hole " + std::to_string(hole));
                }
            }
        }
    }
};

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

}  // namespace

bool point_is_strictly_inside_domain(
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

Options parse_options(int argc, char** argv) {
    Options options;
    if (argc >= 2 && std::string(argv[1]) == "--input") {
        options.input_mode = true;
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
            throw std::invalid_argument("--input requires a JSON path");
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
        << " --input geometry.json [--output output.svg]\n";
}

GeometryInput read_geometry_json(const std::string& input_path) {
    std::ifstream input(input_path, std::ios::binary);
    if (!input) {
        throw std::runtime_error("could not open JSON: " + input_path);
    }
    const std::string text{
        std::istreambuf_iterator<char>(input),
        std::istreambuf_iterator<char>(),
    };
    if (input.bad()) {
        throw std::runtime_error("failed while reading JSON: " + input_path);
    }
    return GeometryJsonParser(text, input_path).parse();
}

GeometryInput load_geometry(const Options& options) {
    if (options.input_mode) {
        return read_geometry_json(options.input_path);
    }
    GeometryInput geometry;
    geometry.points = generate_points(options.point_count, options.seed);
    return geometry;
}

}  // namespace delaunay32_example
