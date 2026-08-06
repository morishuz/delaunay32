// SPDX-License-Identifier: MIT

#include "delaunay32/extras/json.hpp"
#include "internal.hpp"

#include <cctype>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <fstream>
#include <iterator>
#include <limits>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace delaunay32::extras {
namespace {

void validate_geometry(const Geometry& geometry);

class GeometryJsonParser {
public:
    GeometryJsonParser(const std::string& text, std::string input_path)
        : text_(text), input_path_(std::move(input_path)) {}

    Geometry parse() {
        Geometry geometry;
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
                    reject_duplicate(have_points, "points");
                    geometry.points = parse_points();
                    have_points = true;
                } else if (key == "constraints") {
                    reject_duplicate(have_constraints, "constraints");
                    geometry.constraints = parse_constraints();
                    have_constraints = true;
                } else if (key == "polygon") {
                    reject_duplicate(have_polygon, "polygon");
                    geometry.polygon = parse_polygon_domain();
                    have_polygon = true;
                } else if (key == "polygons") {
                    reject_duplicate(have_polygons, "polygons");
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
        if (have_polygons && geometry.polygons.empty()) {
            fail("polygons must contain at least one domain");
        }
        try {
            ::delaunay32::extras::validate_geometry(geometry);
        } catch (const std::invalid_argument& error) {
            fail(error.what());
        }
        return geometry;
    }

private:
    const std::string& text_;
    std::string input_path_;
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

    void reject_duplicate(bool present, const std::string& field) const {
        if (present) {
            fail("duplicate " + field + " field");
        }
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
                        code_point = code_point * 16U +
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
                    reject_duplicate(have_outer, "polygon outer");
                    domain.outer_ring = parse_ring();
                    have_outer = true;
                } else if (key == "holes") {
                    reject_duplicate(have_holes, "polygon holes");
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

};

void validate_geometry(const Geometry& geometry) {
    if (geometry.points.size() < 3) {
        throw std::invalid_argument(
            "points must contain at least three entries");
    }
    if (geometry.points.size() >
        static_cast<std::size_t>(
            std::numeric_limits<std::uint32_t>::max() >> 1U)) {
        throw std::invalid_argument("points contains too many entries");
    }
    for (const Constraint constraint : geometry.constraints) {
        if (constraint.i0 >= geometry.points.size() ||
            constraint.i1 >= geometry.points.size()) {
            throw std::invalid_argument(
                "constraint index is outside the points array");
        }
        if (constraint.i0 == constraint.i1) {
            throw std::invalid_argument(
                "constraint endpoints must use different indices");
        }
    }
    if (geometry.polygon.has_value() && !geometry.polygons.empty()) {
        throw std::invalid_argument(
            "polygon and polygons cannot both be present");
    }
    if (geometry.polygon.has_value()) {
        detail::validate_domain(
            *geometry.polygon, geometry.points.size(), "polygon");
    }
    for (std::size_t i = 0; i < geometry.polygons.size(); ++i) {
        detail::validate_domain(
            geometry.polygons[i],
            geometry.points.size(),
            "polygons domain " + std::to_string(i));
    }
}

void write_ring(
    std::ostream& output,
    const std::vector<std::uint32_t>& ring) {
    output << '[';
    for (std::size_t i = 0; i < ring.size(); ++i) {
        output << (i == 0 ? "" : ", ") << ring[i];
    }
    output << ']';
}

void write_domain(
    std::ostream& output,
    const PolygonDomain& domain,
    const std::string& indent) {
    output << "{\n" << indent << "  \"outer\": ";
    write_ring(output, domain.outer_ring);
    if (!domain.holes.empty()) {
        output << ",\n" << indent << "  \"holes\": [";
        for (std::size_t i = 0; i < domain.holes.size(); ++i) {
            output << (i == 0 ? "\n" : ",\n") << indent << "    ";
            write_ring(output, domain.holes[i]);
        }
        output << '\n' << indent << "  ]";
    }
    output << '\n' << indent << '}';
}

void write_geometry(std::ostream& output, const Geometry& geometry) {
    output << "{\n  \"points\": [";
    for (std::size_t i = 0; i < geometry.points.size(); ++i) {
        const Point& point = geometry.points[i];
        output << (i == 0 ? "\n" : ",\n")
               << "    [" << point.x << ", " << point.y << ']';
    }
    output << "\n  ]";

    if (!geometry.constraints.empty()) {
        output << ",\n  \"constraints\": [";
        for (std::size_t i = 0; i < geometry.constraints.size(); ++i) {
            const Constraint constraint = geometry.constraints[i];
            output << (i == 0 ? "\n" : ",\n")
                   << "    [" << constraint.i0 << ", "
                   << constraint.i1 << ']';
        }
        output << "\n  ]";
    }
    if (geometry.polygon.has_value()) {
        output << ",\n  \"polygon\": ";
        write_domain(output, *geometry.polygon, "  ");
    }
    if (!geometry.polygons.empty()) {
        output << ",\n  \"polygons\": [";
        for (std::size_t i = 0; i < geometry.polygons.size(); ++i) {
            output << (i == 0 ? "\n    " : ",\n    ");
            write_domain(output, geometry.polygons[i], "    ");
        }
        output << "\n  ]";
    }
    output << "\n}\n";
}

}  // namespace

Geometry read_geometry_json(const std::string& input_path) {
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

std::string geometry_to_json(const Geometry& geometry) {
    validate_geometry(geometry);
    std::ostringstream output;
    write_geometry(output, geometry);
    return output.str();
}

void write_geometry_json(
    const std::string& output_path,
    const Geometry& geometry) {
    validate_geometry(geometry);
    std::ofstream output(output_path, std::ios::binary);
    if (!output) {
        throw std::runtime_error("could not create JSON: " + output_path);
    }

    write_geometry(output, geometry);
    if (!output) {
        throw std::runtime_error(
            "failed while writing JSON: " + output_path);
    }
}

}  // namespace delaunay32::extras
