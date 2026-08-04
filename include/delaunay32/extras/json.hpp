// SPDX-License-Identifier: MIT

#pragma once

#include "delaunay32/extras/geometry.hpp"

#include <string>

namespace delaunay32::extras {

// Reads and writes Delaunay32's geometry schema. These functions are not a
// general-purpose JSON API: unknown fields and non-integer geometry values are
// rejected deliberately.
Geometry read_geometry_json(const std::string& input_path);
void write_geometry_json(
    const std::string& output_path,
    const Geometry& geometry);

}  // namespace delaunay32::extras
