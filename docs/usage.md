# Delaunay32 usage guide

This guide describes the public API in
[`delaunay32/delaunay.hpp`](../include/delaunay32/delaunay.hpp), including the
different guarantees for integer and floating-point input.

## Core model

Every triangulation call returns a `std::vector<delaunay32::Triangle>`. A
`Triangle` contains three indices into the caller's original input:

```cpp
struct Triangle {
    std::uint32_t i0;
    std::uint32_t i1;
    std::uint32_t i2;
};
```

The library does not place coordinates inside the returned triangles. This is
important for float input: the library quantizes coordinates internally to
choose the connectivity, but the indices still address the unchanged source
floats.

## Choosing an input API

| Input layout | Coordinates | Call |
|:--|:--|:--|
| `std::vector<Point>` | signed 32-bit integers | `triangulate_int(points)` |
| separate arrays | signed 32-bit integers | `triangulate_int(xs, ys, count)` |
| `std::vector<FloatPoint>` | 32-bit floats | `triangulate_float(points)` |
| `std::vector<FloatPoint>` with mapping report | 32-bit floats | `triangulate_float(points, report)` |
| separate arrays | 32-bit floats | `triangulate_float(xs, ys, count)` |
| separate arrays with mapping report | 32-bit floats | `triangulate_float(xs, ys, count, report)` |

Use integer input when the coordinates are already discrete or when exact
Delaunay decisions on a chosen integer grid are required. Direct float input is
appropriate for most graphics, mapping, visualization, and general meshing
workloads where retaining the source vertices matters more than reproducing
the precise Delaunay edge set of the original floats.

## Public types

All public types are in the `delaunay32` namespace.

### `Point`

```cpp
struct Point {
    std::int32_t x;
    std::int32_t y;
};
```

`Point` supplies exact signed integer coordinates. Large absolute offsets are
allowed; predicate selection depends on the span of the input, not its position
relative to zero.

### `FloatPoint`

```cpp
struct FloatPoint {
    float x;
    float y;
};
```

`FloatPoint` supplies finite IEEE-style float coordinates. The input vector is
read but never modified. Quantization is automatic and affects connectivity
only; use the returned triangle indices to retrieve the original values.

### `Triangle`

`i0`, `i1`, and `i2` index the original input layout. When duplicate or
quantized-coincident points are collapsed, the retained vertex is the lowest
original input index at that grid coordinate.

Triangles are counterclockwise in the coordinates used for triangulation. That
means the original integer coordinates for `triangulate_int()` and the internal
quantized coordinates for `triangulate_float()`. In nearly degenerate float
geometry, restoring the original floats can produce a different winding.

### `Triangulator`

`Triangulator` owns reusable working storage and, when requested, a retained
worker team. It is movable but not copyable.

```cpp
delaunay32::Triangulator serial;       // one thread
delaunay32::Triangulator also_serial(1);
delaunay32::Triangulator automatic(0); // hardware thread count for large input
delaunay32::Triangulator fixed(4);     // request four threads for large input
```

Small inputs automatically use the serial path even when more threads are
requested. Reusing an instance avoids rebuilding its working storage and worker
team. Do not call the same instance concurrently; use separate instances for
independent concurrent calls.

Parallel topology construction uses a fixed arena so workers can safely retain
integer edge indices. If an unusual input exhausts its measured headroom, the
partial parallel topology is discarded and the call transparently restarts
with the growable serial allocator. This recovery can change execution time,
but not the documented geometric guarantees or accepted input domain.

The requested thread count can be inspected or changed between calls:

```cpp
triangulator.set_thread_count(0);
const std::size_t requested = triangulator.thread_count();
```

`thread_count()` reports the requested setting, not the number of workers used
for a particular input.

### `QuantizationReport`

`QuantizationReport` is an optional output report. It does not configure or
perform quantization, and passing it does not change the resulting mesh.

```cpp
delaunay32::QuantizationReport report;
const auto triangles = triangulator.triangulate_float(points, report);
```

Inspect its fields only after a successful triangulation call:

| Field | Meaning |
|:--|:--|
| `origin_x`, `origin_y` | Source-space origin subtracted before scaling; normally the input minima. |
| `scale` | Multiplier from source coordinates to the integer grid. |
| `grid_step` | Source-space distance represented by one integer grid unit; `1 / scale`. |
| `max_coordinate_error` | Largest measured absolute reconstruction error in either coordinate. |
| `unique_points` | Number of distinct grid coordinates triangulated. |
| `collapsed_points` | Input count minus `unique_points`, including exact duplicates and quantization collisions. |

### `PredicateWidth`

`PredicateWidth` describes which exact integer predicate path can support a
pair of coordinate spans:

- `Int64`: the faster 64-bit path;
- `Int128`: the wider path, available when the compiler provides `__int128`;
- `Unsupported`: the requested spans are outside the certified range.

Most callers do not need to select a predicate width. `triangulate_int()` checks
the input and chooses automatically. The public selector is useful for
validating a domain before loading or generating a large point set:

```cpp
const auto width = delaunay32::Triangulator::predicate_width_for_spans(
    x_span, y_span);
if (width == delaunay32::PredicateWidth::Unsupported) {
    // Rescale, tile, or reject the input before triangulation.
}
```

`int64_wide_intermediates_for_spans()` exposes a narrower implementation choice
within the supported wide path. It is primarily useful for diagnostics and
benchmarking; normal callers can ignore it.

## Integer triangulation

### Vector input

```cpp
#include <delaunay32/delaunay.hpp>

#include <vector>

std::vector<delaunay32::Point> points = {
    {0, 0},
    {100, 0},
    {100, 100},
    {0, 100},
    {48, 37},
};

delaunay32::Triangulator triangulator(0);
const std::vector<delaunay32::Triangle> triangles =
    triangulator.triangulate_int(points);
```

The vector is passed by const reference and is not modified.

### Separate coordinate arrays

Use the structure-of-arrays overload when x and y coordinates are already held
in separate contiguous buffers:

```cpp
std::vector<std::int32_t> xs = {0, 100, 100, 0, 48};
std::vector<std::int32_t> ys = {0, 0, 100, 100, 37};

const auto triangles = triangulator.triangulate_int(
    xs.data(), ys.data(), xs.size());
```

Both pointers must be non-null and address at least `point_count` values. The
two integer overloads have identical geometric behavior.

### Exactness and supported spans

Orientation and in-circle decisions are exact for every accepted integer input.
The relevant values are the coordinate spans:

```text
sx = max_x - min_x
sy = max_y - min_y
```

For equal x and y spans, the largest conservatively certified values are:

- `Triangulator::kFastCoordinateSpan`, currently 29,609, for the 64-bit path;
- `Triangulator::kMaxCoordinateSpan`, currently 1,940,470,527 when `__int128`
  is available.

Thin rectangular domains may support a much larger span on one axis. The exact
runtime selector uses `L = sx² + sy²` and the following bounds:

```text
Int64:
  L <= UINT32_MAX
  2*sx*sy <= INT64_MAX
  L*2*max(sx,sy) <= INT64_MAX
  L*6*sx*sy <= INT64_MAX

Int128:
  L <= UINT64_MAX
  2*sx*sy <= INT64_MAX
  L*2*max(sx,sy) <= INT128_MAX
  L*6*sx*sy <= INT128_MAX
```

On compilers without `__int128`, `kMaxCoordinateSpan` equals
`kFastCoordinateSpan`, and inputs requiring the wider path are rejected.

## Float triangulation

### Vector input

Pass floats directly; no caller-side conversion is needed:

```cpp
std::vector<delaunay32::FloatPoint> points = {
    {0.125F, 0.25F},
    {5.5F, 0.1F},
    {6.0F, 4.5F},
    {-1.0F, 5.0F},
};

delaunay32::Triangulator triangulator(0);
const auto triangles = triangulator.triangulate_float(points);

for (const auto& triangle : triangles) {
    const auto& a = points[triangle.i0];
    const auto& b = points[triangle.i1];
    const auto& c = points[triangle.i2];
    // a, b, and c retain the exact float values supplied by the caller.
}
```

Request a report when the mapping precision or grid collisions matter:

```cpp
delaunay32::QuantizationReport report;
const auto triangles = triangulator.triangulate_float(points, report);

if (report.collapsed_points != 0) {
    // Some input points shared a quantized grid coordinate.
}
```

### Separate coordinate arrays

The structure-of-arrays forms accept existing float buffers, with or without a
report:

```cpp
const auto triangles = triangulator.triangulate_float(
    xs.data(), ys.data(), xs.size());

delaunay32::QuantizationReport report;
const auto measured_triangles = triangulator.triangulate_float(
    xs.data(), ys.data(), xs.size(), report);
```

Both pointers must be non-null and address at least `point_count` values.

### How automatic quantization works

Finite floats across the full float range are accepted. The library translates
the input by its minimum bounds and applies one shared scale to both axes:

```text
span = max(max_x - min_x, max_y - min_y)
scale = Triangulator::kMaxCoordinateSpan / span
qx = round((x - min_x) * scale)
qy = round((y - min_y) * scale)
```

One shared scale preserves aspect ratio; the shorter axis is not independently
stretched to fill the grid. The quantized integers are triangulated using the
same exact predicates as direct integer input. On platforms without `__int128`,
the target grid automatically uses the smaller certified 64-bit span.

### What is preserved

- The input float values are never replaced or modified.
- Every triangle index addresses the caller's original input.
- Aspect ratio is preserved by the uniform x/y scale.
- Exact duplicates and quantization collisions retain their lowest input index.

### What can change

The connectivity is the exact Delaunay triangulation of the quantized integer
points, not necessarily the Delaunay triangulation of the original float
coordinates. Quantization can:

- merge points that are closer than the available grid resolution;
- select a different diagonal for nearly cocircular configurations;
- perturb oblique collinearity;
- change winding when the original geometry is nearly degenerate.

Axis-aligned collinearity is preserved. For most graphics, mapping,
visualization, and general meshing applications these differences are
negligible. Use an adaptive-exact triangulator when the precise Delaunay edge
set of the source floats is a requirement.

## Common output behavior

At least three input points and at least three unique triangulation coordinates
are required. Coincident coordinates are collapsed deterministically, retaining
the lowest original input index. Three or more unique but entirely collinear
points are valid input and return an empty triangle vector because no
two-dimensional faces exist.

Do not assume that every input point appears in a returned triangle: duplicates
or float points that collide on the quantization grid can be omitted in favor
of their representative.

For cocircular input, more than one Delaunay triangulation can be valid. Compare
meshes geometrically rather than assuming that another implementation must
choose the same diagonal.

## Invalid input and exceptions

The triangulation functions throw `std::invalid_argument` when:

- fewer than three input points are supplied;
- fewer than three unique integer or quantized coordinates remain;
- a structure-of-arrays pointer is null;
- a float coordinate is NaN or infinite;
- integer coordinate spans exceed the certified predicate range;
- the input count exceeds the internal index or edge-arena range.

As with other allocating C++ interfaces, allocation failures can propagate as
standard library exceptions.

## Complete examples

- [`delaunay_float_example.cpp`](../examples/delaunay_float_example.cpp)
  generates float points, requests a quantization report, and writes the mesh
  using the original float coordinates.
- [`delaunay_svg_example.cpp`](../examples/delaunay_svg_example.cpp) accepts
  integer CSV input or generates random integer points and writes an SVG.

Build both examples with:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target delaunay_float_example delaunay_svg_example
```
