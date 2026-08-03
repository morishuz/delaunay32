# Delaunay32 usage guide

This guide describes the public API in
[`delaunay32/delaunay.hpp`](../include/delaunay32/delaunay.hpp), including the
different guarantees for integer and floating-point input.

## Core model

The basic triangulation calls return a `std::vector<delaunay32::Triangle>`. A
`Triangle` contains three indices into the caller's original input:

```cpp
struct Triangle {
    std::uint32_t i0;
    std::uint32_t i1;
    std::uint32_t i2;
};
```

The library does not place coordinates inside the returned triangles. This is
important for floating-point input: the library quantizes coordinates
internally to choose the connectivity, but the indices still address the
unchanged source values.

The opt-in `triangulate_int_full()` and `triangulate_float_full()` calls return
the same faces plus halfedge adjacency, convex-hull indices, input
representative mappings, and operation metadata.

## Choosing an input API

| Input layout | Coordinates | Call |
|:--|:--|:--|
| `std::vector<Point>` | signed 32-bit integers | `triangulate_int(points)` |
| `std::vector<Point>` plus `std::vector<Constraint>` | integers | `triangulate_constrained_int(points, constraints)` |
| `std::vector<Point>` plus index rings | integers | `triangulate_polygon_int(points, outer, holes)` |
| `std::vector<FloatPoint>` | 32-bit floats | `triangulate_float(points)` |
| `std::vector<Point>` with full output | integers | `triangulate_int_full(points)` |
| `std::vector<FloatPoint>` with full output | floats | `triangulate_float_full(points)` |

Use integer input when the coordinates are already discrete or when exact
Delaunay decisions on a chosen integer grid are required. Direct float input is
appropriate for most graphics, mapping, visualization, and general meshing
workloads where retaining the source vertices matters more than reproducing
the precise Delaunay edge set of the original values.

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
read but never modified. Quantization affects connectivity only; use the
returned triangle indices to retrieve the original values.

### `Triangle`

`i0`, `i1`, and `i2` index the original input layout. When duplicate or
quantized-coincident points are collapsed, the retained vertex is the lowest
original input index at that grid coordinate.

Triangles are counterclockwise in the coordinates used for triangulation. That
means the original integer coordinates for `triangulate_int()` and
`triangulate_constrained_int()` or `triangulate_polygon_int()`, and the
internal quantized coordinates for `triangulate_float()`. In nearly degenerate
source geometry, restoring the original coordinates can produce a different
winding.

### `Constraint`

```cpp
struct Constraint {
    std::uint32_t i0;
    std::uint32_t i1;
};
```

`Constraint` identifies an undirected segment whose endpoints index the
original integer point vector. It is used by
`triangulate_constrained_int()` only.

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

### `QuantizationOptions`

`QuantizationOptions` configures float input. The default value uses
automatic maximum-resolution quantization.

| Field | Meaning |
|:--|:--|
| `mode` | `Automatic`, `GridStep`, or `FixedScale`. |
| `grid_step` | Positive source-space step used by `GridStep`. |
| `origin_x`, `origin_y` | Source-space origin used by `FixedScale`. |
| `scale` | Positive integer units per source unit used by `FixedScale`. |
| `max_coordinate_error` | Optional positive acceptance limit; zero disables the check. |
| `collision_policy` | `Allow` or `Reject`; rejection includes exact duplicates and quantization collisions. |

`GridStep` uses the input minima as its origin. `FixedScale` is the mode for a
stable mapping across batches because its origin and scale do not depend on
each batch's bounds. Every mode uses one scale for both axes.

### `QuantizationReport`

`QuantizationReport` describes the mapping used by a completed full float
triangulation.

```cpp
const auto result = triangulator.triangulate_float_full(points);
const delaunay32::QuantizationReport& report = result.quantization;
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

### `TriangulationResult`

The full-result API is opt-in, so existing triangle-only calls do not construct
the additional vectors:

```cpp
const delaunay32::TriangulationResult result =
    triangulator.triangulate_int_full(points);
```

`TriangulationResult` contains:

| Field | Meaning |
|:--|:--|
| `triangles` | The same faces returned by the corresponding triangle-only call. |
| `halfedges` | One flattened opposite-edge index per triangle edge, or `-1` for a boundary edge. |
| `hull` | Original input indices around the convex hull. |
| `representatives` | One retained representative index for every original input index. |
| `quantization` | The float quantization report; zero-initialized for integer input. |
| `predicate_width` | Exact predicate path selected for the completed operation. |
| `actual_thread_count` | Threads actually used after size thresholds and any serial recovery. |

For triangle `t`, flattened edges `3*t`, `3*t+1`, and `3*t+2` correspond to
`i0->i1`, `i1->i2`, and `i2->i0`. If `result.halfedges[e]` is nonnegative, its
value identifies the same geometric edge in the neighboring triangle with the
opposite direction. Halfedge indices are signed 64-bit values so `-1` remains
an unambiguous boundary sentinel.

For a nondegenerate mesh, `hull` is counterclockwise in the coordinates used
for triangulation and is rotated to start at its lowest original input index.
A collinear result has no triangles or halfedges and contains its two geometric
endpoints in `hull`.

Every `representatives[i]` is the lowest input index at the same integer or
quantized coordinate as input `i`. Unique points map to themselves. This makes
duplicate removal and floating-point quantization collisions explicit without
changing the indices stored in `triangles` or `hull`.

The full calls populate every field. Use the corresponding triangle-only call
when adjacency, hull, representative mappings, and operation metadata are not
needed.

### `PredicateWidth`

`PredicateWidth` describes which exact integer predicate path can support a
pair of coordinate spans:

- `Int64`: the faster 64-bit path;
- `Int128`: the wider path, available when the compiler provides `__int128`;
- `Unsupported`: the requested spans are outside the certified range.

Most callers do not need to select a predicate width. The integer triangulation
calls check the input and choose automatically. The public selector is useful
for validating a domain before loading or generating a large point set:

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

### Constrained Delaunay input

```cpp
std::vector<delaunay32::Constraint> constraints = {
    {0, 2},
    {2, 4},
};

const std::vector<delaunay32::Triangle> triangles =
    triangulator.triangulate_constrained_int(points, constraints);
```

This call first constructs the ordinary integer Delaunay triangulation, then
recovers every constraint with topology-preserving edge flips. A requested
segment appears directly as an edge when possible. If one or more existing
input points lie in its interior, it appears as a chain of consecutive mesh
edges instead. All remaining flippable interior edges are locally Delaunay.

Constraints may share endpoints, overlap along the same line, or meet at an
existing input point. Repeated and reversed copies are accepted. Proper
crossings away from an existing input point are rejected; the implementation
does not insert intersection or Steiner points. Both endpoint indices must be
valid and must resolve to different retained coordinates. Duplicate point
indices resolve to their lowest-index representative.

The returned faces cover the full convex hull of the input. Constraint recovery
is currently serial; initial topology construction and final triangle export
still use the requested worker count for sufficiently large inputs. There is
not yet a constrained full-result overload.

### Polygon domains with holes

```cpp
std::vector<std::uint32_t> outer = {0, 1, 2, 3};
std::vector<std::vector<std::uint32_t>> holes = {
    {4, 5, 6, 7},
    {8, 9, 10, 11},
};

const std::vector<delaunay32::Triangle> triangles =
    triangulator.triangulate_polygon_int(points, outer, holes);
```

Every ring contains indices into `points`; its last-to-first edge is implicit.
A repeated first index at the end is optional. Clockwise and counterclockwise
rings are both accepted and normalized internally. Existing input points on a
boundary segment split that segment into a constrained edge chain without
requiring the points to be named in the ring.

The returned triangles cover the outer-ring interior minus every hole. Input
points outside the outer ring or inside a hole are accepted and omitted from
the result. Ring indices resolve through the same lowest-index coordinate
representatives as the other APIs; repeating one retained coordinate within a
ring is invalid.

Rings must be simple, mutually disjoint, and non-touching. Each hole must lie
strictly inside the outer ring; holes may not overlap, touch, or nest. Boundary
validation is currently quadratic in the total ring-edge count. The method
does not create vertices at intersections and does not insert Steiner points.
Boundary recovery and domain classification are serial; large initial
triangulations and final exports still use the requested worker count. There is
not yet a polygon full-result overload.

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

## Floating-point triangulation

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

Request the full result when the mapping precision or grid collisions matter:

```cpp
const auto result = triangulator.triangulate_float_full(points);
const auto& report = result.quantization;

if (report.collapsed_points != 0) {
    // Some input points shared a quantized grid coordinate.
}
```

### Configuring the mapping

Automatic mode is the default. A caller-requested step can intentionally use a
coarser, application-defined grid:

```cpp
delaunay32::QuantizationOptions options;
options.mode = delaunay32::QuantizationMode::GridStep;
options.grid_step = 0.001;
options.max_coordinate_error = 0.0005;

const auto result = triangulator.triangulate_float_full(points, options);
```

`GridStep` places the origin at the current input minima. To keep the same
mapping across independent batches, supply both origin and scale:

```cpp
delaunay32::QuantizationOptions options;
options.mode = delaunay32::QuantizationMode::FixedScale;
options.origin_x = 500000.0;
options.origin_y = 5800000.0;
options.scale = 100.0; // one-centimetre grid when source units are metres
options.collision_policy =
    delaunay32::QuantizationCollisionPolicy::Reject;
```

A requested step or fixed mapping must still produce signed 32-bit coordinates
whose spans fit an exact predicate path. A positive `max_coordinate_error`
rejects a call when its measured error is larger. `Reject` rejects any shared
grid coordinate, including an exact source duplicate; use a full result with
`Allow` to inspect collision counts and representatives instead.

### How automatic quantization works

Finite floats are accepted. The library translates the input by its minimum
bounds and applies one shared scale to both axes:

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

- The input floating-point values are never replaced or modified.
- Every triangle index addresses the caller's original input.
- Aspect ratio is preserved by the uniform x/y scale.
- Exact duplicates and quantization collisions retain their lowest input index.

### What can change

The connectivity is the exact Delaunay triangulation of the quantized integer
points, not necessarily the Delaunay triangulation of the original
floating-point coordinates. Quantization can:

- merge points that are closer than the available grid resolution;
- select a different diagonal for nearly cocircular configurations;
- perturb oblique collinearity;
- change winding when the original geometry is nearly degenerate.

Axis-aligned collinearity is preserved. For most graphics, mapping,
visualization, and general meshing applications these differences are
negligible. Use an adaptive-exact triangulator when the precise Delaunay edge
set of the source coordinates is a requirement.

## Common output behavior

At least three input points and at least three unique triangulation coordinates
are required. Coincident coordinates are collapsed deterministically, retaining
the lowest original input index. Three or more unique but entirely collinear
points are valid input and return an empty triangle vector because no
two-dimensional faces exist.

Do not assume that every input point appears in a returned triangle: duplicates
or floating-point values that collide on the quantization grid can be omitted
in favor of their representative.

For cocircular input, more than one Delaunay triangulation can be valid. Compare
meshes geometrically rather than assuming that another implementation must
choose the same diagonal.

## Invalid input and exceptions

The triangulation functions throw `std::invalid_argument` when:

- fewer than three input points are supplied;
- fewer than three unique integer or quantized coordinates remain;
- a float coordinate is NaN or infinite;
- quantization options are invalid, exceed an error limit, reject a collision,
  or map outside the supported integer domain;
- integer coordinate spans exceed the certified predicate range;
- the input count exceeds the internal index or edge-arena range.

`triangulate_constrained_int()` also throws `std::invalid_argument` when a
constraint endpoint is outside the point array, when both endpoints resolve to
the same retained coordinate, or when constraints cross away from an existing
input point.

`triangulate_polygon_int()` also throws `std::invalid_argument` for invalid or
repeated ring coordinates, self-intersections, touching or crossing rings,
holes outside the outer ring, and overlapping or nested holes.

As with other allocating C++ interfaces, allocation failures can propagate as
standard library exceptions.

## Complete examples

- [`delaunay_float_example.cpp`](../examples/delaunay_float_example.cpp)
  generates float points, requests a quantization report, and writes the mesh
  using the original float coordinates.
- [`delaunay_svg_example.cpp`](../examples/delaunay_svg_example.cpp) accepts
  integer JSON input or generates random integer points and writes an SVG.
- [`delaunay_constrained_example.cpp`](../examples/delaunay_constrained_example.cpp)
  loads a fixed JSON point/constraint set and writes ordinary and constrained
  meshes side by side.
- [`delaunay_polygon_example.cpp`](../examples/delaunay_polygon_example.cpp)
  loads indexed outer and hole rings, triangulates the domain, and renders
  retained and omitted points.

Build all examples with:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target \
  delaunay_float_example \
  delaunay_svg_example \
  delaunay_constrained_example \
  delaunay_polygon_example
```

The example JSON schema stores points as `[x, y]`, constraints as endpoint-index
pairs, and polygon rings as `polygon.outer` plus `polygon.holes`. The parser
belongs only to the examples; the installed library has no JSON dependency.
