# Usage

The integer triangulation API is declared in
[`delaunay32/delaunay.hpp`](../include/delaunay32/delaunay.hpp). Explicit
floating-point conversion is declared separately in
[`delaunay32/quantization.hpp`](../include/delaunay32/quantization.hpp).
All public types are in the `delaunay32` namespace.

## API at a glance

```cpp
Triangulator();
void set_options(TriangulationOptions options);
void set_points(const std::vector<Point>& points);
void set_constraints(std::vector<Constraint> constraints);
void set_polygons(std::vector<PolygonDomain> polygons);
TriangulationResult triangulate();
```

A problem has one point array and optional constraint and polygon collections.
The normal sequence is:

1. Set options if the defaults are not suitable.
2. Call `set_points()` to begin a problem.
3. Set or replace constraints and polygon domains.
4. Call `triangulate()` once.

`triangulate()` consumes the problem, whether it succeeds or throws.
Another `set_points()` is required before any geometry setter or another
`triangulate()` call can be used.

## Integer points

`Point` stores two signed 32-bit coordinates. `Triangle` stores three
`std::uint32_t` indices into the configured point vector:

```cpp
struct Point {
    std::int32_t x;
    std::int32_t y;
};

struct Triangle {
    std::uint32_t i0;
    std::uint32_t i1;
    std::uint32_t i2;
};
```

Basic triangulation requires only the points:

```cpp
#include <delaunay32/delaunay.hpp>

#include <vector>

const std::vector<delaunay32::Point> points = {
    {0, 0},
    {100, 0},
    {100, 100},
    {0, 100},
    {48, 37},
};

delaunay32::Triangulator triangulator;
triangulator.set_points(points);
const delaunay32::TriangulationResult result =
    triangulator.triangulate();
```

Triangles are counterclockwise and refer to the original point vector.
Coincident points are collapsed deterministically to their lowest original
input index.

Orientation and in-circle decisions are exact for every accepted input. Large
absolute offsets are allowed; the supported range depends on the x and y
spans. Inputs outside the certified predicate range throw
`std::invalid_argument`.

## Problem lifecycle

`set_points()` is the reset boundary. It starts a new problem and clears:

- previous constraints and polygon domains;
- previous output and run state;
- a consumed or failed state from the preceding problem.

`TriangulationOptions`, reusable allocations, and worker threads persist.

`set_constraints()` and `set_polygons()` replace their respective
collections. Passing an empty vector clears that collection:

```cpp
triangulator.set_points(points);
triangulator.set_constraints(first_constraints);
triangulator.set_constraints(replacement_constraints);
triangulator.set_polygons({}); // No polygon clipping for this run.
const auto result = triangulator.triangulate();
```

Calling a geometry setter or `triangulate()` without a ready point set throws
`std::logic_error`. Invalid input values throw `std::invalid_argument`.

Because polygon and constraint indices refer to the configured point vector,
they must be submitted again after every `set_points()` call. This prevents
geometry from one problem being applied accidentally to another point order.

## Options and threading

```cpp
delaunay32::TriangulationOptions options;
options.thread_count = 0;
options.result_detail = delaunay32::ResultDetail::Full;
triangulator.set_options(options);
```

`thread_count = 1` selects the serial path and is the default. Zero selects
the hardware thread count. A larger value requests that many threads, subject
to the implementation limit. Small inputs still use the serial path.

One `Triangulator` is not safe for concurrent calls. Use independent
instances for independent concurrent problems.

## Constraints

A `Constraint` is an undirected segment whose endpoints index the configured
point vector:

```cpp
const std::vector<delaunay32::Constraint> constraints = {
    {0, 2},
    {2, 4},
};

triangulator.set_points(points);
triangulator.set_constraints(constraints);
const auto result = triangulator.triangulate();
```

The result preserves each requested segment as an edge or a chain of edges
through existing collinear points. Constraints may share endpoints but cannot
cross away from an existing point. No intersection or Steiner points are
inserted.

When constraints and polygon domains are configured together, all segments
and polygon boundaries are recovered in the same topology. Constraint portions
outside the retained polygon union are clipped from the exported triangles.
Intersections with polygon boundaries follow the same noncrossing rules.

![Ordinary Delaunay triangulation beside the constrained result, with requested segments highlighted in red](../images/delaunay32_constrained.svg)

## Polygon domains

`PolygonDomain` contains one outer index ring and zero or more hole rings:

```cpp
delaunay32::PolygonDomain domain;
domain.outer_ring = {0, 1, 2, 3};
domain.holes = {
    {4, 5, 6, 7},
};

triangulator.set_points(points);
triangulator.set_polygons({domain});
const auto result = triangulator.triangulate();
```

The closing edge of each ring is implicit, and either winding direction is
accepted. A result contains triangles inside its outer ring and outside every
hole. Other configured points can remain unused.

Multiple domains are processed as one disjoint union:

```cpp
triangulator.set_points(points);
triangulator.set_polygons({first_domain, second_domain, third_domain});
const auto result = triangulator.triangulate();
```

The unconstrained topology is built once, every boundary is recovered
together, and the union is legalized, clipped, and exported once. Returned
triangles do not carry domain labels.

Each ring must be simple. Holes must be strictly inside their outer ring and
cannot overlap, touch, or nest. Outer domains must also be disjoint: touching,
overlapping, and nested outer domains are rejected.

![Constrained Delaunay triangulation of a polygon with three holes](../images/delaunay32_polygon.svg)

## Result detail and reports

Triangle detail is the default:

```cpp
triangulator.set_points(points);
const auto result = triangulator.triangulate();
```

`result.triangles` is populated while `halfedges`, `hull`, and
`representatives` remain empty. The operation report is populated for every
successful run.

Select full detail through persistent options:

```cpp
delaunay32::TriangulationOptions options;
options.result_detail = delaunay32::ResultDetail::Full;
triangulator.set_options(options);

triangulator.set_points(points);
const delaunay32::TriangulationResult result =
    triangulator.triangulate();
```

| Field | Contents |
|:--|:--|
| `triangles` | Counterclockwise faces using original input indices |
| `halfedges` | Opposite flattened edge, or `-1` at a convex-hull or clipped-domain boundary |
| `hull` | Counterclockwise convex hull of the complete unique input set |
| `representatives` | One entry per input, mapping coincident sites to the lowest retained index |
| `report.predicate_width` | Exact predicate path selected for the coordinate spans |
| `report.actual_thread_count` | Threads actually used for the operation |
| `report.input_points` | Number of configured point entries |
| `report.unique_points` | Number of distinct integer coordinates |
| `report.collapsed_points` | Number of coincident entries collapsed |

For triangle `t`, flattened edges `3*t`, `3*t+1`, and `3*t+2`
correspond to `i0->i1`, `i1->i2`, and `i2->i0`.

For polygon output, `hull` remains the hull of the complete unique input set,
including points outside the retained domains. Clipping affects triangles and
halfedge boundaries, not this input-set report.

## Explicit floating-point quantization

The triangulator accepts only integer `Point` values. Convert finite
`FloatPoint` input explicitly before configuring it:

```cpp
#include <delaunay32/quantization.hpp>

const std::vector<delaunay32::FloatPoint> source = {
    {0.125, 0.25},
    {5.5, 0.1},
    {6.0, 4.5},
    {-1.0, 5.0},
};

const delaunay32::QuantizationResult converted =
    delaunay32::quantize(source);

triangulator.set_points(converted.points);
const auto mesh = triangulator.triangulate();
```

`quantize()` preserves the source vector's length and index order. It reports
the selected mapping, maximum measured coordinate error, unique coordinates,
and collisions. It never invokes `Triangulator`.

`FloatPoint` stores double-precision coordinates. Every signed 32-bit integer
coordinate can therefore be converted to `FloatPoint` exactly.

Automatic quantization is the default. Use explicit settings when an
application needs a particular grid or a stable mapping across batches:

```cpp
delaunay32::QuantizationOptions options;
options.mode = delaunay32::QuantizationMode::FixedScale;
options.origin_x = 0.0;
options.origin_y = 0.0;
options.scale = 1000.0;

const auto converted = delaunay32::quantize(source, options);
```

| Mode | Mapping |
|:--|:--|
| `Automatic` | Finest supported uniform grid for the supplied input |
| `GridStep` | Caller-provided `grid_step`, using input minima as the origin |
| `FixedScale` | Caller-provided `origin_x`, `origin_y`, and `scale` |

Automatic mode rejects an input whose coordinate span overflows double
precision, even when every individual coordinate is finite.

`max_coordinate_error` rejects a mapping that exceeds an application limit.
Set `collision_policy` to `QuantizationCollisionPolicy::Reject` to reject
both exact duplicates and distinct source points that land on the same integer
coordinate. With `Allow`, collisions are retained in the converted vector and
subsequently handled by the triangulator's representative mapping.

## Copying and transferring geometry

`set_constraints()` and `set_polygons()` accept vectors by value. Passing
an lvalue copies it, while `std::move` explicitly transfers ownership:

```cpp
triangulator.set_polygons(polygons);            // Keep the caller's copy.
triangulator.set_polygons(std::move(polygons)); // Transfer its allocation.
```

`set_points()` takes a const reference because public points are converted
into the triangulator's internal working representation. Callers also normally
retain the point vector because returned triangles contain indices rather than
coordinate copies.

`Triangulator` itself is movable but not copyable. Moving a configured
instance transfers its pending problem, options, retained allocations, and
worker team.

## Input behavior

At least three point entries and three unique integer coordinates are required.
Three or more collinear coordinates are valid and return no triangles.
Cocircular inputs can have more than one valid Delaunay edge set.

Invalid coordinates, options, constraints, rings, counts, or certified spans
throw `std::invalid_argument`. Allocation failures can propagate as standard
library exceptions. After any failed run, call `set_points()` before retrying.

## Optional extras

Link `delaunay32::extras` for:

- Delaunay32 geometry JSON input and output;
- uniform and polygon-interior point sampling;
- domain queries;
- composable SVG drawing and one-call mesh visualizations.

### SVG

The `Svg` builder records geometry in world coordinates and auto-fits it when
rendered. Points, line widths, and font sizes remain in canvas units, so they
stay legible regardless of the coordinate range:

```cpp
#include <delaunay32/extras/svg.hpp>

delaunay32::extras::Svg svg(1200.0, 800.0);
svg.set_background("#f7f7f5");
svg.set_auto_fit({30.0, 70.0, 30.0, 30.0});

delaunay32::extras::SvgTriangleColorStyle mesh_style;
mesh_style.stroke = "#344044";
mesh_style.stroke_width = 0.75;
svg.draw_colored_triangles(points, result.triangles, mesh_style);

delaunay32::extras::SvgPointStyle point_style;
point_style.fill = "#2479a6";
point_style.radius = 2.0;
svg.draw_points(points, point_style);
svg.draw_text("triangulation", 30.0, 38.0);
svg.render_to_svg("mesh.svg");
```

`draw_point()` and `draw_line()` provide individual primitives.
`draw_polygon(points, ring)` closes one index ring, while
`draw_polygon(points, domain)` emits an even-odd compound path for an outer
ring and all its holes. Integer `Point` and floating-point `FloatPoint` vectors
are both accepted. Draw order is preserved.

`draw_triangles()` applies one `SvgShapeStyle` to the complete triangle batch.
`draw_colored_triangles()` uses a configurable palette and guarantees that
triangles sharing an indexed edge receive different palette entries. Triangles
that meet only at one point are not considered adjacent. A manifold triangle
mesh needs at most four distinct palette colors; the default style provides
the eight-color diagnostic palette used by the one-call exporters. Smaller
palettes are accepted when they can color the supplied mesh; a failure reports
the supplied count and recommends four colors for a guaranteed result.

Autoscaling is the default and uses Cartesian orientation, with larger y
values appearing higher on the canvas. `set_auto_fit()` receives margins in
`{left, top, right, bottom}` order and in the same canvas units as the SVG
width and height. For a `1200` by `440` canvas, margins of
`{28, 68, 28, 28}` leave a `1144` by `344` fitting rectangle. Geometry is
uniformly scaled and centered inside that rectangle at render time. Text
positions, point radii, and stroke widths remain fixed canvas sizes and do not
affect the fitted bounds.

`set_transform(x_scale, y_scale, x_offset, y_offset)` disables autoscaling and
applies that mapping directly to SVG coordinates. Use a negative y scale when
an explicit transform should retain Cartesian orientation.

Fill, stroke, and background colors accept CSS color strings. Transparency can
be encoded directly, for example as `#2479a680` or
`rgb(36 121 166 / 50%)`; these alpha forms are defined by
[CSS Color 4](https://www.w3.org/TR/css-color-4/). The strings are preserved in
the SVG output without introducing separate opacity settings.

`to_svg()` returns the complete document as a string. `render_to_svg()` writes
it to a file. The older `write_mesh_svg()` and `write_polygon_mesh_svg()`
functions remain as opinionated diagnostic shortcuts built from the same
primitives.

### JSON

The geometry JSON helpers use the same `Geometry` model for files and
in-memory serialization:

```cpp
#include <delaunay32/extras/json.hpp>

delaunay32::extras::Geometry geometry;
geometry.points = points;
geometry.constraints = constraints;
geometry.polygons = domains;

const std::string json =
    delaunay32::extras::geometry_to_json(geometry);
delaunay32::extras::write_geometry_json("geometry.json", geometry);

const delaunay32::extras::Geometry restored =
    delaunay32::extras::read_geometry_json("geometry.json");
```

`points` is required. `constraints`, one `polygon`, or a `polygons` array are
written only when present; `polygon` and `polygons` are mutually exclusive.
Coordinates must be signed 32-bit integers and topology entries are unsigned
point indices. The reader rejects duplicate or unknown fields and reports the
source file, line, and column for malformed input. `geometry_to_json()` and
`write_geometry_json()` perform the same geometry validation and produce
identical formatted JSON.

### Point generation

`PointSampler` generates double-precision points inside either rectangular
bounds or indexed polygon domains:

```cpp
#include <delaunay32/extras/sampling.hpp>

delaunay32::extras::PointSampler sampler;
sampler.set_bounds({0.0, 100.0, 0.0, 60.0});

delaunay32::extras::UniformSamplingOptions uniform;
uniform.point_count = 500;
uniform.seed = 42;
const std::vector<delaunay32::FloatPoint> random_points =
    sampler.generate_uniform(uniform);

delaunay32::extras::BlueNoiseSamplingOptions blue_noise;
blue_noise.point_count = 500;
blue_noise.candidates_per_point = 16;
blue_noise.seed = 42;
const std::vector<delaunay32::FloatPoint> spaced_points =
    sampler.generate_blue_noise(blue_noise);

delaunay32::extras::JitteredGridSamplingOptions fast_spaced;
fast_spaced.point_count = 500;
fast_spaced.jitter = 0.75;
fast_spaced.seed = 42;
const std::vector<delaunay32::FloatPoint> demo_points =
    sampler.generate_jittered_grid(fast_spaced);
```

`SamplingBounds` uses `{min_x, max_x, min_y, max_y}` order. All algorithms
operate on either configured region. Calling `set_bounds()` or
`set_polygon_interiors()` replaces the previous region, and generation before
setting a region throws `std::logic_error`. Each call owns an independent
random engine, so the same region, options, and seed produce the same result
regardless of earlier calls.

Polygon regions accept either native integer `Point` arrays or floating-point
`FloatPoint` arrays while reusing the same `PolygonDomain` indices:

```cpp
delaunay32::extras::PointSampler sampler;
sampler.set_polygon_interiors(boundary_points, domains);

delaunay32::extras::BlueNoiseSamplingOptions options;
options.point_count = 600;
const std::vector<delaunay32::FloatPoint> interiors =
    sampler.generate_blue_noise(options);
```

Polygon samples are strictly inside an outer ring and outside its holes.
Multiple domains form one sampling region. Uniform sampling uses rejection
sampling over that region. Best-candidate blue noise maximizes distance from
domain boundaries and previously accepted samples. It uses an exact spatial
nearest-neighbor index, but still evaluates `candidates_per_point` polygon
candidates for every result.

`generate_jittered_grid()` is the fast visual-demo path. It builds a triangular
lattice with a seed-derived random rotation and phase offset, then perturbs
each site in a uniformly random direction inside a disk. Candidates are
clipped to the polygon region, and the sampler returns exactly `point_count`
points. A `jitter` of zero preserves the lattice sites; one permits movement
by half the lattice spacing. The default `0.75` is a compromise between strong
spacing and visible lattice structure. Lattice density is fitted to the
observed clipped count. Small surpluses are removed from locally crowded,
mutually separated sites; small shortages are filled with best-candidate gap
samples. This avoids randomly punching visible holes into a low-jitter
lattice. Use best-candidate blue noise when boundary clearance matters more
than generation time; use the jittered grid when sampling should be cheaper
than triangulation. The attempt limits bound work for very small or sparse
polygon regions.

Uniform sampling does not inject fixed points by default. Set
`include_bounds_corners` to `true` to place the distinct rectangle corners at
the beginning of a bounds result; these corners count toward `point_count`.
That option is invalid for polygon interiors because polygon boundary vertices
are not interior samples.

When generated points will be triangulated with floating-point boundary
geometry, append everything before quantizing. One shared conversion preserves
the boundary indices and applies the same coordinate mapping to every point:

```cpp
std::vector<delaunay32::FloatPoint> all_points = boundary_points;
all_points.insert(all_points.end(), interiors.begin(), interiors.end());

const delaunay32::QuantizationResult converted =
    delaunay32::quantize(all_points);
triangulator.set_points(converted.points);
triangulator.set_polygons(domains);
```

These utilities are separate from the dependency-free core target. Their
headers are under
[`include/delaunay32/extras`](../include/delaunay32/extras). The focused
[logo program](../examples/delaunay32_logo_polygon_example.cpp) demonstrates
JSON input, polygon-interior sampling, one batched multi-domain triangulation,
and SVG output.
