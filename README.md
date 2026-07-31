# Delaunay32

**Fast, exact, parallel 2D Delaunay triangulation for signed 32-bit integer coordinates, with convenient float quantization.**

`Delaunay32` is a C++17 library for triangulating large sets of discrete 2D
points: pixels, raster samples, voxel projections, fixed-point geometry, and
other quantized spatial data. Finite `float` input can also be mapped
automatically onto the certified integer domain while output indices continue
to reference the original coordinates.

It combines exact integer predicates with a Morton-ordered divide-and-conquer
algorithm, compact two-dart topology, and optional multithreading. The result is
a triangulator that is deterministic, robust, and particularly fast on large
point sets.

For large point sets, Delaunay32 is over 10× faster than [delaunator-cpp](https://github.com/delfrrr/delaunator-cpp) and around 3× faster than [Fade2D](https://www.geom.at/products/fade2d/).

![Delaunay32 mesh of 1,000 integer-coordinate points filling a square with 1,994 triangles](images/delaunay32_mesh.svg)

## Highlights

- Exact orientation and in-circle predicates for certified coordinate ranges
- Signed 32-bit integer input, including negative coordinates and large offsets
- Serial and shared-memory parallel execution
- Deterministic handling of duplicate points
- Triangle indices referencing the original input, counterclockwise on the
  triangulation grid
- Efficient struct-of-arrays API for existing integer buffers
- Uniform float quantization with an optional precision and collision report
- MIT licensed and dependency-free for normal library use

## When to use it

`Delaunay32` is intended for data that is already discrete or can tolerate a
high-resolution uniform quantization. Typical examples include image-space
geometry, raster and height-field samples, projected voxel data, fixed-point
maps, graphics, and projected spatial datasets.

For floating-point data, `triangulate_float()` keeps the source coordinates
untouched and returns indices into them, but it chooses edges using quantized
integer coordinates. This is suitable when preserving the original vertices is
important but small topology changes near degeneracies are acceptable. Use a
native floating-point or adaptive-exact triangulator when the Delaunay topology
of the source floats themselves is the required result.

## Performance

The benchmark compares `Delaunay32` with
[Delaunator](https://github.com/mapbox/delaunator), a widely used,
high-performance open-source Delaunay implementation originally developed by
Mapbox. It uses
[delaunator-cpp](https://github.com/delfrrr/delaunator-cpp), whose maintainers
describe it as probably the fastest open-source C++ implementation.

Delaunator-cpp is included only as an optional benchmark submodule. `Delaunay32`
does not wrap it or depend on it as part of the triangulation algorithm.

### Headline results

For one million input points:

| distribution | Delaunay32 1T | Delaunay32 auto | Delaunator | serial speedup | auto speedup |
|:--|--:|--:|--:|--:|--:|
| uniform | 147.909 ms | 51.741 ms | 551.592 ms | 3.73× | 10.66× |
| clustered | 143.195 ms | 47.979 ms | 542.438 ms | 3.79× | 11.31× |
| diagonal | 144.347 ms | 49.935 ms | 545.491 ms | 3.78× | 10.92× |

A same-process benchmark (not included in this repo) at one million points measured Delaunay32 to be:

3.7–4.3× as fast as single-threaded [Fade2D](https://www.geom.at/products/fade2d/).    
2.6–3.5× as fast as automatically threaded [Fade2D](https://www.geom.at/products/fade2d/).   

As always, benchmark results are machine- and workload-specific. Run the
included benchmark on your target hardware before making deployment decisions.

### Full benchmark

These are median end-to-end times from a Release build on an Apple M1
(8 cores, 16 GB) using Apple Clang 21. The deterministic inputs contain unique
points in `[0, 100000)²`.

Input generation and validation are excluded. Input ingestion by `Delaunay32`
and output triangle materialization by both implementations are included. The
default fresh-instance lifecycle also includes `Triangulator` construction,
working-storage teardown, and worker-thread lifetime for every sample.

Delaunator's `std::vector<double>` is prepared once outside the timed region,
giving it the favorable assumption that its input is already available in its
native representation. In the ratio columns, lower is better.

| distribution | points | Delaunay32 1T ms | Delaunay32 auto ms | threads | Delaunator ms | 1T / del | auto / del |
|:--|--:|--:|--:|--:|--:|--:|--:|
| uniform | 1,000 | 0.147 | 0.147 | 1 | 0.295 | 0.498x | 0.498x |
| uniform | 10,000 | 1.369 | 1.363 | 1 | 3.052 | 0.449x | 0.447x |
| uniform | 100,000 | 14.317 | 5.658 | 8 | 39.837 | 0.359x | 0.142x |
| uniform | 1,000,000 | 147.909 | 51.741 | 8 | 551.592 | 0.268x | 0.094x |
| clustered | 1,000 | 0.093 | 0.093 | 1 | 0.221 | 0.419x | 0.420x |
| clustered | 10,000 | 1.343 | 1.355 | 1 | 3.260 | 0.412x | 0.416x |
| clustered | 100,000 | 14.051 | 5.354 | 8 | 39.629 | 0.355x | 0.135x |
| clustered | 1,000,000 | 143.195 | 47.979 | 8 | 542.438 | 0.264x | 0.088x |
| diagonal | 1,000 | 0.097 | 0.095 | 1 | 0.242 | 0.400x | 0.390x |
| diagonal | 10,000 | 1.376 | 1.375 | 1 | 3.159 | 0.436x | 0.435x |
| diagonal | 100,000 | 14.318 | 5.838 | 8 | 39.682 | 0.361x | 0.147x |
| diagonal | 1,000,000 | 144.347 | 49.935 | 8 | 545.491 | 0.265x | 0.092x |

## Quick start

Clone the repository and initialize the optional benchmark dependency:

```sh
git submodule update --init --recursive
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

For a library-only build, Delaunator is not required:

```sh
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DDELAUNAY32_BUILD_BENCHMARKS=OFF \
  -DDELAUNAY32_BUILD_TESTS=OFF \
  -DDELAUNAY32_BUILD_EXAMPLES=OFF
cmake --build build -j
```

Install with:

```sh
cmake --install build
```

The exported CMake target is `delaunay32::delaunay32`.

## Usage

```cpp
#include <delaunay32/delaunay.hpp>

#include <vector>

int main() {
    std::vector<delaunay32::Point> points = {
        {0, 0},
        {100, 0},
        {100, 100},
        {0, 100},
        {48, 37},
    };

    // 1 selects the serial path; 0 selects the hardware thread count.
    delaunay32::Triangulator triangulator(0);
    const std::vector<delaunay32::Triangle> triangles =
        triangulator.triangulate_int(points);

    for (const auto& triangle : triangles) {
        // i0, i1, and i2 index the original point vector in CCW order.
    }
}
```

For existing struct-of-arrays storage, `triangulate_int(xs, ys, count)` avoids
constructing a temporary `Point` vector. Version 0.2 uses the explicit
`triangulate_int` name for both integer layouts; the former
`triangulate(points)` spelling is not retained as an ambiguous compatibility
alias.

Finite float input uses a separate point type and an explicit method so the
quantized contract is visible at the call site:

```cpp
std::vector<delaunay32::FloatPoint> points = {
    {0.125F, 0.25F},
    {5.5F, 0.1F},
    {6.0F, 4.5F},
    {-1.0F, 5.0F},
};

delaunay32::QuantizationReport report;
const std::vector<delaunay32::Triangle> triangles =
    triangulator.triangulate_float(points, report);

// Triangle indices address the unchanged FloatPoint vector. The report shows
// the grid step, measured coordinate error, and any quantized point collisions.
```

The equivalent struct-of-arrays overload is
`triangulate_float(xs, ys, count, report)`. Both input layouts also have an
overload without a report argument.

A `Triangulator` can be reused across calls to retain working storage and worker
threads. A single instance must not be called concurrently; separate instances
are independent.

## Exact integer geometry

Integer coordinates give the library a clear numerical contract:

- orientation and in-circle decisions can be exact;
- no epsilon or floating-point tie policy is required;
- storage and sorting remain compact;
- duplicate coordinates have an unambiguous meaning.

Coordinates may be negative and may include a large fixed offset. Predicate
selection depends on the coordinate spans

```text
sx = max_x - min_x
sy = max_y - min_y
```

rather than on the absolute coordinate values.

The library inspects the input spans and automatically chooses the narrowest
certified predicate path:

- `Int64`: fast 64-bit predicates and 32-bit lifted coordinates
- `Int128`: 128-bit predicates and 64-bit lifts where `__int128` is available
- `Unsupported`: throws `std::invalid_argument` before topology construction

For equal spans, the largest conservatively certified values are 29,609 for
`Int64` and 1,940,470,527 for `Int128`. Thin rectangular domains may be much
wider. The exact runtime test uses `L = sx² + sy²` and requires:

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

Coincident points are collapsed automatically. Every output index refers to the
lowest original input index at that coordinate. Entirely collinear input has no
two-dimensional faces and returns an empty triangle vector.

## Quantized float geometry

`triangulate_float()` accepts finite `float` coordinates across their full
range. It translates the input by its minimum bounds and applies one shared
scale to both axes:

```text
scale = kMaxCoordinateSpan / max(max_x - min_x, max_y - min_y)
qx = round((x - min_x) * scale)
qy = round((y - min_y) * scale)
```

Using one scale preserves aspect ratio; the shorter axis is not stretched to
fill the grid. The resulting integers are triangulated with the same exact
predicates as direct integer input. On platforms without `__int128`, the target
span automatically falls back to the certified 64-bit limit.

The returned `Triangle` values contain only indices, so callers recover the
original float vertices without coordinate loss. The topology is nevertheless
the topology of `(qx, qy)`, not necessarily of `(x, y)`. Quantization can merge
nearby points, change an edge choice, change winding in nearly degenerate
geometry, or perturb oblique collinearity. Axis-aligned collinearity is
preserved. If fewer than three unique grid points remain, triangulation throws
`std::invalid_argument`.

`QuantizationReport` describes the actual mapping through `origin_x`,
`origin_y`, `scale`, and `grid_step`; records the largest source-space error in
either coordinate as `max_coordinate_error`; and reports `unique_points` and
`collapsed_points`. A collapsed count includes both exact duplicate inputs and
distinct float points that map to the same grid coordinate. NaN and infinity
are rejected.

## How it works

At a high level, Delaunay32:

1. translates coordinates by the input minima and certifies predicate widths;
2. generates Morton keys and radix-sorts the sites;
3. constructs small divide-and-conquer leaves using exact orientation and
   in-circle predicates;
4. merges neighboring triangulations through compact primal edge rings;
5. marks the outer face and materializes indices that are counterclockwise in
   the triangulation coordinates.

For large inputs, radix sorting, independent subtrees, merge levels, and
triangle export share a retained worker team. Small inputs stay serial to avoid
synchronization overhead.

The architecture belongs to the established divide-and-conquer Delaunay family,
particularly:

- L. Guibas and J. Stolfi, *Primitives for the Manipulation of General
  Subdivisions and the Computation of Voronoi Diagrams* (1985)
- R. A. Dwyer, *A Faster Divide-and-Conquer Algorithm for Constructing
  Delaunay Triangulations* (1987)

## Running the benchmark

```sh
./build/delaunay_benchmark
./build/delaunay_benchmark --quick
./build/delaunay_benchmark --reuse
./build/delaunay_benchmark --format markdown
./build/delaunay_benchmark --format csv
./build/delaunay_benchmark --sizes 1000,10000,1000000
```

Each case is checked against Delaunator before timing. Validation first attempts
an exact unoriented triangle-set match. Where cocircular points permit a
different valid diagonal, it instead requires equal triangle counts and checks
manifold edge incidence and exact local Delaunay legality.

The benchmark rotates implementation order and reports medians. Dataset
generation, seeds, validation, and timing boundaries are defined in
`benchmarks/benchmark.cpp` and `benchmarks/support.hpp`.

Full runs use 11 samples through 10,000 points, 7 samples at 100,000 points, and
5 samples above 100,000 points. `--quick` uses 3 samples.

The default `--fresh` mode constructs and destroys a `Triangulator` for each
measured sample, including working-storage and thread-pool lifetime. `--reuse`
models repeated triangulations through one retained instance.

## SVG examples

### Random float points

`delaunay_float_example` is a compact end-to-end example of the quantized float
API. It generates 5,000 deterministic random `FloatPoint` values, triangulates
them with a `QuantizationReport`, prints the mapping and precision information,
and writes an SVG using the unchanged source coordinates addressed by the
returned triangle indices.

```sh
cmake -S . -B build-debug -DCMAKE_BUILD_TYPE=Debug
cmake --build build-debug --target delaunay_float_example --parallel
./build-debug/delaunay_float_example float-mesh.svg
```

On macOS, inspect the result with:

```sh
open float-mesh.svg
```

The example keeps its point count, seed, and rectangular float bounds as named
constants near the top of `examples/delaunay_float_example.cpp`, making them
easy to change while keeping the API flow uncluttered.

### Integer points and CSV input

The dependency-free example accepts arbitrary signed integer points from a
two-column CSV file:

```csv
x,y
0,0
500,0
1000,0
500,500
```

It reads the points exactly as supplied and scales the SVG to the input bounds
while preserving its aspect ratio.

```sh
./build/delaunay_svg_example \
  --input examples/data/delaunay32.csv \
  --output delaunay32.svg
```

With no CSV input, the original random mode remains available:

```sh
./build/delaunay_svg_example 100000 mesh.svg 42
```

Random mode generates unique interior points and inserts the four square
corners. Its point count includes those corners. In CSV mode, the file must
contain every desired point, including any boundary and corner points; the
example does not add or remove points.

## Development

Performance changes should pass the validation suite, improve the full
benchmark geometric mean, and avoid serious regressions in any large case.
Failed experimental kernels should remain outside the main branch so the
production path stays readable.

## Scope

Included:

- signed 32-bit integer coordinates
- finite float coordinates through automatic uniform quantization
- exact certified predicates
- deterministic duplicate handling
- serial and shared-memory parallel construction
- triangle indices that are counterclockwise on the triangulation grid
- up to `2^31 - 1` input entries, subject to available memory

Not included:

- exact predicates on the unquantized floating-point coordinates
- constrained edges, holes, or polygon clipping
- dynamic insertion or deletion
- Voronoi output

## License

The library and repository-owned utilities are MIT licensed. The optional
Delaunator benchmark submodule has its own permissive notices; see
[THIRD_PARTY.md](THIRD_PARTY.md).
