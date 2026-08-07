# Changelog

Notable changes to Delaunay32 are documented here. The project follows
[Semantic Versioning](https://semver.org/).

## 0.6.1 - 2026-08-08

### Added

- A deterministic, randomly rotated and phase-shifted triangular-lattice
  sampler with isotropic jitter for fast blue-noise-like point sets in
  rectangular bounds and indexed polygon interiors.

### Changed

- Floating-point polygon containment and sampling now use binary64 arithmetic
  throughout, avoiding software-emulated IEEE-128 operations in WebAssembly.
- Best-candidate sampling now uses per-domain bounds and an exact spatial
  nearest-neighbor index instead of scanning every accepted sample.
- Jittered-grid sampling fits lattice density to the clipped region and uses
  crowding-aware trimming or gap filling for small count corrections, avoiding
  the visible holes caused by random removal from a fixed safety surplus.

## 0.6.0 - 2026-08-06

### Added

- A stateful integer-only `Triangulator` configured with `set_options()`,
  `set_points()`, `set_constraints()`, and `set_polygons()`, followed by one
  owning `TriangulationResult` from `triangulate()`.
- `TriangulationOptions`, `ResultDetail`, and `TriangulationReport`. Triangle
  output remains the default; full detail additionally returns halfedges, the
  complete input hull, and duplicate representatives.
- A standalone `delaunay32/quantization.hpp` API. `quantize()` explicitly maps
  `FloatPoint` input to an index-preserving `std::vector<Point>` and returns a
  separate error and collision report.
- Batched disjoint polygon domains and support for standalone constraints in
  the same polygon-clipped run.
- A composable `extras::Svg` document with auto-fit or explicit transforms,
  styled point and line primitives, indexed polygons with holes, triangle
  batches, deterministic edge-adjacent triangle coloring, text, in-memory
  rendering, and file output. The existing one-call mesh exporters now use
  the same drawing API.
- `geometry_to_json()` for validated in-memory serialization using the same
  formatted geometry schema as `write_geometry_json()`.
- A stateful float-only `extras::PointSampler` with deterministic uniform and
  boundary-aware blue-noise generation over rectangular bounds or indexed
  polygon interiors.
- Five focused examples, built through `DELAUNAY32_BUILD_EXAMPLES`, covering
  basic meshes, quantization, constraints, polygons with holes, and a batched
  multi-domain logo.
- Lifecycle, result-detail, report, quantization, moved-instance, batched
  polygon, mixed-constraint, installed-package, stress, fallback, and example
  smoke-test coverage for the new contract.

### Changed

- The core triangulator now accepts signed 32-bit integer points exclusively.
  Floating-point conversion is an explicit preprocessing operation and never
  a triangulator overload.
- A run builds the unconstrained topology once, recovers standalone
  constraints and all polygon boundaries together, legalizes once, clips the
  disjoint domain union once, and exports one combined triangle vector.
- `PolygonDomain` is now a core type. Touching, overlapping, and nested outer
  domains are rejected; existing per-domain hole validation remains in force.
- Each configured problem is consumed by `triangulate()`, including failed
  runs. `set_points()` starts the next problem and clears geometry and run
  state while preserving options, allocations, and worker threads.
- Full polygon results use `-1` halfedges across every clipped boundary while
  preserving `hull` as the convex hull of the complete unique input set.
- `FloatPoint` coordinates now use double precision, allowing every signed
  32-bit integer coordinate to be represented exactly before quantization.
- Uniform sampling no longer injects rectangle corners by default; callers can
  request them explicitly with `include_bounds_corners`.
- The benchmark now measures Delaunay32 one-thread and automatic modes across
  its three point distributions without an embedded competitor baseline.
- Randomized validation uses exact mesh invariants and serial/parallel result
  agreement without an external triangulation implementation.
- Package metadata and installed-package consumers now target version 0.6.0.

### Removed

- The overload-based `triangulate_int()`, `triangulate_float()`,
  `triangulate_int_full()`, `triangulate_float_full()`,
  `triangulate_constrained_int()`, and `triangulate_polygon_int()` entry
  points. There are no compatibility wrappers for this pre-1.0 breaking
  release.
- The former integer and floating-point free sampling functions and their
  option structs. Sampling is now configured through `PointSampler` and always
  returns `FloatPoint`.
- The Delaunator adapter, Git submodule, build target, and benchmark columns.

### Fixed

- Automatic quantization now rejects coordinate sets whose double-precision
  span overflows instead of silently collapsing the complete input.
- Embedded CMake builds no longer force the parent project to use the Release
  build type.

## 0.5.1 - 2026-08-04

### Added

- An optional, separately linked `delaunay32::extras` companion library with
  unique integer point generation, float point generation, boundary-aware
  best-candidate polygon sampling, exact integer domain queries, Delaunay32
  geometry JSON input/output, and diagnostic SVG export.
- Installed-package and round-trip validation for the extras target, including
  sampling, geometry JSON, and ordinary and polygon SVG output.

### Changed

- Constrained triangulation now legalizes only the neighborhoods changed by
  constraint-recovery flips instead of scanning the complete mesh.
- Constraint endpoint indexing now initializes only requested endpoints and
  shares the retained worker team in automatic mode.
- One-million-point constrained workloads are about 30% faster with one
  thread and 50–58% faster in automatic mode on the reference Apple M1 system,
  without a measurable ordinary-Delaunay regression.
- The examples now delegate reusable sampling, geometry, JSON, and generic SVG
  work to `delaunay32::extras`; only CLI handling and bespoke comparison/logo
  presentation remain example-private.

## 0.5.0 - 2026-08-04

### Added

- `Constraint` and `triangulate_constrained_int()`, providing exact
  constrained Delaunay triangulation of integer point sets with noncrossing
  segments, collinear segment chains, duplicate representatives, and
  constraints meeting at existing sites.
- `triangulate_polygon_int()`, providing constrained Delaunay triangulation of
  simple integer polygon domains with multiple holes, automatic winding
  normalization, boundary chains through existing sites, and domain clipping.
- A fixed constrained-triangulation comparison example and an extensible,
  dependency-free example JSON schema for points, constraints, and polygon
  rings. Polygon-with-holes examples exercise domain filtering with both a
  fixed domain and a multi-domain Delaunay32 logo.

### Changed

- Existing CSV example fixtures now use JSON.

### Fixed

- The multi-domain logo fixture stays within the certified 64-bit predicate
  range, allowing the polygon example to run on MSVC and other platforms
  without `__int128`.

## 0.4.0 - 2026-08-03

### Added

- Configurable floating-point quantization with automatic, requested-grid-step,
  and fixed-origin/scale modes, plus error limits and collision rejection.

### Fixed

- Parallel topology and export now abort their phase barriers when a worker
  throws, preventing an exception from stranding another worker while retaining
  the existing single-run fast path.

## 0.3.0 - 2026-08-02

### Added

- `triangulate_int_full()` and `triangulate_float_full()`, returning triangles,
  halfedge adjacency, the convex hull, duplicate representative mappings,
  predicate width, actual thread count, and float quantization metadata.
- Direct full-result extraction from the compact internal topology, avoiding
  hash-based adjacency reconstruction.
- Deterministic recovery from parallel edge-arena exhaustion by rebuilding on
  the growable serial path.
- Stress, concurrency, forced-fallback, installed-package, ThreadSanitizer, and
  forced-64-bit-predicate validation.
- Windows MSVC builds and installed-package checks in GitHub Actions.

### Changed

- Full-result extraction is parallelized and reuses topology information
  produced during triangulation.
- Documentation now covers Windows builds, complete-result semantics, arena
  fallback behavior, and the streamlined input API.

### Removed

- Raw pointer overloads accepting separate x/y coordinate arrays. Callers
  should pass `std::vector<Point>` or `std::vector<FloatPoint>`.
- `triangulate_float(points, QuantizationReport&)`. Quantization metadata is
  available as `triangulate_float_full(points).quantization`.

These removals are intentional pre-1.0 breaking API changes and are the reason
for the minor-version increment from 0.2 to 0.3.

## 0.2.0 - 2026-07-31

- Initial public release with exact signed-integer predicates, direct float
  quantization, deterministic duplicate handling, parallel triangulation,
  CMake package installation, examples, and validated performance benchmarks.
