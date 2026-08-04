# Changelog

Notable changes to Delaunay32 are documented here. The project follows
[Semantic Versioning](https://semver.org/).

## Unreleased

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
