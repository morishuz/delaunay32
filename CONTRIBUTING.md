# Contributing to Delaunay32

Bug reports, focused pull requests, and reproducible performance findings are
welcome. Please search the existing issues before opening a new one.

## Reporting a bug

Include the smallest point set that reproduces the problem, the API call and
thread count, the compiler and platform, and the observed error or invalid
mesh. For floating-point input, include the `QuantizationOptions` and
`QuantizationReport` when available.

Do not report security vulnerabilities in a public issue. Follow
[`SECURITY.md`](SECURITY.md) instead.

## Building and validating

Clone the benchmark submodule, configure a release build, and run the complete
test suite:

```sh
git submodule update --init --recursive
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

Changes to the triangulation hot path should also run the full benchmark. New
public behavior should include validation coverage, documentation, and an
installed-package consumer check when relevant.

## Pull requests

Keep changes narrowly scoped and avoid adding a runtime dependency to the
library target. Explain the behavior and motivation, list the checks run, and
call out API, performance, portability, or exactness implications. All source
files use the MIT SPDX header and compile as C++17 with the warnings configured
in `CMakeLists.txt`.
