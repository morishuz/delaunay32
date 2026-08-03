# Release process

Delaunay32 follows [Semantic Versioning](https://semver.org/). Before version
1.0, a minor version may contain breaking API changes; patch versions should
remain backward compatible within that minor release.

## Version source of truth

The project version is declared once in the top-level `CMakeLists.txt`:

```cmake
project(delaunay32 VERSION 0.4.0 LANGUAGES CXX)
```

The generated CMake package configuration uses this value. Do not maintain a
second version file unless another packaging system eventually requires it.

Published releases use annotated Git tags named `vMAJOR.MINOR.PATCH`, for
example `v0.4.0`. A release tag identifies an immutable release commit; create
it only after the release changes are merged and CI is green.

## Release checklist

1. Decide the next semantic version.
2. Update `project(... VERSION ...)` in `CMakeLists.txt`.
3. Update user-facing documentation and release notes.
4. Build all targets with warnings enabled and run the complete validation
   suite.
5. Run the full benchmark when algorithmic or performance-sensitive code has
   changed.
6. Install to a temporary prefix and compile an external consumer using
   `find_package(delaunay32 CONFIG REQUIRED)`.
7. Merge the release commit to `main` and confirm CI is green.
8. Create and push an annotated tag:

   ```sh
   git switch main
   git pull --ff-only github main
   git tag -a v0.4.0 -m "Delaunay32 0.4.0"
   git push github v0.4.0
   ```

9. Create a GitHub Release from the tag, with a concise summary of API changes,
   behavior changes, fixes, and known limitations.

Do not move or reuse a published version tag. If a release needs a correction,
make a new patch release instead.
