# Third-party software

The `delaunay32` library target has no third-party source dependency.

The optional benchmark and validation targets use
[`delaunator-cpp`](https://github.com/delfrrr/delaunator-cpp) as a Git
submodule at `third_party/delaunator-cpp`. It is not compiled into, linked
into, or installed with the `delaunay32` library. The release is pinned to
commit `c1521f6e879881232dcddabd6c2ddb6187e8714b`.

`delaunator-cpp` is distributed under the MIT License and includes code
derived from Mapbox Delaunator, which is distributed under the ISC License.
The complete license and notice texts remain in the submodule:

- `third_party/delaunator-cpp/LICENSE`
- `third_party/delaunator-cpp/THIRD-PARTY-NOTICES`
