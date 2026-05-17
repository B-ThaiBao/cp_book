# `.ai/modules/geo.md` — computational geometry

Files in [`geo/`](../../geo/):

| File | Top-level name | One-liner |
|---|---|---|
| [`pair_point.hpp`](../../geo/pair_point.hpp) | `template <typename T> struct pair_point` | 2D point with arithmetic, dot, cross helpers. The canonical scalar is integer (long long); floating point only when forced. |
| [`pair_line.hpp`](../../geo/pair_line.hpp) | `pair_line<T>` and helpers | 2D line / segment primitives that compose with `pair_point`. |
| [`convex_hull.hpp`](../../geo/convex_hull.hpp) | `convex_hull(pts, order)` free fn | Andrew monotone chain. Returns `std::array<std::vector<int>, 2>`: `[0]` = ccw vertex sequence, `[1]` = role of each input point (vertex / same as another vertex / colinear on edge / outside). The caller pre-sorts and passes `order`. |
| [`convex_inclusion.hpp`](../../geo/convex_inclusion.hpp) | `convex_inclusion(pts, p, rt)` | Binary-search point-in-convex-polygon. Returns `{inclusion, idx}` with inclusion ∈ {-1=inside, 0=on edge, 1=outside}. |
| [`half_plane_intersect.hpp`](../../geo/half_plane_intersect.hpp) | `half_plane_intersect(lines)` | Returns edge indices of the half-plane intersection in CCW order. Detect unbounded by `size < 3 \|\| cross(front, back) >= 0`. |
| [`minkowski_sum.hpp`](../../geo/minkowski_sum.hpp) | `minkowski_sum(A, B, start)` | Minkowski sum of two convex polygons starting from a given point. |
| [`polygon.hpp`](../../geo/polygon.hpp) | `polygon_area2(...)`, `polygon_centroid(...)`, etc. | Doubled oriented area; centroid for arbitrary (possibly non-convex) polygons. |

## Module-specific conventions

- Default scalar is **integer** (`long long`). Use `__int128` for cross
  products when coordinates exceed ~2³². Avoid `double` unless the
  algorithm genuinely needs it (half-plane intersection in some
  variants does).
- All comparisons go through `cross3(a, b, c)` (or `cross`) — never use
  raw subtraction. The codebase consistently treats `> 0`, `< 0`, `== 0`
  as ccw / cw / colinear.
- Convex polygons are stored **ccw** with no duplicate first/last
  vertex.
- Functions returning ranges (like `convex_hull` and
  `half_plane_intersect`) return **indices into the input**, not the
  values themselves, so the caller controls memory.
