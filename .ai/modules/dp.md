# `.ai/modules/dp.md` — DP optimisations

Files in [`dp/`](../../dp/):

| File | Top-level name | One-liner |
|---|---|---|
| [`convex_container.hpp`](../../dp/convex_container.hpp) | (no doc block — read source) | Monotonic / convex container helper used by CHT-style DPs. |
| [`line_multiset.hpp`](../../dp/line_multiset.hpp) | `struct line<T>`, `struct line_multiset<Line, Comp>` | Dynamic CHT via `std::multiset`. Query with `*lower_bound(x)`; `(*it)[0] * x + (*it)[1]`. For min, negate slope+intercept. |

## Module-specific conventions

- `line<T>` inherits from `std::array<T, 2>` and adds a **mutable** `x`
  storing the intersection abscissa with the next line. Treat that field
  as a private implementation detail when reading — but it is what makes
  the `lower_bound(x)` query branch-free.
- Insertion is done via `add_line(a, b)`. Do not insert raw lines —
  the `intersect()` helper relies on the field being initialised to 0.
- `line_multiset` chooses integer or floating-point `divide()` via
  `std::enable_if` on the value type. Keep both overloads if you add a
  new numeric backing type.

## What to copy when adding a new dp optimisation

- Mirror the `line<T>` / `line_multiset<...>` split: the **element**
  type carries the comparison hooks (`friend operator <` both
  directions), and the **container** wraps an STL one via inheritance.
- Provide `add_*` and let `lower_bound`/`upper_bound` come for free
  from the base STL container.
