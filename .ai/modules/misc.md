# `.ai/modules/misc.md` — utilities

Files in [`misc/`](../../misc/):

| File | Top-level name | One-liner |
|---|---|---|
| [`compressor.hpp`](../../misc/compressor.hpp) | `compressor` | Compress a vector of values to 0-indexed ids. **Caller must sort first** — the compressor does not sort. |
| [`io.hpp`](../../misc/io.hpp) | Custom IO helpers | Fast read / print. |
| [`reverse_args.hpp`](../../misc/reverse_args.hpp) | `reverse_args(f)` | Wraps a callable to reverse its argument order. |
| [`ska_sort.hpp`](../../misc/ska_sort.hpp) | `ska_sort` | Malte Skarupke's radix sort — fast for integer/composite keys. |
| [`tensor.hpp`](../../misc/tensor.hpp) | `tensor<T, RANK>` | N-dimensional dense array helper with arbitrary rank. Useful for DP states. |
| [`y_combinator.hpp`](../../misc/y_combinator.hpp) | `std::y_combinator` | Recursive lambda wrapper (in `namespace std`, deliberately). |

## Module-specific conventions

- `y_combinator` is **injected into `namespace std`** (anti-pattern in
  normal C++ but consistent with competitive-programming convention
  here). Do not "fix" it.
- IO helpers in `io.hpp` are not used by every problem; many templates
  in `template/` opt for plain `std::cin` with `sync_with_stdio(false)`.
- `ska_sort` is opt-in — `std::sort` is fine for the typical N ≤ 10⁶
  workload.
