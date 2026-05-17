# `.ai/modules/string.md` — string algorithms

Files in [`string/`](../../string/):

| File | Top-level name | One-liner |
|---|---|---|
| [`prefix_function.hpp`](../../string/prefix_function.hpp) | `prefix_function(s)` + `prefix_automaton` | KMP failure function and the matching automaton. |
| [`z_function.hpp`](../../string/z_function.hpp) | `z_function(s)` | Z-array. |
| [`suffix_array.hpp`](../../string/suffix_array.hpp) | `struct suffix_array : std::vector<int32_t>` | SAIS (induced sort) in O(N). **Caller must compress** input to small non-negative ints first. Construct with `suffix_array sa(s);`. Pair with `range_min_query` for LCP queries (`tests/string/suffix_array.cpp` includes `ds/range_min_query.hpp`). |
| [`manacher.hpp`](../../string/manacher.hpp) | `struct manacher` | Palindromic radius array of length 2N-1. Build via `manacher m(s);` or `manacher m; m.build(s);`. `m.is_palindrome(l, r)` answers the standard interval query. Iterate raw radii via `m.begin()/m.end()` (the struct inherits from a vector). |
| [`prefix_tree.hpp`](../../string/prefix_tree.hpp) | `struct prefix_tree` | Trie. Construct with `prefix_tree p(num_childs);`, then `p.make_root(); p.insert(s); p.find(s); p.erase(s);`. |
| [`prefix_power_hash.hpp`](../../string/prefix_power_hash.hpp) | `template <typename T> struct prefix_power_hash`, plus `base_power<T>::base` and aggregate types `pairnum<T>` / `arraynum<T, K>` | Polynomial rolling hash. `T` is typically `modnum<…>` or a `pairnum`/`arraynum` of those. **Set `base_power<T>::base = T(some_value);` once before use**; clear/reset if you switch base mid-program. |

## Module-specific conventions

- All string algorithms take **arbitrary `Vector`-like input** (the
  caller can pass `std::string`, `std::vector<int>`, `std::vector<int8_t>`, …)
  via template type erasure. Don't hard-code `std::string`.
- For hashing, the codebase prefers **multi-modulus pairs** (e.g.
  `pairnum<hashnum>`) over a single 64-bit modulus, because the
  multiplier policy is parameterised. The `prefix_power_hash` header
  enumerates the supported `T` instantiations at the top.
- `suffix_array` requires the alphabet to start at a non-negative
  number and be densely packed. Use `compressor` (see misc/) or a
  manual subtract-min.
- `manacher` exposes both the raw radius array (vector-like iteration)
  and the high-level `is_palindrome(l, r)` query — prefer the latter
  unless you really need radii.
