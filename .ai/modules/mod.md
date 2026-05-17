# `.ai/modules/mod.md` — modular arithmetic

Files in [`mod/`](../../mod/):

| File | Top-level name | One-liner |
|---|---|---|
| [`modnum.hpp`](../../mod/modnum.hpp) | `template <typename T, const T num> struct constant`, `template <typename T> struct inconstant`, `template <typename T> struct naive_multiplier`, `template <typename T> struct barrett_multiplier`, `template <typename T, typename multiplier> struct modnum` | Policy-based modular integer. `T` = compile-time `constant<int, MOD>` or runtime `inconstant<int>`; `multiplier` = `naive_multiplier<int>` or `barrett_multiplier<int>`. |

## Canonical type aliases

```cpp
// Compile-time prime, naive 64-bit mulmod (fastest for fixed MOD):
using mint  = modnum<constant<int, 1000000007>, naive_multiplier<int>>;
using mint2 = modnum<constant<int,  998244353>, naive_multiplier<int>>;

// Runtime prime, Barrett reduction:
using dmint = modnum<inconstant<int>, barrett_multiplier<int>>;
// Set with: inconstant<int>::set_value(some_runtime_mod);
```

## Module-specific conventions

- **No Montgomery multiplier** is included; the header explicitly
  comments where to copy one from if you really need it
  (NyaanNyaan/library).
- The `multiply()` policy is selected via `std::enable_if` on
  `sizeof(K) == 4` vs `sizeof(K) == 8`. Keep both overloads when you
  add a new multiplier.
- Extract the underlying value with `mint()` (function-call operator)
  or `static_cast<int>(m)`. There is **no `.value()` method**.
- `is_zero(m)`, `abs(m)`, `to_string(m)`, and `__print(m)` are provided
  as friend overloads — they are picked up by ADL by the rest of the
  library (e.g. tests, debug.hpp).
- Streaming operators `<<` / `>>` are templated on the stream type, so
  they work with both `std::ostream` and the custom `misc/io.hpp`
  writers.
