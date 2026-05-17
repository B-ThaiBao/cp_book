# `.ai/modules/math.md` — number theory & matrix

Files in [`math/`](../../math/):

| File | Top-level name | One-liner |
|---|---|---|
| [`combinatorics.hpp`](../../math/combinatorics.hpp) | `template <typename num> struct combinatorics : std::vector<std::array<num, 2>>` | Factorial + inverse factorial table over a `modnum`-like field. `C.fact(n)`, `C.inv_fact(n)`, `C.choose(n,k)`, `C.permute(n,k)`, `C.catalan(...)`. |
| [`derangement.hpp`](../../math/derangement.hpp) | `derangement(N)` | !N table. |
| [`diophantine.hpp`](../../math/diophantine.hpp) | `extended_gcd`, `diophantine` | Solves `a*x + b*y = c`. |
| [`div.hpp`](../../math/div.hpp) | `divisors(n)`, etc. | Divisor enumeration. |
| [`euler_totient_function.h`](../../math/euler_totient_function.h) | `phieuler(n)`, `range_phieuler(n)` | Euler totient (single value O(√n) / range O(n log log n)). **Header has `.h` extension, not `.hpp`** — leave it alone, tests already handle it. |
| [`factorizer.hpp`](../../math/factorizer.hpp) | `namespace factorizer`: `linear_sieve(N)`, `factorize<T>(x)`, `is_prime<T>(x)` | Pollard-rho on top of a precomputed linear sieve. **Must call `factorizer::linear_sieve(N)` once** before `factorizer::factorize<uint64_t>(x)`. |
| [`legendre_power.hpp`](../../math/legendre_power.hpp) | `legendre_power(n, p)` | Legendre's formula: largest k with p^k | n!. |
| [`matrix.hpp`](../../math/matrix.hpp) | `template <typename T, size_t ROW, size_t COL> struct matrix : std::array<T, ROW*COL>` | **Fixed-size** matrix. Dimensions are template parameters, not runtime. Supports `+`, `-`, `*`, scalar ops, `inverse`, etc. |
| [`miller_rabin.hpp`](../../math/miller_rabin.hpp) | `namespace miller_rabin`: `is_prime<T>(N)` | Deterministic Miller–Rabin over 64-bit ints. Free function `miller_rabin::is_prime<uint64_t>(N)`. |
| [`sieve.hpp`](../../math/sieve.hpp) | `eratos_sieve(N)`, `block_sieve(K)`, `struct linear_sieve` | Three separate APIs: classic O(N log log N) returning `vector<bool>`; cache-friendly block sieve returning `vector<int>` of primes (fastest); linear sieve struct exposing both `primes` and SPF table `ls[x]`. `range_sieve` in this file is **broken** (initialises `is_primes` inverted) and intentionally not tested. |
| [`sqrt.hpp`](../../math/sqrt.hpp) | `isqrt(n)`, `tonelli_shanks(...)` | Integer sqrt + modular sqrt (Tonelli–Shanks). |

## Module-specific conventions

- `combinatorics<num>(N)` precomputes `fact[0..N-1]` and `inv_fact[0..N-1]`.
  Always pass `N` larger than the max `n` you'll query — there is a
  `_GLIBCXX_DEBUG`-only bound check.
- `matrix<T, R, C>` requires R and C at compile time. Need runtime
  dims? Don't try to extend this header — write a separate
  `dynamic_matrix` (none currently exists).
- Factorizer is namespace-style; Miller-Rabin is namespace-style.
  Combinatorics, matrix, sieve are struct/free-function style. Don't
  mix the patterns.
- The `.h` extension on `euler_totient_function.h` is an oddity. Don't
  rename to `.hpp` without a reason.
