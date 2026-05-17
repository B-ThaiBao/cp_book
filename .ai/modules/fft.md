# `.ai/modules/fft.md` — FFT / NTT / linear recurrence

Files in [`fft/`](../../fft/):

| File | Top-level name | One-liner |
|---|---|---|
| [`fft.hpp`](../../fft/fft.hpp) | `namespace fft`: `naive_multiplier`, `fft_complex`, `fft_complex_multiplier`, `fft_complex_double_multiplier`, `fft_mod_multiplier`, `fft_numeric_multiplier` (NTT), `multiply_inverser`, `fft_numeric_inverser`, and free fns `naive_multiply`, `fft_complex_multiply`, `fft_complex_double_multiply`, `fft_complex_mod_multiply`, `fft_numeric_multiply`, plus matching `*_inverse`. | Polynomial multiplication & inversion. Auto-selects scratch space via static `scratch_a` / `scratch_b`. |
| [`berlekamp_massey.hpp`](../../fft/berlekamp_massey.hpp) | `berlekamp_massey(A)` | Shortest linear recurrence for a sequence. Returns coefficient vector `C` with `A[i] = C[0]*A[i-1] + ... + C[L-1]*A[i-L]`. |
| [`kitamasa.hpp`](../../fft/kitamasa.hpp) | `kitamasa(A, C, k)` | k-th term of a linear recurrence in O(N² log k). Pair with `berlekamp_massey` for the typical workflow. |

## Module-specific conventions

- Everything sits inside `namespace fft`. Closing `} // namespace fft`.
- `complex<num>` is a hand-rolled struct (faster than
  `std::complex<double>`). Don't replace.
- Roots-of-unity tables are cached in `static std::vector<cnum> roots`
  inside `fft_complex<num>`. `init(k)` is idempotent — safe to call
  before every multiply.
- NTT-side primes & primitive roots live in `fft_numeric_multiplier`.
  Adding a new NTT prime: add it to the prime list and update the
  triple-mod CRT in `fft_numeric_multiply` if you also need
  `fft_complex_mod_multiply`-style precision.
