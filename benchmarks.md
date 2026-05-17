# Performance Benchmarks — cp_book

End-to-end micro-benchmarks for the library, written with Catch2's
`BENCHMARK` macro. All numbers below were captured on a single Linux
x86_64 host (GCC 11, `-O2 -DNDEBUG`, 5 samples per benchmark,
`--rng-seed=1`). Numbers will differ on your machine, but the **relative**
ordering between competing implementations is the interesting signal.

## How to run

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_FLAGS_RELEASE="-O2 -DNDEBUG"
cmake --build build -j

# Run a single benchmark binary:
./build/tests/bench/ds          "[!benchmark]"
./build/tests/bench/ds_advanced "[!benchmark]"
./build/tests/bench/math        "[!benchmark]"
./build/tests/bench/string      "[!benchmark]"
./build/tests/bench/graph       "[!benchmark]"
./build/tests/bench/heap        "[!benchmark]"
./build/tests/bench/fft         "[!benchmark]"
./build/tests/bench/flow        "[!benchmark]"
./build/tests/bench/geo         "[!benchmark]"
./build/tests/bench/dp          "[!benchmark]"
```

Benchmarks are tagged `[!benchmark]` (hidden + non-default) so they are
**not** triggered by a plain `ctest` run. They only fire when you ask for
them explicitly with the tag selector above.

---

## Results (5 samples, mean time per iteration)

### `tests/bench/ds.cpp` — core data structures

| Benchmark | Mean |
|---|---:|
| `seg_tree` (in_order, sum): build + 1e5 point updates + 1e5 range queries on N=1e5 | **9.20 ms** |
| `binary_indexed_tree`: 1e6 point updates + 1e6 prefix queries on N=1e6 | **21.4 ms** |
| `sparse_table<long long>` min: build + 1e6 range-min queries on N=1e6 | **71.6 ms** |
| `range_min_query<int>` (mask + sparse): N=1e6, 1e6 queries | **21.8 ms** |

> `range_min_query` (block + sparse-table-of-block-mins) is ~3.3× faster
> than the plain `sparse_table` on the same workload because each query
> touches only one cache-friendly mask word in the common case.

### `tests/bench/ds_advanced.cpp` — associative & advanced structures

| Benchmark | Mean |
|---|---:|
| `hash_map<long long,int>`: 1e6 insert + 1e6 lookup | **100 ms** |
| `std::unordered_map<long long,int>`: same workload | 190 ms |
| `order_statistic_set<int>` (pb_ds): 5e5 insert + 5e5 `find_by_order` | 657 ms |
| `indexed_set` (bit-trie, universe = 1e6): 1e6 insert + 1e6 contains | **3.23 ms** |
| `std::set<int>`: same workload | 1.43 s |
| `van_emde_boas_tree<20>`: 5e5 insert + 5e5 `find_next`, universe 2²⁰ | **2.75 ms** |
| `disjoint_set_size`: 1e6 nodes, 1e6 random unions | **17.0 ms** |

> The big win is `indexed_set` vs `std::set`: **~440× faster** on a
> bounded-universe contains workload, because the 64-ary bitset trie
> stays entirely in cache and uses popcount/ctz instead of pointer
> chasing. `hash_map` is also a clean ~1.9× over `std::unordered_map`.

### `tests/bench/math.cpp` — number theory & combinatorics

| Benchmark | Mean |
|---|---:|
| `miller_rabin::is_prime<uint64_t>`: 5×10⁴ ~62-bit primality tests | **14.8 ms** (~300 ns / test) |
| `factorizer::factorize<uint64_t>`: 10⁴ ~60-bit Pollard-rho factorizations | **129 ms** (~13 µs / number) |
| `eratos_sieve(1e7)` | 18.5 ms |
| `block_sieve(1e7)` | **2.97 ms** |
| `linear_sieve(1e7)` (also fills SPF table) | 50.8 ms |
| `combinatorics<mint>(1e6+1)` + 1e6 `choose(n,k)` queries | **15.2 ms** |

> `block_sieve` is the clear winner for "I just need the prime list up to
> N" — it's ~6× faster than the classic sieve and ~17× faster than the
> linear sieve (which pays for also producing smallest-prime-factor).

### `tests/bench/string.cpp` — string algorithms

| Benchmark | Mean |
|---|---:|
| `suffix_array` (SAIS) on N=10⁶ over alphabet 26 | **48.0 ms** |
| `z_function` on N=5×10⁶ | **10.7 ms** |
| `prefix_function` (KMP failure table) on N=5×10⁶ | **5.63 ms** |
| `manacher.build` on N=2×10⁶ | **17.5 ms** |

### `tests/bench/graph.cpp` — graph algorithms

| Benchmark | Mean |
|---|---:|
| Tarjan SCC on 10⁵ nodes / 5×10⁵ random directed edges | **26.7 ms** |
| `topo_sort` on a random DAG (10⁵ / 5×10⁵) | **4.92 ms** |
| Kruskal MST on 10⁵ nodes / 5×10⁵ weighted edges | **42.8 ms** |
| `dfs_matching` on 10⁴+10⁴ bipartite with 8×10⁴ edges | **4.08 ms** |
| `hopcroft_karp_matching` on the same instance | 6.10 ms |

> On these random dense instances `dfs_matching` (Kuhn) is actually a
> bit faster than Hopcroft–Karp; HK starts to win on harder, more
> adversarial inputs where the augmenting-path length matters.

### `tests/bench/heap.cpp` — priority queues

| Benchmark | Mean |
|---|---:|
| `radix_heap<int32_t,int>`: 1e6 monotone push + 1e6 pop | **38.9 ms** |
| `std::priority_queue` (min-heap): same workload | 89.7 ms |

> ~2.3× speedup for monotone workloads (Dijkstra-style).

### `tests/bench/fft.cpp` — polynomial multiplication

| Benchmark | Mean |
|---|---:|
| `fft::naive_multiply` 4096 × 4096 | 6.93 ms |
| `fft::fft_complex_multiply` 4096 × 4096 | **120 µs** |
| `fft::fft_complex_multiply` 2¹⁸ × 2¹⁸ | **11.6 ms** |

> ~58× speedup at N=4096 for the complex-FFT path; the 2¹⁸ × 2¹⁸ result
> (~262k × 262k coefficients) finishes in 12 ms.

### `tests/bench/flow.cpp` — max flow

| Benchmark | Mean |
|---|---:|
| `dinic_max_flow` on dense DAG (N=400, avg_deg=20) | 248 µs |
| `hlpp_max_flow` on the same | **241 µs** |
| `dinic_max_flow` on sparse DAG (N=2000, avg_deg=5) | 556 µs |
| `hlpp_max_flow` on the same | **432 µs** |

> Both are essentially instantaneous on these sizes; HLPP edges out
> Dinic on the larger sparse instance.

### `tests/bench/geo.cpp` — computational geometry

| Benchmark | Mean |
|---|---:|
| `convex_hull` on 10⁶ uniformly random points in a square | **188 ms** |
| `convex_hull` on 10⁵ points on a circle (worst case, all on hull) | 7.18 ms |

### `tests/bench/dp.cpp` — DP optimisations

| Benchmark | Mean |
|---|---:|
| `line_multiset` (CHT): 10⁵ line inserts + 10⁵ max queries | **4.01 ms** |

---

## Summary highlights

- **bit-trie `indexed_set` vs `std::set<int>`** on a contains workload:
  ~440× faster (3.2 ms vs 1.4 s).
- **`hash_map` vs `std::unordered_map`**: 1.9× faster.
- **`radix_heap` vs `std::priority_queue`** on monotone keys: 2.3×
  faster.
- **`block_sieve` vs `linear_sieve`** for prime-list generation:
  ~17× faster (because linear sieve also fills the SPF table).
- **`fft_complex_multiply` vs naive convolution** at N=4096: ~58×
  faster.
- **`range_min_query` vs `sparse_table`** on a min query: ~3.3× faster
  thanks to the mask-of-block-mins trick.
