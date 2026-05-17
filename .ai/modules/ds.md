# `.ai/modules/ds.md` — data structures (root level)

Files in [`ds/`](../../ds/) (subfolders documented separately:
[ds-bst.md](ds-bst.md), [ds-heap.md](ds-heap.md)).

| File | Top-level name | One-liner |
|---|---|---|
| [`binary_indexed_tree.hpp`](../../ds/binary_indexed_tree.hpp) | `struct binary_indexed_tree<T> : std::vector<T>` | Generic Fenwick. Update: `for (auto& x : bit.suffix(i)) x += v;`. Query: `for (auto& x : bit.prefix(r)) s += x;`. Also: `build_prefix`, `build_suffix`, `lower_bound`, `upper_bound`, and `range_update_sum_query` (two-tree). |
| [`cartesian_tree.hpp`](../../ds/cartesian_tree.hpp) | `cartesian_tree(...)` free fn | Returns parent array. Comparator chooses leftmost/rightmost min/max tree (`std::less<>`, `std::less_equal<>`, …). |
| [`disjoint_set.hpp`](../../ds/disjoint_set.hpp) | `struct disjoint_set_size`, `struct disjoint_set_ancestor` | DSU with size or with ancestor pointer. `merge(a, b)` returns true on real merge. **There is no plain `disjoint_set`** — anyone writing `disjoint_set` means `disjoint_set_size`. |
| [`hash_map.hpp`](../../ds/hash_map.hpp) | `hash_map<K, V>` | pb_ds `gp_hash_table` with a custom anti-collision hash. Linux-only. |
| [`hilbert_curve.hpp`](../../ds/hilbert_curve.hpp) | `hilbert_order(...)` | Hilbert-curve order for Mo's algorithm. |
| [`indexed_set.hpp`](../../ds/indexed_set.hpp) | `struct indexed_set` | 64-ary bitset trie over a bounded `int` universe. `insert/erase/find/find_next/find_prev`. Built via `build(M)` or `build(M, predicate)`. Extremely fast vs `std::set` (~440× on contains workloads). |
| [`order_statistic.hpp`](../../ds/order_statistic.hpp) | `order_statistic_set<K>`, `order_statistic_map<K,V>` | pb_ds `tree<>` with `find_by_order` / `order_of_key`. Linux-only. |
| [`range_min_query.hpp`](../../ds/range_min_query.hpp) | `struct range_min_query<T, Comp = std::less<T>>` | Mask-of-block-mins + sparse table on block mins. O(N)/O(1). `rmq.range_query(l, r)` returns `{idx, val}` for INCLUSIVE range. ~3× faster than `sparse_table` for min. |
| [`seg_tree.hpp`](../../ds/seg_tree.hpp) | `namespace seg_tree`: `in_order_tree`, `circular_tree`, `point_t`, `range_t` | Layout-based segment tree. User defines their own `std::vector<Node>` and uses `layout.point(i)` / `layout.range(l, r).for_each(...)` / `for_ancestor_up` / `for_ancestor_down`. Supports lazy propagation by writing `update` + `downdate` callbacks. See header doc-block for the four standard patterns. |
| [`sparse_table.hpp`](../../ds/sparse_table.hpp) | `struct sparse_table<T>`, `struct disjoint_sparse_table<T>` | Generic associative-op sparse table. `st(0, i) = v;` for leaves, then `st.build(merge)`. Iterate with `st.range(l, r)` or `st.for_range(l, r, fn)`. |
| [`swag.hpp`](../../ds/swag.hpp) | SWAG variants (deque-based and sqrt-decomposed) | Sliding window aggregation: O(1) amortised push/pop with a monoid op. |
| [`van_emde_boas_tree.hpp`](../../ds/van_emde_boas_tree.hpp) | `template <int LOG_U> struct van_emde_boas_tree` | Static-universe successor structure. Universe = `1 << LOG_U`. O(log log U) ops. |
| [`wavelet.hpp`](../../ds/wavelet.hpp) | `wavelet_tree` (user-extended node) | User defines a node struct with `std::array<int, 2> ch` and a `bit_vector`. Library provides the build / query skeleton. |

## Module-specific conventions

- BIT/sparse-table/segment-tree are all **iterator/callback driven**.
  Do not write a "do this with an int loop" wrapper — use the provided
  `range`, `prefix`, `suffix`, `for_each`, `for_range` callbacks.
- Many structs inherit from `std::vector` to inherit constructors and
  iteration. Keep that pattern for any new vector-shaped DS.
- Bounds asserts go behind `#ifdef _GLIBCXX_DEBUG`. Never plain
  `assert()` in a hot path.
- `seg_tree::in_order_tree` uses 2·N internal storage; the user
  allocates the data vector externally.
