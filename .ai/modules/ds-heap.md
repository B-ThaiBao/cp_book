# `.ai/modules/ds-heap.md` — heaps

Files in [`ds/heap/`](../../ds/heap/):

| File | Top-level name | One-liner |
|---|---|---|
| [`radix_heap.hpp`](../../ds/heap/radix_heap.hpp) | `template <typename key_t, typename val_t> struct radix_heap` | Monotone radix heap. `h.emplace(k, v)` / `h.push({k, v})`. `h.top()` and `h.pop()` return `std::pair<key_t, val_t>`. Used inside the included Dijkstra reference impl at the bottom of the file. Keys MUST be monotonically non-decreasing across pops. Use `int32_t` or `int64_t` keys — `log_2` is only overloaded for those (mixing in `uint32_t` causes an ambiguity). |
| [`skew_heap.hpp`](../../ds/heap/skew_heap.hpp) | `template <typename T> struct skew_heap_node_base`, `struct skew_heap_node : skew_heap_node_base<…>`, `struct dijkstra_skew_heap_node` | Intrusive meldable heap. Users derive a node struct, allocate nodes themselves, and call `merge(a, b)` / `pop(root)`. **Not** an STL-shaped container — there is no `push_back`. |

## Module-specific conventions

- Radix heap is a 2-template-argument template
  (`<key_t, val_t>`). It does **not** take a 3rd "container" argument
  like `std::priority_queue`. Use `int32_t` keys for monotone Dijkstra
  workloads where keys fit; `int64_t` otherwise.
- Skew heap is intrusive — never wrap it in a generic
  `priority_queue`-like façade; that defeats the purpose. Look at the
  embedded `dijkstra_skew_heap_node` example before adding a new use
  site.
