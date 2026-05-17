# `.ai/modules/graph.md` — graph algorithms

Files in [`graph/`](../../graph/):

| File | Top-level name | One-liner |
|---|---|---|
| [`graph.hpp`](../../graph/graph.hpp) | `template <typename T> struct graph`, `undigraph<T>`, `digraph<T>`, `condi_undigraph<T,F>`, `condi_digraph<T,F>` | Edge-list-based graph templates. User defines edge struct (`struct E { int from, to; … }`). `g.add_edge(u, v, …)` forwards trailing args to the edge constructor. `condi_*` adds a per-edge ignore predicate. |
| [`ancestor.hpp`](../../graph/ancestor.hpp) | `rmq_ancestor` | LCA via Euler tour + O(N)/O(1) RMQ. Plug-in for path queries. |
| [`sparse_ancestor.hpp`](../../graph/sparse_ancestor.hpp) | `sparse_ancestor` | Binary-lifting LCA + path aggregate. Build via `combine(i, u)` nested loops. |
| [`build_bipartite.hpp`](../../graph/build_bipartite.hpp) | `build_bipartite(g, side)` | Two-colour a graph if bipartite, otherwise return false. |
| [`centroid_tree.hpp`](../../graph/centroid_tree.hpp) | `centroid_tree` | Builds centroid decomposition. Height O(log N); `path(A, B)` = `A → lca(A, B) → B`. |
| [`euler_tour.hpp`](../../graph/euler_tour.hpp) | `euler_tour_edge` | DFS Euler tour over edges. `et.dfs(g, u)` populates `euler` (visited-edge sequence) and `loc` (first index per node). |
| [`eulerian.hpp`](../../graph/eulerian.hpp) | Eulerian circuit / path utilities | (see source). |
| [`matching.hpp`](../../graph/matching.hpp) | `dfs_matching`, `hopcroft_karp_matching` | Bipartite matching. Build with `mat(N, M)` and call `mat.max_match(g)` on a `digraph<E>` with left→right edges. |
| [`min_vertex_cover.hpp`](../../graph/min_vertex_cover.hpp) | `min_vertex_cover` | **KNOWN BUG**: loop bounds use `mat.mate[1][i]` indexed by `i < mat.N` — incorrect for `N != M`. Tagged `[.skip]` in tests. Do not depend on it. |
| [`scc_comp.hpp`](../../graph/scc_comp.hpp) | `scc_comp(g, cnt)` | Tarjan SCC. Returns `std::vector<int>` component id per node; sets `cnt` to total SCC count. **Component ids are emitted in reverse-topological order — id 0 is a sink.** |
| [`span_forest.hpp`](../../graph/span_forest.hpp) | `dfs_span_forest` | DFS span tree (span edges + back edges). Exposes `par`, etc. |
| [`span_tree.hpp`](../../graph/span_tree.hpp) | `build_span_tree(g, sum)` | Kruskal MST. Returns chosen edge ids; writes total weight into `sum`. **Patch history**: `disjoint_set` → `disjoint_set_size`. |
| [`topo_sort.hpp`](../../graph/topo_sort.hpp) | `topo_sort(g)` | Kahn topo sort; empty vector if not a DAG. |
| [`tree_diameter.hpp`](../../graph/tree_diameter.hpp) | `tree_diameter` | Two-BFS diameter. After `bfs()`, you get `{x, y}`, plus distances and parent-edges rooted at one endpoint. |
| [`twosat_graph.hpp`](../../graph/twosat_graph.hpp) | `twosat_graph` | 2-SAT solver. Negated variables = `~ x`. `sat()` returns assignment or empty vector. **Patch history**: SCC-id comparison direction fixed (`<` is correct given this header's emission-order SCC numbering with id 0 = sink). |

## Module-specific conventions

- Edge struct convention: `struct E { int from, to; … }`. `add_edge`
  forwards trailing args via parameter pack.
- Undirected `g(u, i)` returns "the *other* endpoint of edge i when
  coming from u" via `u ^ from ^ to`. Use this — don't manually pattern
  match.
- SCC component ids in `scc_comp` are in **reverse topological order**
  (sink first). This is unusual and is the source of the `twosat` `<`
  vs `>` patch — be sure when writing new code on top of `scc_comp` to
  account for it.
- For weighted edges, use `struct WE { int from, to; long long cost; }`.
