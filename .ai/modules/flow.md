# `.ai/modules/flow.md` — flows

Files in [`flow/`](../../flow/):

| File | Top-level name | One-liner |
|---|---|---|
| [`flow_graph.hpp`](../../flow/flow_graph.hpp) | `template <typename flow_t> struct flow_graph` | Edge-list graph with paired forward/back edges (`id` / `id ^ 1`). `add_edge(u, v, cap)`. Tracks `eps` for floating-point flows. All max-flow / min-cost-flow algorithms in this folder operate on this template. |
| [`max_flow.hpp`](../../flow/max_flow.hpp) | `dinic_max_flow<flow_t>`, `hlpp_max_flow<flow_t>` | Dinic (O(V²·E)) and HLPP (O(V²·√E)). Both expose `max_flow(g, s, t)` and a `min_cut(g)` reachability vector. |
| [`dijkstra_min_cost_flow.hpp`](../../flow/dijkstra_min_cost_flow.hpp) | `dijkstra_min_cost_flow<flow_t, cost_t>` | Min-cost max-flow via potentials + Dijkstra. **Requires non-negative initial costs** (or use SPFA potential init). Uses pb_ds `__gnu_pbds::priority_queue` for fast decrease-key. |
| [`network_simplex.hpp`](../../flow/network_simplex.hpp) | `network_simplex<flow_t, cost_t>` | Network simplex for min-cost flow. **Patch history**: `random_device{}()` → `std::random_device{}()` (line 281). |

## Module-specific conventions

- All max-flow / mcmf algorithms are **templated over the user's flow
  graph type** (parameter name `Network`), not over `flow_graph`
  directly. So `bfs(const Network& g, …)` works with any flow-graph-like
  struct with the same shape.
- `g.eps` is the comparison tolerance. For integer `flow_t` this is
  zero; for `double` set it before use (e.g. `g.eps = 1e-9;`).
- Re-build the `flow_graph` between max-flow runs unless you genuinely
  want a residual graph. The benchmarks in [`tests/bench/flow.cpp`](../../tests/bench/flow.cpp)
  follow this pattern.
- `min_cut` only works **after** `max_flow` (it inspects `depth[]`).
