# Usage.md — `cp_book`

A practical, header-by-header guide to using every component of this
library. Every snippet assumes you already `#include`d the matching
header and are inside `int main()`. All examples are derived from the
test suite under [`tests/`](tests/), which you can consult for richer
end-to-end scenarios.

> **Minimum build flags**: `g++ -std=c++20 -O2 -I.` (run from the repo
> root). A few headers (`hash_map.hpp`, `order_statistic.hpp`,
> `dijkstra_min_cost_flow.hpp`) depend on `__gnu_pbds`, so they only
> build on GCC / Linux.
>
> **Type convention**: throughout this guide we prefer the
> fixed-width aliases `int32_t` / `int64_t` (defined in `<cstdint>`)
> over `int` / `long long`. The library itself accepts any integral
> type; the choice is the caller's.

---

## Table of contents

- [dp/](#dp) — DP optimisations
- [ds/](#ds) — core data structures
- [ds/bst/](#dsbst) — balanced BSTs
- [ds/heap/](#dsheap) — priority queues
- [fft/](#fft) — FFT / NTT / linear recurrences
- [flow/](#flow) — max-flow and min-cost flow
- [geo/](#geo) — computational geometry
- [graph/](#graph) — graph algorithms
- [math/](#math) — number theory and matrices
- [misc/](#misc) — utilities
- [mod/](#mod) — modular arithmetic
- [string/](#string) — string algorithms
- [template/](#template) — contest skeletons
- [General tips](#general-tips)

---

## dp/

### [`dp/line_multiset.hpp`](dp/line_multiset.hpp) — dynamic convex hull trick

```cpp
#include "dp/line_multiset.hpp"

// Maintain a set of lines y = a*x + b and query the upper envelope.
line_multiset<line<int64_t>> hull;
hull.add_line(2, 5);              // y = 2x + 5
hull.add_line(-1, 10);            // y = -x + 10
hull.add_line(3, -4);
int64_t best_max = hull.query(7); // max over all added lines at x = 7

// For a min hull, negate both slope and intercept on insert,
// then negate the query result.
```

> `line<T>` inherits `std::array<T, 2>` (slope is `[0]`, intercept is
> `[1]`); always insert through `add_line(a, b)` — the intersect helper
> relies on private internal state.

### [`dp/convex_container.hpp`](dp/convex_container.hpp)

Currently empty (placeholder for a monotonic CHT). Use
`line_multiset` above for the dynamic case.

---

## ds/

### [`ds/binary_indexed_tree.hpp`](ds/binary_indexed_tree.hpp) — Fenwick tree

```cpp
#include "ds/binary_indexed_tree.hpp"

binary_indexed_tree<int64_t> bit(N);   // N elements, 0-indexed

// Point update at index i:
for (auto& x : bit.suffix(i)) x += v;

// Prefix sum over [0, r):
int64_t s = 0;
for (auto& x : bit.prefix(r)) s += x;

// Build from an existing array in O(N) instead of N * O(log N):
std::vector<int64_t> a = {1, 2, 3, 4, 5};
binary_indexed_tree<int64_t> bit2(a.begin(), a.end());
bit2.build_prefix();                   // or build_suffix(), depending on direction

// Find the smallest index whose prefix sum is >= k (requires build_prefix):
int32_t idx = bit2.lower_bound(k);

// Range update + range sum (two-tree trick) — use a dedicated wrapper:
range_update_sum_query<int64_t> rsq(N);
rsq.range_update(l, r, v);             // add v to [l, r)
int64_t t = rsq.range_query(l, r);
```

### [`ds/cartesian_tree.hpp`](ds/cartesian_tree.hpp)

```cpp
#include "ds/cartesian_tree.hpp"

std::vector<int32_t> a = {3, 1, 4, 1, 5, 9, 2, 6};
// par[i] = parent of i in the Cartesian tree; par[root] = -1.
auto par      = cartesian_tree(a);                          // leftmost-min tree (std::less)
auto par_le   = cartesian_tree(a, std::less_equal<>());     // duplicate-friendly
auto par_max  = cartesian_tree(a, std::greater<>());        // max tree
```

### [`ds/disjoint_set.hpp`](ds/disjoint_set.hpp) — DSU

```cpp
#include "ds/disjoint_set.hpp"

disjoint_set_size dsu(N);              // there is NO bare `disjoint_set`
bool merged   = dsu.merge(u, v);       // true if the two roots actually merged
int32_t root  = dsu.find(u);
int32_t sz    = dsu.size(root);

// Variant that maintains an ancestor pointer (useful for offline LCA):
disjoint_set_ancestor dsa(N);
dsa.merge(u, v);

// Bipartite DSU (detect odd cycles, check side):
disjoint_set_bipartite dsb(N);
dsb.merge(u, v);
bool bip      = dsb.is_bipartite();
bool same_lr  = dsb.is_same(u, v);
```

### [`ds/hash_map.hpp`](ds/hash_map.hpp) — pb_ds hash with anti-collision *(Linux-only)*

```cpp
#include "ds/hash_map.hpp"

hash_map<int64_t, int32_t> mp;
mp[42] = 7;
if (mp.find(42) != mp.end()) { /* ... */ }

hash_set<int32_t> st;
st.insert(3);
```

### [`ds/hilbert_curve.hpp`](ds/hilbert_curve.hpp) — Mo's ordering

```cpp
#include "ds/hilbert_curve.hpp"

// Sort offline queries (l, r) by Hilbert order for Mo's algorithm.
struct query { int32_t l, r, id; int64_t ord; };
for (auto& q : qs) q.ord = hilbert_order(q.l, q.r, LOG_N, 0);
std::sort(qs.begin(), qs.end(),
          [](const auto& a, const auto& b) { return a.ord < b.ord; });
```

### [`ds/indexed_set.hpp`](ds/indexed_set.hpp) — 64-ary bitset trie

```cpp
#include "ds/indexed_set.hpp"

indexed_set s(M);                  // universe = [0, M)
s.insert(7);
s.erase(7);
bool has    = s.find(7);
int32_t nxt = s.find_next(10);     // smallest i >= 10 in s, or -1
int32_t prv = s.find_prev(10);

// Build from a predicate in O(M / 64):
s.build(M, [](int32_t i) { return i % 2 == 0; });

// Ascending iteration:
s.for_each([](int32_t i) { /* ... */ });
```

### [`ds/order_statistic.hpp`](ds/order_statistic.hpp) — pb_ds order-statistic tree *(Linux-only)*

```cpp
#include "ds/order_statistic.hpp"

order_statistic_set<int32_t> st;
st.insert(5); st.insert(2); st.insert(8);
auto it       = st.find_by_order(0);     // 0-th smallest -> 2
int32_t rank  = st.order_of_key(5);      // number of elements < 5 -> 1

order_statistic_map<int32_t, std::string> mp;
mp[10] = "ten";
```

### [`ds/range_min_query.hpp`](ds/range_min_query.hpp)

```cpp
#include "ds/range_min_query.hpp"

std::vector<int32_t> a = {3, 1, 4, 1, 5, 9, 2, 6};
range_min_query<int32_t> rmq(a.begin(), a.end());         // leftmost min
auto [idx, val] = rmq.range_query(l, r);                  // INCLUSIVE [l, r]

// Rightmost min: pass std::less_equal<>.
range_min_query<int32_t, std::less_equal<>> rmq2(a.begin(), a.end());

// Range max: pass std::greater<>.
range_min_query<int32_t, std::greater<>> mx(a.begin(), a.end());
```

### [`ds/seg_tree.hpp`](ds/seg_tree.hpp) — generic segment tree

You own the `std::vector<Node>` and the `update` / `downdate`
callbacks. The header only provides the layout and iterators. Canonical
pattern:

```cpp
#include "ds/seg_tree.hpp"

// 1) Point update + range sum.
struct sum_seg {
	seg_tree::in_order_tree layout;
	std::vector<int64_t> val;                  // size 2*N
	sum_seg(const int32_t& N) : layout(N), val(2 * N, 0) {}
	void point_update(const int32_t& i, const int64_t& v) {
		seg_tree::point_t p = layout.point(i);
		val[p] = v;
		p.for_ancestor_up([&](seg_tree::point_t q) {
			val[q] = val[q.c(0)] + val[q.c(1)];
		});
	}
	int64_t range_sum(const int32_t& l, const int32_t& r) {   // [l, r)
		int64_t s = 0;
		layout.range(l, r).for_each([&](seg_tree::point_t p) { s += val[p]; });
		return s;
	}
};

// 2) Range update + range sum (with lazy propagation): see
//    tests/ds/seg_tree.cpp for a full lazy_seg implementation
//    (apply_node / push / pull / range_add / range_sum).

// circular_tree: a different layout that is faster for power-of-two N.
```

### [`ds/sparse_table.hpp`](ds/sparse_table.hpp)

```cpp
#include "ds/sparse_table.hpp"

sparse_table<int32_t> st(N);
for (int32_t i = 0; i < N; ++ i) st(0, i) = a[i];
st.build([](const int32_t& x, const int32_t& y) { return std::min(x, y); });
int32_t m = st.range(l, r).value;             // [l, r), idempotent operations only

// disjoint_sparse_table: works for non-idempotent operations (sum, product, ...).
disjoint_sparse_table<int64_t> dst(N);
```

### [`ds/swag.hpp`](ds/swag.hpp) — sliding-window aggregation

```cpp
#include "misc/reverse_args.hpp"                  // MUST be included BEFORE swag.hpp
#include "ds/swag.hpp"

auto f_min = [](const int32_t& a, const int32_t& b) { return std::min(a, b); };

// Stack flavour (push_back / pop_back):
auto sk = make_swag_stack<int32_t>(f_min);
sk.push_back(5); sk.push_back(2); int32_t m = sk.query();

// Queue flavour (push_back / pop_front):
auto q = make_swag_queue<int32_t>(f_min);

// Deque flavour (push and pop on both ends):
auto d = make_swag_deque<int32_t>(f_min);
d.push_back(3); d.push_front(1); int32_t cur = d.query();
```

### [`ds/van_emde_boas_tree.hpp`](ds/van_emde_boas_tree.hpp)

```cpp
#include "ds/van_emde_boas_tree.hpp"

van_emde_boas_tree<20> veb;        // universe = 1 << 20
veb.insert(123);
veb.erase(123);
bool has    = veb.find(123);
int32_t nxt = veb.find_next(100);  // successor, or -1
int32_t prv = veb.find_prev(100);
```

### [`ds/wavelet.hpp`](ds/wavelet.hpp) — wavelet tree (user-extended node)

```cpp
#include "ds/wavelet.hpp"

// The bit-vector primitive is reusable on its own:
bit_vector bv(N);
bv.set(i);                         // does not auto-build
bv.build();                        // O(N / 64)
int32_t ones = bv.rank(l, r);

// For the wavelet tree itself, define a node struct holding
// std::array<int32_t, 2> ch and a bit_vector b — see the header
// doc-block for the build / query skeleton.
```

---

## ds/bst/

Every BST in this library is **intrusive**: derive your own node from
the provided base and call free functions such as `merge` / `split` /
`link` / `cut`.

### [`ds/bst/treap.hpp`](ds/bst/treap.hpp) — implicit treap

```cpp
#include "ds/bst/treap.hpp"

struct my_node : public treap_node_base<my_node> {
	int64_t val = 0, sum = 0;
	void do_downdate() {}              // no extra lazy data
	void do_update() {
		sum = val;
		if (c[0]) sum += c[0]->sum;
		if (c[1]) sum += c[1]->sum;
	}
};

std::vector<my_node> pool(N);
my_node* root = nullptr;
for (int32_t i = 0; i < N; ++ i) {
	pool[i].val = a[i];
	root = merge(root, &pool[i]);      // amortised push_back
}

// Implicit (by-index) access:
my_node* x = find_implicit(root, k);

// Split into [0, k) | [k, end):
auto [l, r] = split(root, k);
root = merge(l, r);

// To reverse [l, r) you need a flip flag on the node — see
// tests/ds/bst/treap.cpp.
// IMPORTANT: helpers like linearize / collect_inorder must call
// u -> downdate() (NOT do_downdate()) so that lazy state propagates.
```

### [`ds/bst/splay_tree.hpp`](ds/bst/splay_tree.hpp)

```cpp
#include "ds/bst/splay_tree.hpp"

struct my_splay : public splay_node_base<my_splay> {
	int32_t val = 0;
	void do_downdate() {}
	void do_update() {}
};

std::vector<my_splay> pool(N);
my_splay* root = build(pool.data(), pool.data() + N);

my_splay* p = find_implicit(root, k);
auto [l, r] = split(root, k);
root = merge(l, r);
root = erase(root, k);
```

### [`ds/bst/link_cut_tree.hpp`](ds/bst/link_cut_tree.hpp)

```cpp
#include "ds/bst/link_cut_tree.hpp"

struct lct : public link_cut_tree_node_base<lct> {
	int32_t id = 0;
	void do_downdate() {}
	void do_update() {}
};

std::vector<lct> pool(N);
for (int32_t i = 0; i < N; ++ i) pool[i].id = i;

link(&pool[u], &pool[v]);              // connect if not already connected
cut(&pool[u], &pool[v]);
auto* r = find_lct_root(&pool[u]);     // root of the component
make_lct_root(&pool[u]);               // evert
```

---

## ds/heap/

### [`ds/heap/radix_heap.hpp`](ds/heap/radix_heap.hpp) — monotone radix heap

```cpp
#include "ds/heap/radix_heap.hpp"

radix_heap<int64_t, int32_t> h;        // key (int32_t or int64_t), payload
h.emplace(10, 5);
h.push({3, 1});
auto [k, v] = h.top();
h.pop();
// Monotone requirement: each popped key must be >= the previously popped key.
```

### [`ds/heap/skew_heap.hpp`](ds/heap/skew_heap.hpp) — meldable, intrusive

```cpp
#include "ds/heap/skew_heap.hpp"

struct min_node : public skew_heap_node_base<min_node> {
	int32_t key;
	void do_downdate() {}
};
auto cmp = [](const min_node* a, const min_node* b) {
	return a->key <= b->key;           // min-heap
};

std::vector<min_node> pool(N);
min_node* root = nullptr;
for (int32_t i = 0; i < N; ++ i) {
	pool[i].key = a[i];
	root = merge(root, &pool[i], cmp);
}
min_node* mn = root;
root = pop(root, cmp);
```

---

## fft/

Everything lives inside `namespace fft`.

### [`fft/fft.hpp`](fft/fft.hpp)

```cpp
#include "fft/fft.hpp"

// Double-complex polynomial multiply (precision ~1e-8 for coefficients <= 1e6):
std::vector<int64_t> a = {1, 2, 3}, b = {4, 5};
auto c = fft::fft_complex_multiply<int64_t>(a, b);          // {4, 13, 22, 15}

// NTT multiply over Z_p (p must be an NTT prime):
using mint = modnum<constant<int32_t, 998244353>, naive_multiplier<int32_t>>;
std::vector<mint> p = {mint(1), mint(2)}, q = {mint(3), mint(4)};
auto r = fft::fft_numeric_multiply<mint>(p, q);

// Multiply modulo an arbitrary prime (triple-mod CRT):
auto s = fft::fft_complex_mod_multiply<int64_t>(a, b, 1000000007);

// Polynomial inverse modulo x^n: fft::*_inverse(...).
```

### [`fft/berlekamp_massey.hpp`](fft/berlekamp_massey.hpp)

```cpp
#include "fft/berlekamp_massey.hpp"

// Given a sequence A, return the shortest C such that
//   A[i] = C[0]*A[i-1] + C[1]*A[i-2] + ... + C[L-1]*A[i-L].
std::vector<mint> A = { /* ... */ };
auto C = berlekamp_massey(A);
```

### [`fft/kitamasa.hpp`](fft/kitamasa.hpp)

```cpp
#include "fft/kitamasa.hpp"

// k-th term of a linear recurrence in O(L^2 * log k).
std::vector<mint> A = {mint(1), mint(1)};                   // a[0]=1, a[1]=1
std::vector<mint> C = {mint(1), mint(1)};                   // a[i] = a[i-1] + a[i-2]
mint fib_k = kitamasa(A, C, /*k=*/int64_t(1'000'000'000));
```

---

## flow/

### [`flow/flow_graph.hpp`](flow/flow_graph.hpp)

```cpp
#include "flow/flow_graph.hpp"

flow_graph<int64_t> g(N);              // N vertices
g.add_edge(u, v, /*cap=*/10);          // back-edge with cap = 0 added automatically
g.add_edge(u, v, 10, /*back_cap=*/10); // undirected edge
// Each edge has an id; forward/back pair is (id, id ^ 1).
// For flow_t = double, set g.eps = 1e-9 before solving.
```

### [`flow/max_flow.hpp`](flow/max_flow.hpp)

```cpp
#include "flow/max_flow.hpp"

dinic_max_flow<int64_t> dn;
int64_t F = dn.max_flow(g, /*src=*/0, /*sink=*/N - 1);
auto cut  = dn.min_cut(g);             // vector<bool>: cut[i] = true iff i is on src side

hlpp_max_flow<int64_t> hp;
int64_t F2 = hp.max_flow(g, 0, N - 1);
```

### [`flow/dijkstra_min_cost_flow.hpp`](flow/dijkstra_min_cost_flow.hpp)

```cpp
#include "flow/dijkstra_min_cost_flow.hpp"

flow_graph<int64_t> g(N);
g.add_edge(u, v, cap, /*cost=*/c);     // cost passed as a trailing constructor arg

dijkstra_min_cost_flow<int64_t, int64_t> mcmf;
auto [flow, cost] = mcmf.min_cost_flow(g, src, snk);
// Requires non-negative initial costs (otherwise initialise potentials via SPFA).
```

### [`flow/network_simplex.hpp`](flow/network_simplex.hpp)

```cpp
#include "flow/network_simplex.hpp"

network_simplex<int64_t, int64_t> ns(N);
ns.add(u, v, /*lower=*/0, /*upper=*/cap, /*cost=*/c);
ns.add_supply(node, +5);
ns.add_supply(node, -5);
auto opt = ns.solve();                 // std::optional<cost>: nullopt if infeasible
```

---

## geo/

### [`geo/pair_point.hpp`](geo/pair_point.hpp)

```cpp
#include "geo/pair_point.hpp"

pair_point<int64_t> A(1, 2), B(3, 4);
auto    C   = A + B;                   // vector sum
int64_t dot = A * B;                   // dot product
int64_t crs = A ^ B;                   // cross product
int64_t crs3 = cross3(A, B, C);        // (B-A) x (C-A)
// Default scalar is int64_t; promote to __int128 when coordinates exceed ~2^32.
```

### [`geo/pair_line.hpp`](geo/pair_line.hpp)

```cpp
#include "geo/pair_line.hpp"

pair_line<int64_t> ln(A, B);           // segment / line through A, B
// Helpers: intersect, on_segment, signed_area2, ... — see the header.
```

### [`geo/convex_hull.hpp`](geo/convex_hull.hpp)

```cpp
#include "geo/convex_hull.hpp"

std::vector<pair_point<int64_t>> pts = { /* ... */ };
std::vector<int32_t> order(pts.size());
std::iota(order.begin(), order.end(), 0);
std::sort(order.begin(), order.end(),
          [&](int32_t a, int32_t b) { return pts[a] < pts[b]; });

auto res = convex_hull(pts, order);
// res[0]: indices of hull vertices in CCW order (no duplicated first vertex).
// res[1][i]: role of input point i (vertex / duplicate of another vertex /
//            colinear on an edge / strictly outside).
```

### [`geo/convex_inclusion.hpp`](geo/convex_inclusion.hpp)

```cpp
#include "geo/convex_inclusion.hpp"

// hull is already CCW; root_index is the index of the bottom-left vertex.
auto [inc, idx] = convex_inclusion(hull, query_point, root_index);
// inc = -1 inside, 0 on edge, +1 outside.
```

### [`geo/half_plane_intersect.hpp`](geo/half_plane_intersect.hpp)

```cpp
#include "geo/half_plane_intersect.hpp"

std::vector<pair_line<int64_t>> lines = { /* each line bounds a "left-side" half-plane */ };
auto poly = half_plane_intersect(lines);   // line indices in CCW order
bool unbounded = poly.size() < 3 ||
                 cross(lines[poly.front()].dir(), lines[poly.back()].dir()) >= 0;
```

### [`geo/minkowski_sum.hpp`](geo/minkowski_sum.hpp)

```cpp
#include "geo/minkowski_sum.hpp"

auto S = minkowski_sum(A, B, /*start_index_in_A=*/0);
// A and B are convex polygons in CCW order; S is the CCW Minkowski sum.
```

### [`geo/polygon.hpp`](geo/polygon.hpp)

```cpp
#include "geo/polygon.hpp"

int64_t s2 = polygon_area2(poly);      // signed area * 2
auto    cn = polygon_centroid(poly);   // pair_point<double> centroid
```

---

## graph/

### [`graph/graph.hpp`](graph/graph.hpp) — graph templates

```cpp
#include "graph/graph.hpp"

struct E {
	int32_t from, to;
	int64_t w;
	E(const int32_t& f, const int32_t& t, const int64_t& c) : from(f), to(t), w(c) {}
};

undigraph<E> g(N);
g.add_edge(u, v, /*w=*/int64_t(5));    // trailing args forwarded to E's ctor

// On an undirected graph, g(u, eid) returns "the other endpoint of edge eid
// when standing at u" — use this instead of a manual branch.
for (int32_t eid : g[u]) {
	int32_t v = g(u, eid);
	int64_t w = g.edges[eid].w;
}

digraph<E> dg(N);
condi_undigraph<E, std::function<bool(int32_t)>> cg(N, [&](int32_t eid){ return alive[eid]; });
```

### [`graph/ancestor.hpp`](graph/ancestor.hpp) — RMQ-based LCA

```cpp
#include "graph/ancestor.hpp"

rmq_ancestor anc;
anc.build(g, /*root=*/0);
int32_t l = anc.lca(u, v);
int32_t dist = anc.depth[u] + anc.depth[v] - 2 * anc.depth[l];
```

### [`graph/sparse_ancestor.hpp`](graph/sparse_ancestor.hpp) — binary lifting

```cpp
#include "graph/sparse_ancestor.hpp"

sparse_ancestor anc(g, /*root=*/0);
int32_t up = anc.kth_ancestor(u, k);
int32_t l  = anc.lca(u, v);
// Path aggregates can be layered on top of combine() — see the header.
```

### [`graph/build_bipartite.hpp`](graph/build_bipartite.hpp)

```cpp
#include "graph/build_bipartite.hpp"

std::vector<int32_t> side(N);
bool ok = build_bipartite(g, side);    // side[i] in {0, 1}; false if an odd cycle exists
```

### [`graph/centroid_tree.hpp`](graph/centroid_tree.hpp)

```cpp
#include "graph/centroid_tree.hpp"

centroid_tree ct;
ct.build(g);
int32_t p = ct.par[u];
// A → B in the original tree passes through at most O(log N) centroids;
// see the header doc-block for the canonical query pattern.
```

### [`graph/euler_tour.hpp`](graph/euler_tour.hpp)

```cpp
#include "graph/euler_tour.hpp"

euler_tour_edge et;
et.dfs(g, /*root=*/0);
// et.euler: sequence of visited edge ids (in/out)
// et.loc[u]: first index in `euler` when u is entered
```

### [`graph/eulerian.hpp`](graph/eulerian.hpp)

```cpp
#include "graph/eulerian.hpp"

// Returns the sequence of edge ids for an Eulerian path / circuit,
// or an empty vector if none exists.
auto seq = eulerian_circuit(g);
```

### [`graph/matching.hpp`](graph/matching.hpp) — bipartite matching

```cpp
#include "graph/matching.hpp"

digraph<simple_edge> g(N + M);         // only add left → right edges (left ids in [0, N))
for (auto [u, v] : pairs) g.add_edge(u, N + v);

dfs_matching mat(N, M);
int32_t sz = mat.max_match(g);
// mat.mate[0][u] = right vertex matched to left u, or -1.

hopcroft_karp_matching hk(N, M);
int32_t sz2 = hk.max_match(g);         // faster on dense edge sets
```

### [`graph/min_vertex_cover.hpp`](graph/min_vertex_cover.hpp)

> **Do not use**: the header has a known bug when `N != M` (loop bounds
> use the wrong matching side). Marked `[.skip]` in the test suite —
> roll your own König if you need it.

### [`graph/scc_comp.hpp`](graph/scc_comp.hpp) — Tarjan SCC

```cpp
#include "graph/scc_comp.hpp"

int32_t cnt = 0;
auto comp = scc_comp(g, cnt);          // comp[u] in [0, cnt)
// IMPORTANT: ids are emitted in REVERSE topological order — id 0 is a SINK.
```

### [`graph/span_forest.hpp`](graph/span_forest.hpp)

```cpp
#include "graph/span_forest.hpp"

dfs_span_forest sf;
sf.build(g, /*root=*/0);
// sf.par[u]: parent edge id; sf.back: list of non-tree edges.
```

### [`graph/span_tree.hpp`](graph/span_tree.hpp) — Kruskal MST

```cpp
#include "graph/span_tree.hpp"

int64_t sum = 0;
auto used = build_span_tree(g, sum);   // selected edge ids; sum = total weight
```

### [`graph/topo_sort.hpp`](graph/topo_sort.hpp)

```cpp
#include "graph/topo_sort.hpp"

auto order = topo_sort(g);             // empty if g has a cycle
```

### [`graph/tree_diameter.hpp`](graph/tree_diameter.hpp)

```cpp
#include "graph/tree_diameter.hpp"

tree_diameter td;
td.bfs(g);
int32_t u = td.x, v = td.y;
int64_t  D = td.dist[v];
```

### [`graph/twosat_graph.hpp`](graph/twosat_graph.hpp) — 2-SAT

```cpp
#include "graph/twosat_graph.hpp"

twosat_graph ts(N);                    // N variables
ts.add_clause(x, ~ y);                 // x OR NOT y  (use unary ~ for negation)
auto val = ts.sat();                   // vector<bool> of size N, or empty if UNSAT
```

---

## math/

### [`math/combinatorics.hpp`](math/combinatorics.hpp)

```cpp
#include "math/combinatorics.hpp"

combinatorics<mint> C(MAXN + 1);       // precompute fact / inv_fact over [0, MAXN]
mint c   = C.choose(n, k);
mint p   = C.permute(n, k);
mint cat = C.catalan(n);
mint f   = C.fact(n);
mint i   = C.inv_fact(n);
```

### [`math/derangement.hpp`](math/derangement.hpp)

```cpp
#include "math/derangement.hpp"
auto D = derangement<mint>(N);         // D[i] = !i
```

### [`math/diophantine.hpp`](math/diophantine.hpp)

```cpp
#include "math/diophantine.hpp"

int64_t g, x, y;
g = extended_gcd(a, b, x, y);          // a*x + b*y = g
bool ok = diophantine(a, b, c, x, y);  // a single solution to a*x + b*y = c
```

### [`math/div.hpp`](math/div.hpp)

```cpp
#include "math/div.hpp"

auto ds = divisors(n);                 // every divisor of n (NOT pre-sorted)
```

### [`math/euler_totient_function.h`](math/euler_totient_function.h)

```cpp
#include "math/euler_totient_function.h"   // note: extension is .h, not .hpp

int64_t p = phieuler(n);                  // O(sqrt n)
auto    phi = range_phieuler(N);          // phi[i] for i in [0, N]
```

### [`math/factorizer.hpp`](math/factorizer.hpp) — Pollard rho

```cpp
#include "math/factorizer.hpp"

factorizer::linear_sieve(1'000'000);                            // MUST call once before factorize
auto pfs = factorizer::factorize<uint64_t>(1234567890123ULL);
bool prime = factorizer::is_prime<uint64_t>(x);
```

### [`math/legendre_power.hpp`](math/legendre_power.hpp)

```cpp
#include "math/legendre_power.hpp"
int64_t k = legendre_power(n, p);      // largest k with p^k | n!
```

### [`math/matrix.hpp`](math/matrix.hpp) — fixed-size matrix

```cpp
#include "math/matrix.hpp"

matrix<mint, 3, 3> M;                  // dimensions are template parameters, NOT runtime
M(0, 0) = mint(1); /* ... */
auto P  = M * M;
auto I  = M.inverse();
auto Mk = pow(M, /*k=*/int64_t(1'000'000'000));   // fast exponentiation
```

### [`math/miller_rabin.hpp`](math/miller_rabin.hpp)

```cpp
#include "math/miller_rabin.hpp"
bool p = miller_rabin::is_prime<uint64_t>(n);     // deterministic for 64-bit integers
```

### [`math/sieve.hpp`](math/sieve.hpp)

```cpp
#include "math/sieve.hpp"

auto sv = eratos_sieve(N);             // vector<bool> is_prime[0..N]
auto pr = block_sieve(N);              // vector<int32_t> of primes <= N (cache-friendly, fastest)

linear_sieve ls(N);
auto&   primes = ls.primes;
int32_t spf    = ls.ls[x];             // smallest prime factor of x

// range_sieve in this file is BROKEN — do not use.
```

### [`math/sqrt.hpp`](math/sqrt.hpp)

```cpp
#include "math/sqrt.hpp"
int64_t r = isqrt(n);
auto    s = tonelli_shanks(a, p);      // modular square root (p an odd prime)
```

---

## misc/

### [`misc/compressor.hpp`](misc/compressor.hpp)

```cpp
#include "misc/compressor.hpp"

std::vector<int32_t> a = {30, 10, 20, 10};
auto sorted = a;                                                                   // the CALLER sorts
std::sort(sorted.begin(), sorted.end());
sorted.erase(std::unique(sorted.begin(), sorted.end()), sorted.end());
compressor<int32_t> cp(sorted);
int32_t id = cp(20);                   // -> 1
```

### [`misc/io.hpp`](misc/io.hpp)

```cpp
#include "misc/io.hpp"
// Fast reader / writer — see the header for exact prototypes.
```

### [`misc/reverse_args.hpp`](misc/reverse_args.hpp)

```cpp
#include "misc/reverse_args.hpp"

auto f = [](int32_t a, int32_t b, int32_t c) { return a - b - c; };
auto r = std::reverse_args_t(f);
r(1, 2, 3);                            // invokes f(3, 2, 1)
```

### [`misc/ska_sort.hpp`](misc/ska_sort.hpp)

```cpp
#include "misc/ska_sort.hpp"

std::vector<uint32_t> a = { /* ... */ };
ska_sort(a.begin(), a.end());

// Composite key via an extractor:
struct rec { int32_t k; std::string v; };
std::vector<rec> rs = { /* ... */ };
ska_sort(rs.begin(), rs.end(), [](const rec& r) { return r.k; });
```

### [`misc/tensor.hpp`](misc/tensor.hpp)

```cpp
#include "misc/tensor.hpp"

tensor<int64_t, 3> dp({N, M, K}, /*init=*/int64_t(0));
dp(i, j, k) = /* ... */;
```

### [`misc/y_combinator.hpp`](misc/y_combinator.hpp)

```cpp
#include "misc/y_combinator.hpp"

auto fib = std::y_combinator([](auto self, int32_t n) -> int32_t {
	return n < 2 ? n : self(n - 1) + self(n - 2);
});
int32_t x = fib(10);
```

---

## mod/

### [`mod/modnum.hpp`](mod/modnum.hpp)

```cpp
#include "mod/modnum.hpp"

// Compile-time prime modulus, default multiplier:
using mint  = modnum<constant<int32_t, 1000000007>, naive_multiplier<int32_t>>;
using mint2 = modnum<constant<int32_t,  998244353>, naive_multiplier<int32_t>>;

// Runtime prime modulus with Barrett reduction:
using dmint = modnum<inconstant<int32_t>, barrett_multiplier<int32_t>>;
inconstant<int32_t>::set_value(MOD);   // call once before any arithmetic

mint a(5), b(3);
mint c = a + b * a.pow(10);
int32_t raw = int32_t(c);              // there is NO .value() — use a cast
```

---

## string/

### [`string/prefix_function.hpp`](string/prefix_function.hpp)

```cpp
#include "string/prefix_function.hpp"

auto pi = prefix_function(s);          // s may be std::string, std::vector<int32_t>, ...
prefix_automaton aut(s);
int32_t nxt = aut.go(state, c);
```

### [`string/z_function.hpp`](string/z_function.hpp)

```cpp
#include "string/z_function.hpp"
auto z = z_function(s);                // z[0] = |s|
```

### [`string/suffix_array.hpp`](string/suffix_array.hpp)

```cpp
#include "string/suffix_array.hpp"
#include "ds/range_min_query.hpp"

suffix_array sa(s);                    // s must be COMPRESSED (alphabet in [0, sigma))
// sa[k]: starting position of the suffix that ranks k-th.
auto& rnk = sa.rank;                   // rnk[i] = rank of the suffix starting at i
auto& lcp = sa.lcp;                    // lcp[k] = LCP(sa[k], sa[k-1])

// LCP of arbitrary suffixes: combine sa.rank with an RMQ over sa.lcp.
range_min_query<int32_t> rmq(sa.lcp.begin() + 1, sa.lcp.end());
```

### [`string/manacher.hpp`](string/manacher.hpp)

```cpp
#include "string/manacher.hpp"

manacher m(s);
bool pal = m.is_palindrome(l, r);      // [l, r] INCLUSIVE
// Or iterate raw radii (length 2N - 1): for (int32_t r : m) { /* ... */ }
```

### [`string/prefix_tree.hpp`](string/prefix_tree.hpp) — trie

```cpp
#include "string/prefix_tree.hpp"

prefix_tree pt(26);                    // alphabet size
pt.make_root();
pt.insert(s);                          // s is viewed as a sequence in [0, 26)
bool has = pt.find(s);
pt.erase(s);
```

### [`string/prefix_power_hash.hpp`](string/prefix_power_hash.hpp)

```cpp
#include "string/prefix_power_hash.hpp"

using hnum  = modnum<constant<int32_t, 1'000'000'007>, naive_multiplier<int32_t>>;
using hpair = pairnum<hnum>;                                  // double-mod to fend off collisions

// Set the base ONCE before constructing the first hash:
base_power<hpair>::base = hpair(hnum(131), hnum(137));

prefix_power_hash<hpair> ph(s);
hpair h = ph.range_hash(l, r);         // hash of s[l..r)
```

---

## template/

The skeletons in [`template/`](template/) are starting points to drop
into a contest:

- [`temp.cpp`](template/temp.cpp): default Codeforces / generic
  contest starter.
- [`fhc.cpp`](template/fhc.cpp): Facebook Hacker Cup multi-test driver
  that writes to `output.txt`.
- [`leet.cpp`](template/leet.cpp): LeetCode local harness — class
  skeleton plus a `main` driver.
- [`Makefile`](template/Makefile): `make foo` builds `foo.cpp` with the
  project's stock flags. Mirror its CXXFLAGS in any external build glue.

These templates are **intentionally monolithic** — they `#include
<bits/stdc++.h>` and pull in helpers (debug, IO, `y_combinator`) the
author commonly uses in contests. Don't refactor them.

---

## General tips

- No header in the library includes `<bits/stdc++.h>`. The caller
  (test or contest file) is expected to include `<bits/stdc++.h>` (or
  a reasonable subset of the STL) before any library header.
- Linux-only headers (depend on `__gnu_pbds`):
  [`ds/hash_map.hpp`](ds/hash_map.hpp),
  [`ds/order_statistic.hpp`](ds/order_statistic.hpp),
  [`flow/dijkstra_min_cost_flow.hpp`](flow/dijkstra_min_cost_flow.hpp).
- Bound-check `assert`s are wrapped in `#ifdef _GLIBCXX_DEBUG`. Build
  with `-D_GLIBCXX_DEBUG` to enable them.
- If you hit `reverse_args_t not a member of std` while compiling
  `swag.hpp`, include `misc/reverse_args.hpp` *before* `ds/swag.hpp`.
- Benchmarks under [`tests/bench/`](tests/bench/) are hidden behind
  the Catch2 tag `[!benchmark]`. Run them explicitly with
  `./build/tests/bench/<module> "[!benchmark]"`.
- Known library bugs are catalogued in
  [`library_fixes.md`](library_fixes.md); current benchmark numbers
  live in [`benchmarks.md`](benchmarks.md).
- Pass inputs by `const T&` (the library does the same for primitives —
  e.g. `const int32_t&` is intentional). When in doubt, mirror the
  signature you see in the matching test under [`tests/`](tests/).
