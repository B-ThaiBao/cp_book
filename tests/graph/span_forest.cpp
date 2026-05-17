#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "graph/graph.hpp"
#include "graph/span_forest.hpp"

namespace {
struct edge_t { int from, to; };

std::vector<bool> naive_bridges(const undigraph<edge_t>& g) {
	int N = g.V, M = (int)g.edges.size();
	std::vector<bool> bri(M, false);
	auto count_comps = [&](const int& skip) {
		std::vector<int> par(N); std::iota(par.begin(), par.end(), 0);
		std::function<int(int)> find = [&](int x) { return par[x] == x ? x : par[x] = find(par[x]); };
		for (int i = 0; i < M; ++ i) {
			if (i == skip) continue;
			int a = find(g.edges[i].from), b = find(g.edges[i].to);
			if (a != b) par[a] = b;
		}
		int c = 0;
		for (int i = 0; i < N; ++ i) if (find(i) == i) ++ c;
		return c;
	};
	int base = count_comps(- 1);
	for (int i = 0; i < M; ++ i) if (count_comps(i) > base) bri[i] = true;
	return bri;
}

std::vector<bool> naive_cuts(const undigraph<edge_t>& g) {
	int N = g.V;
	std::vector<bool> cut(N, false);
	auto count_comps = [&](const int& skip) {
		std::vector<int> col(N, - 1); int c = 0;
		for (int s = 0; s < N; ++ s) {
			if (s == skip || col[s] != - 1) continue;
			col[s] = c;
			std::vector<int> q = {s};
			while (!q.empty()) {
				int u = q.back(); q.pop_back();
				for (const int& id : g.adj[u]) {
					int v = g(u, id);
					if (v == skip || col[v] != - 1) continue;
					col[v] = c; q.push_back(v);
				}
			}
			++ c;
		}
		return c;
	};
	int base = count_comps(- 1);
	for (int i = 0; i < N; ++ i) if (count_comps(i) > base) cut[i] = true;
	return cut;
}
} // namespace

TEST_CASE("dfs_span_forest: bridges and cutpoints stress", "[span_forest][stress]") {
	std::mt19937 rng(91);
	for (int t = 0; t < 40; ++ t) {
		int N = 3 + int(rng() % 6u);
		undigraph<edge_t> g(N);
		std::vector<int> perm(N); std::iota(perm.begin(), perm.end(), 0);
		std::shuffle(perm.begin(), perm.end(), rng);
		for (int i = 1; i < N; ++ i) g.add_edge(perm[i - 1], perm[i]);
		int M = int(rng() % 5u);
		std::set<std::pair<int, int>> seen;
		for (int i = 1; i < N; ++ i) {
			int a = std::min(perm[i - 1], perm[i]), b = std::max(perm[i - 1], perm[i]);
			seen.insert({a, b});
		}
		for (int i = 0; i < M; ++ i) {
			int a = int(rng() % unsigned(N)), b = int(rng() % unsigned(N));
			if (a == b) continue;
			if (a > b) std::swap(a, b);
			if (seen.insert({a, b}).second) g.add_edge(a, b);
		}
		dfs_span_forest dsf(N);
		dsf.dfs(g);
		REQUIRE(dsf.find_bridge(g) == naive_bridges(g));
		REQUIRE(dsf.find_cutpoint(g) == naive_cuts(g));
	}
}

TEST_CASE("bfs_span_forest: BFS levels on tree", "[span_forest][bfs]") {
	undigraph<edge_t> g(5);
	g.add_edge(0, 1);
	g.add_edge(0, 2);
	g.add_edge(1, 3);
	g.add_edge(2, 4);
	bfs_span_forest bsf(5);
	bsf.bfs(g, 0);
	REQUIRE(bsf.depth[0] == 0);
	REQUIRE(bsf.depth[1] == 1);
	REQUIRE(bsf.depth[2] == 1);
	REQUIRE(bsf.depth[3] == 2);
	REQUIRE(bsf.depth[4] == 2);
}

TEST_CASE("dfs_span_forest: tree has all bridges, no cuts at leaves", "[span_forest]") {
	int N = 5;
	undigraph<edge_t> g(N);
	for (int i = 0; i < N - 1; ++ i) g.add_edge(i, i + 1);
	dfs_span_forest dsf(N);
	dsf.dfs(g);
	auto bri = dsf.find_bridge(g);
	for (const bool& b : bri) REQUIRE(b);
	auto cut = dsf.find_cutpoint(g);
	REQUIRE_FALSE(cut[0]);
	REQUIRE_FALSE(cut[N - 1]);
	for (int i = 1; i < N - 1; ++ i) REQUIRE(cut[i]);
}

TEST_CASE("dfs_span_forest: cycle has no bridges or cuts", "[span_forest]") {
	int N = 4;
	undigraph<edge_t> g(N);
	for (int i = 0; i < N; ++ i) g.add_edge(i, (i + 1) % N);
	dfs_span_forest dsf(N);
	dsf.dfs(g);
	auto bri = dsf.find_bridge(g);
	for (const bool& b : bri) REQUIRE_FALSE(b);
	auto cut = dsf.find_cutpoint(g);
	for (const bool& c : cut) REQUIRE_FALSE(c);
}

TEST_CASE("bfs_span_forest: disconnected components", "[span_forest][bfs]") {
	undigraph<edge_t> g(4);
	g.add_edge(0, 1);
	g.add_edge(2, 3);
	bfs_span_forest bsf(4);
	bsf.bfs(g);
	REQUIRE(bsf.depth[0] == 0);
	REQUIRE(bsf.depth[1] == 1);
	REQUIRE(bsf.depth[2] == 0);
	REQUIRE(bsf.depth[3] == 1);
}
