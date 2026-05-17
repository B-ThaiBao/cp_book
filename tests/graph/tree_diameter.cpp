#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "graph/graph.hpp"
#include "graph/tree_diameter.hpp"

namespace {
struct wedge_t { int from, to; long long cost; };
} // namespace

TEST_CASE("tree_diameter: path graph", "[diameter]") {
	undigraph<wedge_t> g(5);
	for (int i = 0; i < 4; ++ i) g.add_edge(i, i + 1, 1LL);
	tree_diameter<long long> td(5);
	auto [x, y] = td.bfs(g, 0);
	REQUIRE(td.dist[y] == 4);
	REQUIRE(((x == 0 && y == 4) || (x == 4 && y == 0)));
}

TEST_CASE("tree_diameter: star", "[diameter]") {
	undigraph<wedge_t> g(5);
	for (int i = 1; i < 5; ++ i) g.add_edge(0, i, 3LL);
	tree_diameter<long long> td(5);
	auto [x, y] = td.bfs(g, 0);
	REQUIRE(td.dist[y] == 6);
	REQUIRE(x != 0);
	REQUIRE(y != 0);
}

TEST_CASE("tree_diameter: weighted skewed", "[diameter]") {
	undigraph<wedge_t> g(4);
	g.add_edge(0, 1, 10LL);
	g.add_edge(1, 2, 1LL);
	g.add_edge(1, 3, 100LL);
	tree_diameter<long long> td(4);
	auto [x, y] = td.bfs(g, 0);
	REQUIRE(td.dist[y] == 110);
}

TEST_CASE("tree_diameter: single node", "[diameter]") {
	undigraph<wedge_t> g(1);
	tree_diameter<long long> td(1);
	auto [x, y] = td.bfs(g, 0);
	REQUIRE(x == 0);
	REQUIRE(y == 0);
	REQUIRE(td.dist[0] == 0);
}

TEST_CASE("tree_diameter: stress against double BFS naive", "[diameter][stress]") {
	std::mt19937 rng(311);
	for (int t = 0; t < 30; ++ t) {
		int N = 2 + int(rng() % 20u);
		undigraph<wedge_t> g(N);
		std::vector<std::vector<std::pair<int, long long>>> adj(N);
		for (int i = 1; i < N; ++ i) {
			int p = int(rng() % unsigned(i));
			long long c = 1 + (long long)(rng() % 20u);
			g.add_edge(p, i, c);
			adj[p].push_back({i, c});
			adj[i].push_back({p, c});
		}
		tree_diameter<long long> td(N);
		auto [x, y] = td.bfs(g, 0);
		// naive: all-pairs BFS, take max distance
		long long best = 0;
		for (int s = 0; s < N; ++ s) {
			std::vector<long long> d(N, - 1);
			d[s] = 0;
			std::queue<int> q; q.push(s);
			while (!q.empty()) {
				int u = q.front(); q.pop();
				for (const auto& [v, c] : adj[u]) {
					if (d[v] == - 1) { d[v] = d[u] + c; q.push(v); }
				}
			}
			for (const auto& x_ : d) best = std::max(best, x_);
		}
		REQUIRE(td.dist[y] == best);
	}
}
