#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "graph/graph.hpp"
#include "graph/scc_comp.hpp"

namespace {
struct edge_t { int from, to; };

std::vector<int> brute_scc(const int& N, const std::vector<std::pair<int, int>>& es) {
	std::vector<std::vector<bool>> reach(N, std::vector<bool>(N, false));
	for (int i = 0; i < N; ++ i) reach[i][i] = true;
	for (const auto& p : es) reach[p.first][p.second] = true;
	for (int k = 0; k < N; ++ k)
		for (int i = 0; i < N; ++ i)
			for (int j = 0; j < N; ++ j)
				if (reach[i][k] && reach[k][j]) reach[i][j] = true;
	std::vector<int> c(N, - 1);
	int cnt = 0;
	for (int i = 0; i < N; ++ i) {
		if (c[i] != - 1) continue;
		c[i] = cnt;
		for (int j = i + 1; j < N; ++ j)
			if (c[j] == - 1 && reach[i][j] && reach[j][i]) c[j] = cnt;
		++ cnt;
	}
	return c;
}
} // namespace

TEST_CASE("scc_comp: single SCC", "[scc]") {
	digraph<edge_t> g(3);
	g.add_edge(0, 1);
	g.add_edge(1, 2);
	g.add_edge(2, 0);
	int cnt;
	auto c = scc_comp(g, cnt);
	REQUIRE(cnt == 1);
	REQUIRE(c == std::vector<int>({0, 0, 0}));
}

TEST_CASE("scc_comp: chain in topo order", "[scc]") {
	digraph<edge_t> g(4);
	g.add_edge(0, 1);
	g.add_edge(1, 2);
	g.add_edge(2, 3);
	int cnt;
	auto c = scc_comp(g, cnt);
	REQUIRE(cnt == 4);
	for (const auto& e : g.edges) REQUIRE(c[e.from] <= c[e.to]);
}

TEST_CASE("scc_comp: two SCCs connected", "[scc]") {
	digraph<edge_t> g(5);
	g.add_edge(0, 1); g.add_edge(1, 2); g.add_edge(2, 0);
	g.add_edge(3, 4); g.add_edge(4, 3);
	g.add_edge(2, 3);
	int cnt;
	auto c = scc_comp(g, cnt);
	REQUIRE(cnt == 2);
	REQUIRE(c[0] == c[1]);
	REQUIRE(c[1] == c[2]);
	REQUIRE(c[3] == c[4]);
	REQUIRE(c[0] < c[3]);
}

TEST_CASE("scc_comp: isolated nodes", "[scc]") {
	digraph<edge_t> g(4);
	int cnt;
	auto c = scc_comp(g, cnt);
	REQUIRE(cnt == 4);
	std::set<int> s(c.begin(), c.end());
	REQUIRE(s.size() == 4);
}

TEST_CASE("scc_comp: stress vs brute", "[scc][stress]") {
	std::mt19937 rng(91);
	for (int t = 0; t < 40; ++ t) {
		int N = 2 + int(rng() % 6u);
		digraph<edge_t> g(N);
		std::vector<std::pair<int, int>> es;
		int M = int(rng() % 12u);
		for (int i = 0; i < M; ++ i) {
			int a = int(rng() % unsigned(N)), b = int(rng() % unsigned(N));
			g.add_edge(a, b);
			es.push_back({a, b});
		}
		int cnt;
		auto c = scc_comp(g, cnt);
		auto bc = brute_scc(N, es);
		// Same equivalence classes
		for (int i = 0; i < N; ++ i)
			for (int j = 0; j < N; ++ j)
				REQUIRE((c[i] == c[j]) == (bc[i] == bc[j]));
		// Topo property
		for (const auto& e : g.edges) REQUIRE(c[e.from] <= c[e.to]);
	}
}
