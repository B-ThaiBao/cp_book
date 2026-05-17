#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "graph/graph.hpp"
#include "graph/build_bipartite.hpp"

namespace {
struct edge_t { int from, to; };
} // namespace

TEST_CASE("build_bipartite: even cycle is bipartite", "[bipartite]") {
	undigraph<edge_t> g(4);
	g.add_edge(0, 1);
	g.add_edge(1, 2);
	g.add_edge(2, 3);
	g.add_edge(3, 0);
	auto side = build_bipartite(g);
	REQUIRE(side.size() == 4);
	for (const auto& e : g.edges) REQUIRE(side[e.from] != side[e.to]);
}

TEST_CASE("build_bipartite: odd cycle not bipartite", "[bipartite]") {
	undigraph<edge_t> g(3);
	g.add_edge(0, 1);
	g.add_edge(1, 2);
	g.add_edge(2, 0);
	auto side = build_bipartite(g);
	REQUIRE(side.empty());
}

TEST_CASE("build_bipartite: forest is bipartite", "[bipartite]") {
	undigraph<edge_t> g(6);
	g.add_edge(0, 1);
	g.add_edge(0, 2);
	g.add_edge(1, 3);
	g.add_edge(4, 5);
	auto side = build_bipartite(g);
	REQUIRE(side.size() == 6);
	for (const auto& e : g.edges) REQUIRE(side[e.from] != side[e.to]);
}

TEST_CASE("build_bipartite: single node", "[bipartite]") {
	undigraph<edge_t> g(1);
	auto side = build_bipartite(g);
	REQUIRE(side.size() == 1);
	REQUIRE(side[0] == 0);
}

TEST_CASE("build_bipartite: empty graph", "[bipartite]") {
	undigraph<edge_t> g(5);
	auto side = build_bipartite(g);
	REQUIRE(side.size() == 5);
	for (const auto& s : side) REQUIRE(s == 0);
}

TEST_CASE("build_bipartite: stress on random trees", "[bipartite][stress]") {
	std::mt19937 rng(133);
	for (int t = 0; t < 50; ++ t) {
		int N = 2 + int(rng() % 20u);
		undigraph<edge_t> g(N);
		for (int i = 1; i < N; ++ i) g.add_edge(int(rng() % unsigned(i)), i);
		auto side = build_bipartite(g);
		REQUIRE(side.size() == size_t(N));
		for (const auto& e : g.edges) REQUIRE(side[e.from] != side[e.to]);
	}
}
