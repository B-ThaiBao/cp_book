#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "graph/graph.hpp"
#include "graph/euler_tour.hpp"

namespace {
struct edge_t { int from, to; };
} // namespace

TEST_CASE("euler_tour_edge: small tree", "[euler_tour]") {
	undigraph<edge_t> g(4);
	g.add_edge(0, 1);
	g.add_edge(0, 2);
	g.add_edge(1, 3);
	euler_tour_edge et(4);
	et.dfs(g, 0);
	for (int u = 0; u < 4; ++ u) REQUIRE(et.euler[et.loc[u]] == u);
	REQUIRE(et.euler.back() == - 1);
	REQUIRE((int)et.euler.size() == 2 * 4);
}

TEST_CASE("euler_tour_edge: single node", "[euler_tour]") {
	undigraph<edge_t> g(1);
	euler_tour_edge et(1);
	et.dfs(g, 0);
	REQUIRE(et.loc[0] == 0);
	REQUIRE(et.euler.size() == 2);
	REQUIRE(et.euler[0] == 0);
	REQUIRE(et.euler[1] == - 1);
}

TEST_CASE("euler_tour_edge: path", "[euler_tour]") {
	undigraph<edge_t> g(5);
	for (int i = 0; i < 4; ++ i) g.add_edge(i, i + 1);
	euler_tour_edge et(5);
	et.dfs(g, 0);
	REQUIRE(et.loc[0] == 0);
	REQUIRE(et.loc[4] == 4);
	REQUIRE((int)et.euler.size() == 2 * 5);
}

TEST_CASE("euler_tour_edge: star", "[euler_tour]") {
	undigraph<edge_t> g(5);
	for (int i = 1; i < 5; ++ i) g.add_edge(0, i);
	euler_tour_edge et(5);
	et.dfs(g, 0);
	REQUIRE(et.loc[0] == 0);
	for (int u = 0; u < 5; ++ u) REQUIRE(et.euler[et.loc[u]] == u);
}

TEST_CASE("euler_tour_edge: stress N matches 2N", "[euler_tour][stress]") {
	std::mt19937 rng(77);
	for (int t = 0; t < 30; ++ t) {
		int N = 2 + int(rng() % 30u);
		undigraph<edge_t> g(N);
		for (int i = 1; i < N; ++ i) g.add_edge(int(rng() % unsigned(i)), i);
		euler_tour_edge et(N);
		et.dfs(g, 0);
		REQUIRE((int)et.euler.size() == 2 * N);
		for (int u = 0; u < N; ++ u) REQUIRE(et.euler[et.loc[u]] == u);
	}
}
