#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "graph/graph.hpp"
#include "graph/topo_sort.hpp"

namespace {
struct edge_t { int from, to; };

bool is_valid_topo(const digraph<edge_t>& g, const std::vector<int>& order) {
	if ((int)order.size() != g.V) return false;
	std::vector<int> pos(g.V, - 1);
	for (int i = 0; i < (int)order.size(); ++ i) pos[order[i]] = i;
	for (const auto& e : g.edges)
		if (pos[e.from] > pos[e.to]) return false;
	return true;
}
} // namespace

TEST_CASE("topo_sort: simple DAG", "[topo]") {
	digraph<edge_t> g(4);
	g.add_edge(0, 1);
	g.add_edge(0, 2);
	g.add_edge(1, 3);
	g.add_edge(2, 3);
	auto order = topo_sort(g);
	REQUIRE(is_valid_topo(g, order));
	REQUIRE(order[0] == 0);
	REQUIRE(order.back() == 3);
}

TEST_CASE("topo_sort: empty graph", "[topo]") {
	digraph<edge_t> g(5);
	auto order = topo_sort(g);
	REQUIRE((int)order.size() == 5);
}

TEST_CASE("topo_sort: cycle returns empty", "[topo]") {
	digraph<edge_t> g(3);
	g.add_edge(0, 1);
	g.add_edge(1, 2);
	g.add_edge(2, 0);
	auto order = topo_sort(g);
	REQUIRE(order.empty());
}

TEST_CASE("topo_sort: self-loop cycle", "[topo]") {
	digraph<edge_t> g(2);
	g.add_edge(0, 0);
	auto order = topo_sort(g);
	REQUIRE(order.empty());
}

TEST_CASE("topo_sort: chain", "[topo]") {
	digraph<edge_t> g(5);
	for (int i = 0; i < 4; ++ i) g.add_edge(i, i + 1);
	auto order = topo_sort(g);
	REQUIRE(order == std::vector<int>({0, 1, 2, 3, 4}));
}

TEST_CASE("topo_sort: stress random DAG", "[topo][stress]") {
	std::mt19937 rng(42);
	for (int t = 0; t < 60; ++ t) {
		int N = 2 + int(rng() % 10u);
		digraph<edge_t> g(N);
		std::vector<int> perm(N);
		std::iota(perm.begin(), perm.end(), 0);
		std::shuffle(perm.begin(), perm.end(), rng);
		int M = int(rng() % 15u);
		for (int i = 0; i < M; ++ i) {
			int a = int(rng() % unsigned(N)), b = int(rng() % unsigned(N));
			if (a == b) continue;
			if (a > b) std::swap(a, b);
			g.add_edge(perm[a], perm[b]);
		}
		auto order = topo_sort(g);
		REQUIRE(is_valid_topo(g, order));
	}
}

TEST_CASE("topo_sort: stress detect cycles", "[topo][stress]") {
	std::mt19937 rng(43);
	for (int t = 0; t < 30; ++ t) {
		int N = 3 + int(rng() % 6u);
		digraph<edge_t> g(N);
		for (int i = 0; i < N; ++ i) g.add_edge(i, (i + 1) % N);
		auto order = topo_sort(g);
		REQUIRE(order.empty());
	}
}
