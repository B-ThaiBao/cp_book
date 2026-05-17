#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "ds/disjoint_set.hpp"

TEST_CASE("disjoint_set_size: basic", "[dsu]") {
	disjoint_set_size dsu(6);
	REQUIRE(dsu.size() == 6);
	REQUIRE(dsu.size(0) == 1);
	REQUIRE(dsu.merge(0, 1));
	REQUIRE(dsu.merge(2, 3));
	REQUIRE(dsu.merge(0, 2));
	REQUIRE_FALSE(dsu.merge(1, 3));
	REQUIRE(dsu.size() == 3);
	REQUIRE(dsu.find_par(0) == dsu.find_par(3));
	REQUIRE(dsu.find_par(4) != dsu.find_par(5));
	REQUIRE(dsu.size(1) == 4);
	REQUIRE(dsu.size(4) == 1);
}

TEST_CASE("disjoint_set_size: all singletons", "[dsu]") {
	disjoint_set_size dsu(5);
	REQUIRE(dsu.size() == 5);
	for (int i = 0; i < 5; ++ i) REQUIRE(dsu.size(i) == 1);
}

TEST_CASE("disjoint_set_size: chain merges into one", "[dsu]") {
	disjoint_set_size dsu(10);
	for (int i = 0; i < 9; ++ i) dsu.merge(i, i + 1);
	REQUIRE(dsu.size() == 1);
	REQUIRE(dsu.size(0) == 10);
}

TEST_CASE("disjoint_set_size: stress vs naive", "[dsu][stress]") {
	std::mt19937 rng(101);
	int N = 30;
	disjoint_set_size dsu(N);
	std::vector<int> par(N); std::iota(par.begin(), par.end(), 0);
	std::function<int(int)> rp = [&](int x) { return par[x] == x ? x : par[x] = rp(par[x]); };
	for (int t = 0; t < 1000; ++ t) {
		int op = rng() % 2;
		int a = rng() % N, b = rng() % N;
		if (op == 0) {
			bool naive_diff = (rp(a) != rp(b));
			bool got = dsu.merge(a, b);
			if (naive_diff) par[rp(a)] = rp(b);
			REQUIRE(got == naive_diff);
		} else {
			REQUIRE((dsu.find_par(a) == dsu.find_par(b)) == (rp(a) == rp(b)));
		}
	}
}

TEST_CASE("bipartite_graph: tree is bipartite", "[dsu][bipartite]") {
	bipartite_graph g(5);
	g.add_edge(0, 1);
	g.add_edge(1, 2);
	g.add_edge(2, 3);
	REQUIRE(g.is_bipartite(0));
	REQUIRE(g.is_bipartite(3));
}

TEST_CASE("bipartite_graph: odd cycle not bipartite", "[dsu][bipartite]") {
	bipartite_graph h(3);
	h.add_edge(0, 1);
	h.add_edge(1, 2);
	h.add_edge(2, 0);
	REQUIRE_FALSE(h.is_bipartite(0));
}

TEST_CASE("bipartite_graph: is_same query", "[dsu][bipartite]") {
	bipartite_graph g(4);
	g.add_edge(0, 1);
	g.add_edge(1, 2);
	REQUIRE(g.is_same(0, 1));
	REQUIRE(g.is_same(0, 2));
	REQUIRE_FALSE(g.is_same(0, 3));
}

TEST_CASE("disjoint_set_ancestor: basic", "[dsu][ancestor]") {
	disjoint_set_ancestor dsa(5);
	REQUIRE(dsa.size() == 5);
	REQUIRE(dsa.merge(0, 1));
	REQUIRE(dsa.merge(1, 2));
	REQUIRE_FALSE(dsa.merge(0, 2));
	REQUIRE(dsa.size() == 3);
	REQUIRE(dsa.find_par(0) == dsa.find_par(2));
	REQUIRE(dsa.find_par(3) != dsa.find_par(0));
}
