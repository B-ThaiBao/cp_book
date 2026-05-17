#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "ds/bst/link_cut_tree.hpp"

TEST_CASE("link_cut_tree: link/cut basic", "[lct]") {
	int N = 5;
	std::vector<link_cut_node> nodes(N);
	REQUIRE(link(&nodes[0], &nodes[1]));
	REQUIRE(link(&nodes[1], &nodes[2]));
	REQUIRE(link(&nodes[2], &nodes[3]));
	REQUIRE(same_comp(&nodes[0], &nodes[3]));
	REQUIRE_FALSE(same_comp(&nodes[0], &nodes[4]));
	REQUIRE(cut(&nodes[1], &nodes[2]));
	REQUIRE_FALSE(same_comp(&nodes[0], &nodes[3]));
	REQUIRE(same_comp(&nodes[2], &nodes[3]));
}

TEST_CASE("link_cut_tree: link prevents cycles", "[lct]") {
	int N = 3;
	std::vector<link_cut_node> nodes(N);
	REQUIRE(link(&nodes[0], &nodes[1]));
	REQUIRE(link(&nodes[1], &nodes[2]));
	REQUIRE_FALSE(link(&nodes[0], &nodes[2]));
}

TEST_CASE("link_cut_tree: find_lct_root + make_lct_root", "[lct]") {
	int N = 4;
	std::vector<link_cut_node> nodes(N);
	link(&nodes[0], &nodes[1]);
	link(&nodes[1], &nodes[2]);
	link(&nodes[2], &nodes[3]);
	make_lct_root(&nodes[0]);
	REQUIRE(find_lct_root(&nodes[3]) == &nodes[0]);
	make_lct_root(&nodes[3]);
	REQUIRE(find_lct_root(&nodes[0]) == &nodes[3]);
}

TEST_CASE("link_cut_tree: cut non-existent edge returns false", "[lct]") {
	std::vector<link_cut_node> nodes(3);
	link(&nodes[0], &nodes[1]);
	REQUIRE_FALSE(cut(&nodes[0], &nodes[2]));
}

TEST_CASE("link_cut_tree: stress vs DSU connectivity", "[lct][stress]") {
	std::mt19937 rng(2024);
	int N = 30;
	std::vector<link_cut_node> nodes(N);
	std::vector<int> par(N); std::iota(par.begin(), par.end(), 0);
	std::function<int(int)> rp = [&](int x) { return par[x] == x ? x : par[x] = rp(par[x]); };
	std::set<std::pair<int, int>> edges;
	for (int t = 0; t < 400; ++ t) {
		int u = int(rng() % unsigned(N));
		int v = int(rng() % unsigned(N));
		if (u == v) continue;
		auto e = std::minmax(u, v);
		int op = int(rng() % 3);
		if (op == 0) {
			bool naive_can = (rp(u) != rp(v));
			if (naive_can) {
				REQUIRE(link(&nodes[u], &nodes[v]));
				par[rp(u)] = rp(v);
				edges.insert({e.first, e.second});
			} else {
				REQUIRE_FALSE(link(&nodes[u], &nodes[v]));
			}
		} else if (op == 1) {
			if (edges.count({e.first, e.second})) {
				REQUIRE(cut(&nodes[u], &nodes[v]));
				edges.erase({e.first, e.second});
				std::iota(par.begin(), par.end(), 0);
				for (auto [a, b] : edges) par[rp(a)] = rp(b);
			}
		} else {
			REQUIRE(same_comp(&nodes[u], &nodes[v]) == (rp(u) == rp(v)));
		}
	}
}
