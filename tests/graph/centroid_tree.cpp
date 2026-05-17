#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "graph/graph.hpp"
#include "graph/centroid_tree.hpp"

namespace {
struct edge_t { int from, to; };

int tree_height(const std::vector<int>& par, const int& N) {
	std::vector<std::vector<int>> ch(N);
	int root = - 1;
	for (int i = 0; i < N; ++ i) {
		if (par[i] == - 1) root = i;
		else ch[par[i]].push_back(i);
	}
	int h = 0;
	std::function<void(int, int)> dfs = [&](int u, int d) {
		h = std::max(h, d);
		for (int v : ch[u]) dfs(v, d + 1);
	};
	dfs(root, 0);
	return h;
}
} // namespace

TEST_CASE("centroid_tree: chain of 5", "[centroid]") {
	undigraph<edge_t> g(5);
	for (int i = 0; i < 4; ++ i) g.add_edge(i, i + 1);
	centroid_tree ct(5);
	ct.build_tree(g, 0);
	int root = - 1;
	for (int i = 0; i < 5; ++ i) if (ct.par[i] == - 1) root = i;
	REQUIRE(root == 2);
}

TEST_CASE("centroid_tree: single node", "[centroid]") {
	undigraph<edge_t> g(1);
	centroid_tree ct(1);
	ct.build_tree(g, 0);
	REQUIRE(ct.par[0] == - 1);
}

TEST_CASE("centroid_tree: star root is center", "[centroid]") {
	int N = 7;
	undigraph<edge_t> g(N);
	for (int i = 1; i < N; ++ i) g.add_edge(0, i);
	centroid_tree ct(N);
	ct.build_tree(g, 0);
	int root = - 1;
	for (int i = 0; i < N; ++ i) if (ct.par[i] == - 1) root = i;
	REQUIRE(root == 0);
}

TEST_CASE("centroid_tree: log height stress", "[centroid][stress]") {
	std::mt19937 rng(31);
	for (int t = 0; t < 30; ++ t) {
		int N = 1 + int(rng() % 50u);
		undigraph<edge_t> g(N);
		for (int i = 1; i < N; ++ i) g.add_edge(int(rng() % unsigned(i)), i);
		centroid_tree ct(N);
		ct.build_tree(g, 0);
		int roots = 0;
		for (int i = 0; i < N; ++ i) if (ct.par[i] == - 1) ++ roots;
		REQUIRE(roots == 1);
		int h = tree_height(ct.par, N);
		int bound = 0; while ((1 << bound) < N) ++ bound; bound += 1;
		REQUIRE(h <= bound);
	}
}

TEST_CASE("centroid_tree: all nodes reachable via par", "[centroid][stress]") {
	std::mt19937 rng(32);
	for (int t = 0; t < 20; ++ t) {
		int N = 2 + int(rng() % 30u);
		undigraph<edge_t> g(N);
		for (int i = 1; i < N; ++ i) g.add_edge(int(rng() % unsigned(i)), i);
		centroid_tree ct(N);
		ct.build_tree(g, 0);
		// every node eventually ends at root
		for (int i = 0; i < N; ++ i) {
			int u = i, steps = 0;
			while (ct.par[u] != - 1) { u = ct.par[u]; if (++ steps > N) break; }
			REQUIRE(steps <= N);
		}
	}
}
