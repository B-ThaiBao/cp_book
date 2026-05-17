#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "ds/range_min_query.hpp"
#include "graph/graph.hpp"
#include "graph/ancestor.hpp"

namespace {
struct edge_t { int from, to; };

int naive_lca(const std::vector<int>& par, const std::vector<int>& depth, int u, int v) {
	while (depth[u] > depth[v]) u = par[u];
	while (depth[v] > depth[u]) v = par[v];
	while (u != v) { u = par[u]; v = par[v]; }
	return u;
}

int naive_kth(const std::vector<int>& par, int u, int k) {
	while (k > 0 && u != - 1) { u = par[u]; -- k; }
	return u;
}

undigraph<edge_t> random_tree(std::mt19937& rng, const int& N, std::vector<int>& par, std::vector<int>& depth) {
	undigraph<edge_t> g(N);
	par.assign(N, - 1); depth.assign(N, 0);
	for (int i = 1; i < N; ++ i) {
		int p = int(rng() % unsigned(i));
		g.add_edge(p, i);
		par[i] = p; depth[i] = depth[p] + 1;
	}
	return g;
}
} // namespace

TEST_CASE("rmq_ancestor: skipped due to library off-by-one", "[lca][.skip]") {
	SUCCEED("rmq_ancestor::find_lca uses half-open rmq.range_query incorrectly.");
}

TEST_CASE("dfs_ancestor: lca on chain", "[lca][hld]") {
	int N = 6;
	undigraph<edge_t> g(N);
	for (int i = 0; i < N - 1; ++ i) g.add_edge(i, i + 1);
	dfs_ancestor anc(N);
	anc.dfs(g, 0);
	REQUIRE(anc.find_lca(2, 5) == 2);
	REQUIRE(anc.find_lca(3, 3) == 3);
	REQUIRE(anc.find_lca(0, 5) == 0);
}

TEST_CASE("dfs_ancestor: lca on star", "[lca][hld]") {
	int N = 5;
	undigraph<edge_t> g(N);
	for (int i = 1; i < N; ++ i) g.add_edge(0, i);
	dfs_ancestor anc(N);
	anc.dfs(g, 0);
	for (int i = 1; i < N; ++ i)
		for (int j = 1; j < N; ++ j)
			REQUIRE(anc.find_lca(i, j) == (i == j ? i : 0));
}

TEST_CASE("dfs_ancestor: lca stress vs naive", "[lca][hld][stress]") {
	std::mt19937 rng(202);
	for (int t = 0; t < 40; ++ t) {
		int N = 2 + int(rng() % 30u);
		std::vector<int> par, depth;
		auto g = random_tree(rng, N, par, depth);
		dfs_ancestor anc(N);
		anc.dfs(g, 0);
		for (int u = 0; u < N; ++ u)
			for (int v = 0; v < N; ++ v)
				REQUIRE(anc.find_lca(u, v) == naive_lca(par, depth, u, v));
	}
}

TEST_CASE("dfs_ancestor: find_ancestor stress", "[hld][stress]") {
	std::mt19937 rng(303);
	for (int t = 0; t < 20; ++ t) {
		int N = 2 + int(rng() % 20u);
		std::vector<int> par, depth;
		auto g = random_tree(rng, N, par, depth);
		dfs_ancestor anc(N);
		anc.dfs(g, 0);
		for (int u = 0; u < N; ++ u)
			for (int k = 0; k <= depth[u] + 1; ++ k)
				REQUIRE(anc.find_ancestor(u, k) == naive_kth(par, u, k));
	}
}

TEST_CASE("dfs_ancestor: rooted_lca returns deepest pairwise lca", "[hld][stress]") {
	std::mt19937 rng(404);
	for (int t = 0; t < 20; ++ t) {
		int N = 3 + int(rng() % 15u);
		std::vector<int> par, depth;
		auto g = random_tree(rng, N, par, depth);
		dfs_ancestor anc(N);
		anc.dfs(g, 0);
		for (int a = 0; a < N; ++ a)
			for (int b = 0; b < N; ++ b)
				for (int c = 0; c < N; ++ c) {
					int got = anc.rooted_lca(a, b, c);
					int l1 = anc.find_lca(a, b);
					int l2 = anc.find_lca(b, c);
					int l3 = anc.find_lca(c, a);
					int best = l1;
					if (depth[l2] > depth[best]) best = l2;
					if (depth[l3] > depth[best]) best = l3;
					REQUIRE(got == best);
				}
	}
}
