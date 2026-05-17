#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "graph/sparse_ancestor.hpp"

namespace {
int naive_kth(const std::vector<int>& par, int u, int k) {
	while (k > 0 && u != - 1) { u = par[u]; -- k; }
	return u;
}

int naive_lca(const std::vector<int>& par, const std::vector<int>& depth, int u, int v) {
	while (depth[u] > depth[v]) u = par[u];
	while (depth[v] > depth[u]) v = par[v];
	while (u != v) { u = par[u]; v = par[v]; }
	return u;
}
} // namespace

TEST_CASE("sparse_ancestor: kth ancestor stress", "[sparse_ancestor][stress]") {
	std::mt19937 rng(11);
	for (int t = 0; t < 40; ++ t) {
		int N = 2 + int(rng() % 30u);
		std::vector<int> par(N, - 1), depth(N, 0);
		for (int i = 1; i < N; ++ i) {
			par[i] = int(rng() % unsigned(i));
			depth[i] = depth[par[i]] + 1;
		}
		sparse_ancestor sa(par);
		for (int u = 0; u < N; ++ u)
			for (int k = 1; k <= depth[u] + 1; ++ k)
				REQUIRE(sa.find_ancestor(u, k) == naive_kth(par, u, k));
	}
}

TEST_CASE("sparse_ancestor: lca stress", "[sparse_ancestor][lca][stress]") {
	std::mt19937 rng(22);
	for (int t = 0; t < 40; ++ t) {
		int N = 2 + int(rng() % 30u);
		std::vector<int> par(N, - 1), depth(N, 0);
		for (int i = 1; i < N; ++ i) {
			par[i] = int(rng() % unsigned(i));
			depth[i] = depth[par[i]] + 1;
		}
		sparse_ancestor sa(par);
		for (int u = 0; u < N; ++ u)
			for (int v = 0; v < N; ++ v)
				REQUIRE(sa.find_lca(u, v, depth) == naive_lca(par, depth, u, v));
	}
}

TEST_CASE("sparse_ancestor: chain", "[sparse_ancestor]") {
	int N = 10;
	std::vector<int> par(N, - 1), depth(N, 0);
	for (int i = 1; i < N; ++ i) { par[i] = i - 1; depth[i] = i; }
	sparse_ancestor sa(par);
	REQUIRE(sa.find_ancestor(9, 1) == 8);
	REQUIRE(sa.find_ancestor(9, 5) == 4);
	REQUIRE(sa.find_ancestor(9, 9) == 0);
	REQUIRE(sa.find_lca(3, 7, depth) == 3);
	REQUIRE(sa.find_lca(5, 5, depth) == 5);
}

TEST_CASE("sparse_ancestor: star", "[sparse_ancestor]") {
	int N = 6;
	std::vector<int> par(N, - 1), depth(N, 0);
	for (int i = 1; i < N; ++ i) { par[i] = 0; depth[i] = 1; }
	sparse_ancestor sa(par);
	for (int i = 1; i < N; ++ i) REQUIRE(sa.find_ancestor(i, 1) == 0);
	for (int i = 1; i < N; ++ i)
		for (int j = 1; j < N; ++ j) REQUIRE(sa.find_lca(i, j, depth) == (i == j ? i : 0));
}
