// Benchmarks for graph/* (SCC, MST, topo, bipartite matching).
#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>

#include "ds/disjoint_set.hpp"
#include "graph/graph.hpp"
#include "graph/scc_comp.hpp"
#include "graph/span_tree.hpp"
#include "graph/topo_sort.hpp"
#include "graph/matching.hpp"

namespace {
struct E { int from, to; };
struct WE { int from, to; long long cost; };
}

TEST_CASE("scc_comp: N=1e5, M=5e5 (random digraph)", "[!benchmark][scc]") {
	std::mt19937_64 rng(911);
	digraph<E> g(100'000);
	for (int i = 0; i < 500'000; ++i) {
		g.add_edge(int(rng() % 100'000ULL), int(rng() % 100'000ULL));
	}
	BENCHMARK("Tarjan SCC (1e5 nodes / 5e5 edges)") {
		int cnt = 0;
		auto comp = scc_comp(g, cnt);
		return cnt;
	};
}

TEST_CASE("topo_sort: random DAG (N=1e5 / M~5e5)", "[!benchmark][topo]") {
	std::mt19937_64 rng(913);
	digraph<E> g(100'000);
	for (int i = 0; i < 500'000; ++i) {
		int u = int(rng() % 100'000ULL);
		int v = int(rng() % 100'000ULL);
		if (u == v) continue;
		if (u > v) std::swap(u, v);
		g.add_edge(u, v);
	}
	BENCHMARK("topo_sort (DAG, ~5e5 edges)") {
		auto order = topo_sort(g);
		return order.size();
	};
}

TEST_CASE("build_span_tree: Kruskal on N=1e5 / M=5e5", "[!benchmark][mst]") {
	std::mt19937_64 rng(917);
	undigraph<WE> g(100'000);
	for (int i = 0; i < 500'000; ++i) {
		int u = int(rng() % 100'000ULL);
		int v = int(rng() % 100'000ULL);
		g.add_edge(u, v, (long long)(rng() & 0xffff));
	}
	BENCHMARK("Kruskal MST (1e5 / 5e5)") {
		long long sum = 0;
		auto res = build_span_tree(g, sum);
		return res.size();
	};
}

TEST_CASE("matching: dfs vs hopcroft_karp (N=M=1e4, |E|=8e4)", "[!benchmark][matching]") {
	constexpr int N = 10'000;
	constexpr int M = 80'000;
	std::mt19937_64 rng(919);
	std::vector<std::pair<int, int>> edges(M);
	for (auto& [u, v] : edges) {
		u = int(rng() % unsigned(N));
		v = int(rng() % unsigned(N));
	}

	BENCHMARK("dfs_matching (1e4 + 1e4, 8e4 edges)") {
		digraph<E> g(N);
		for (auto& [u, v] : edges) g.add_edge(u, v);
		dfs_matching mat(N, N);
		return mat.max_match(g);
	};

	BENCHMARK("hopcroft_karp_matching (1e4 + 1e4, 8e4 edges)") {
		digraph<E> g(N);
		for (auto& [u, v] : edges) g.add_edge(u, v);
		hopcroft_karp_matching mat(N, N);
		return mat.max_match(g);
	};
}
