#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "ds/disjoint_set.hpp"
#include "graph/graph.hpp"
#include "graph/span_tree.hpp"

namespace {
struct wedge_t { int from, to; long long cost; };
} // namespace

TEST_CASE("build_span_tree: classic example", "[mst]") {
	undigraph<wedge_t> g(4);
	g.add_edge(0, 1, 1LL);
	g.add_edge(1, 2, 2LL);
	g.add_edge(2, 3, 3LL);
	g.add_edge(0, 3, 100LL);
	long long sum = 0;
	auto ids = build_span_tree(g, sum);
	REQUIRE(ids.size() == 3);
	REQUIRE(sum == 6);
}

TEST_CASE("build_span_tree: disconnected forest", "[mst]") {
	undigraph<wedge_t> g(4);
	g.add_edge(0, 1, 5LL);
	g.add_edge(2, 3, 7LL);
	long long sum = 0;
	auto ids = build_span_tree(g, sum);
	REQUIRE(ids.size() == 2);
	REQUIRE(sum == 12);
}

TEST_CASE("build_span_tree: single node", "[mst]") {
	undigraph<wedge_t> g(1);
	long long sum = 0;
	auto ids = build_span_tree(g, sum);
	REQUIRE(ids.empty());
	REQUIRE(sum == 0);
}

TEST_CASE("build_span_tree: maximum spanning tree via greater", "[mst]") {
	undigraph<wedge_t> g(3);
	g.add_edge(0, 1, 1LL);
	g.add_edge(1, 2, 10LL);
	g.add_edge(0, 2, 5LL);
	long long sum = 0;
	auto ids = build_span_tree(g, sum, std::greater<>());
	REQUIRE(ids.size() == 2);
	REQUIRE(sum == 15);
}

TEST_CASE("build_span_tree: stress vs brute", "[mst][stress]") {
	std::mt19937 rng(7);
	for (int t = 0; t < 50; ++ t) {
		int N = 3 + int(rng() % 6u);
		undigraph<wedge_t> g(N);
		std::vector<int> perm(N); std::iota(perm.begin(), perm.end(), 0);
		std::shuffle(perm.begin(), perm.end(), rng);
		for (int i = 1; i < N; ++ i)
			g.add_edge(perm[i - 1], perm[i], (long long)(rng() % 100u));
		int M = int(rng() % 8u);
		for (int i = 0; i < M; ++ i) {
			int a = int(rng() % unsigned(N)), b = int(rng() % unsigned(N));
			if (a != b) g.add_edge(a, b, (long long)(rng() % 100u));
		}
		long long mst_sum = 0;
		auto ids = build_span_tree(g, mst_sum);
		long long ref_sum = 0; int ref_cnt = 0;
		std::vector<int> ord((int)g.edges.size()); std::iota(ord.begin(), ord.end(), 0);
		std::sort(ord.begin(), ord.end(), [&](int a, int b) { return g.edges[a].cost < g.edges[b].cost; });
		disjoint_set_size d(N);
		for (int id : ord) {
			const auto& e = g.edges[id];
			if (d.merge(e.from, e.to)) { ref_sum += e.cost; ++ ref_cnt; }
		}
		REQUIRE(ref_cnt == (int)ids.size());
		REQUIRE(ref_sum == mst_sum);
	}
}
