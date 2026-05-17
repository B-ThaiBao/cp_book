#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "flow/flow_graph.hpp"
#include "flow/dijkstra_min_cost_flow.hpp"

TEST_CASE("dijkstra_min_cost_flow: small example", "[mcmf]") {
	cost_flow_graph<long long, long long> g(4);
	g.add_edge(0, 1, 2, 0, 1);
	g.add_edge(0, 2, 1, 0, 2);
	g.add_edge(1, 2, 1, 0, 1);
	g.add_edge(1, 3, 1, 0, 3);
	g.add_edge(2, 3, 2, 0, 1);
	dijkstra_min_cost_flow<long long, long long> mcmf(4);
	auto [f, c] = mcmf.min_cost_max_flow(g, 0, 3);
	REQUIRE(f == 3);
	REQUIRE(c == 10);
}

TEST_CASE("dijkstra_min_cost_flow: direct edge", "[mcmf]") {
	cost_flow_graph<long long, long long> g(2);
	g.add_edge(0, 1, 5, 0, 3);
	dijkstra_min_cost_flow<long long, long long> mcmf(2);
	auto [f, c] = mcmf.min_cost_max_flow(g, 0, 1);
	REQUIRE(f == 5);
	REQUIRE(c == 15);
}

TEST_CASE("dijkstra_min_cost_flow: disconnected", "[mcmf]") {
	cost_flow_graph<long long, long long> g(4);
	g.add_edge(0, 1, 5, 0, 1);
	g.add_edge(2, 3, 5, 0, 1);
	dijkstra_min_cost_flow<long long, long long> mcmf(4);
	auto [f, c] = mcmf.min_cost_max_flow(g, 0, 3);
	REQUIRE(f == 0);
	REQUIRE(c == 0);
}

TEST_CASE("dijkstra_min_cost_flow: parallel edges different costs", "[mcmf]") {
	cost_flow_graph<long long, long long> g(2);
	g.add_edge(0, 1, 3, 0, 5);
	g.add_edge(0, 1, 2, 0, 1);
	dijkstra_min_cost_flow<long long, long long> mcmf(2);
	auto [f, c] = mcmf.min_cost_max_flow(g, 0, 1);
	REQUIRE(f == 5);
	REQUIRE(c == 1 * 2 + 5 * 3);
}

TEST_CASE("dijkstra_min_cost_flow: conservation stress", "[mcmf][stress]") {
	std::mt19937 rng(99);
	for (int t = 0; t < 20; ++ t) {
		int N = 3 + int(rng() % 3u);
		cost_flow_graph<long long, long long> g(N);
		int M = 1 + int(rng() % 6u);
		for (int i = 0; i < M; ++ i) {
			int u = int(rng() % unsigned(N)), v = int(rng() % unsigned(N));
			if (u == v) continue;
			long long cap = 1 + (long long)(rng() % 4u);
			long long cost = 1 + (long long)(rng() % 5u);
			g.add_edge(u, v, cap, 0, cost);
		}
		dijkstra_min_cost_flow<long long, long long> mcmf(N);
		auto [f, c] = mcmf.min_cost_max_flow(g, 0, N - 1);
		REQUIRE(f >= 0);
		REQUIRE(c >= 0);
		std::vector<long long> bal(N, 0);
		for (size_t i = 0; i < g.edges.size(); i += 2) {
			bal[g.edges[i].from] += g.edges[i].flow;
			bal[g.edges[i].to] -= g.edges[i].flow;
		}
		REQUIRE(bal[0] == f);
		REQUIRE(bal[N - 1] == - f);
		for (int v = 1; v < N - 1; ++ v) REQUIRE(bal[v] == 0);
		long long check = 0;
		for (size_t i = 0; i < g.edges.size(); i += 2)
			check += g.edges[i].flow * g.edges[i].cost;
		REQUIRE(check == c);
	}
}
