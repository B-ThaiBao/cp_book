#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "flow/network_simplex.hpp"

TEST_CASE("network_simplex: transportation problem", "[simplex]") {
	network_simplex<long long, long long> ns(4);
	ns.nodes[0].supply = 5;
	ns.nodes[1].supply = 3;
	ns.nodes[2].supply = - 4;
	ns.nodes[3].supply = - 4;
	ns.add_edge(0, 2, 0, 10, 1);
	ns.add_edge(0, 3, 0, 10, 2);
	ns.add_edge(1, 2, 0, 10, 2);
	ns.add_edge(1, 3, 0, 10, 1);
	ns.simplex();
	REQUIRE(ns.min_cost());
	long long total = 0;
	for (int i = 0; i < 4; ++ i) total += ns.edges[i].flow * ns.edges[i].cost;
	REQUIRE(total == 9);
}

TEST_CASE("network_simplex: infeasible (no path)", "[simplex]") {
	network_simplex<long long, long long> ns(2);
	ns.nodes[0].supply = 5;
	ns.nodes[1].supply = - 5;
	ns.simplex();
	REQUIRE_FALSE(ns.min_cost());
}

TEST_CASE("network_simplex: zero supply trivial", "[simplex]") {
	network_simplex<long long, long long> ns(3);
	ns.add_edge(0, 1, 0, 10, 1);
	ns.add_edge(1, 2, 0, 10, 1);
	ns.simplex();
	REQUIRE(ns.min_cost());
}

TEST_CASE("network_simplex: min cost max flow via INF supply", "[simplex][mcmf]") {
	constexpr long long INF = 1LL << 30;
	network_simplex<long long, long long> ns(4);
	ns.nodes[0].supply = INF;
	ns.nodes[3].supply = - INF;
	ns.add_edge(0, 1, 0, 3, 1);
	ns.add_edge(0, 2, 0, 2, 2);
	ns.add_edge(1, 3, 0, 2, 1);
	ns.add_edge(2, 3, 0, 3, 2);
	ns.simplex();
	long long mf = ns.min_cost_flow();
	REQUIRE(mf == 4);
}

TEST_CASE("network_simplex: single direct edge feasibility", "[simplex]") {
	network_simplex<long long, long long> ns(2);
	ns.nodes[0].supply = 3;
	ns.nodes[1].supply = - 3;
	ns.add_edge(0, 1, 0, 10, 7);
	ns.simplex();
	REQUIRE(ns.min_cost());
	REQUIRE(ns.edges[0].flow == 3);
}

TEST_CASE("network_simplex: two-path optimal", "[simplex]") {
	network_simplex<long long, long long> ns(3);
	ns.nodes[0].supply = 4;
	ns.nodes[2].supply = - 4;
	ns.add_edge(0, 2, 0, 2, 5);
	ns.add_edge(0, 1, 0, 3, 1);
	ns.add_edge(1, 2, 0, 3, 1);
	ns.simplex();
	REQUIRE(ns.min_cost());
	long long total = 0;
	for (int i = 0; i < 3; ++ i) total += ns.edges[i].flow * ns.edges[i].cost;
	REQUIRE(total == 1 * 3 * 2 + 1 * 5);
}
