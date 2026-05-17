#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "flow/flow_graph.hpp"

TEST_CASE("flow_graph: add_edge creates forward + reverse", "[flow_graph]") {
	flow_graph<int> g(3);
	int e = g.add_edge(0, 1, 5);
	REQUIRE(e == 0);
	REQUIRE(g.edges.size() == 2);
	REQUIRE(g.edges[0].from == 0);
	REQUIRE(g.edges[0].to == 1);
	REQUIRE(g.edges[0].cap == 5);
	REQUIRE(g.edges[1].from == 1);
	REQUIRE(g.edges[1].to == 0);
	REQUIRE(g.edges[1].cap == 0);
	REQUIRE(g.adj[0] == std::vector<int>({0}));
	REQUIRE(g.adj[1] == std::vector<int>({1}));
}

TEST_CASE("flow_graph: multiple edges", "[flow_graph]") {
	flow_graph<int> g(4);
	g.add_edge(0, 1, 5);
	g.add_edge(1, 2, 3);
	g.add_edge(2, 3, 7);
	REQUIRE(g.edges.size() == 6);
	REQUIRE(g.adj[0].size() == 1);
	REQUIRE(g.adj[1].size() == 2);
	REQUIRE(g.adj[2].size() == 2);
	REQUIRE(g.adj[3].size() == 1);
}

TEST_CASE("flow_graph: clear_flow", "[flow_graph]") {
	flow_graph<int> g(2);
	g.add_edge(0, 1, 5);
	g.edges[0].flow = 3;
	g.edges[1].flow = - 3;
	g.clear_flow();
	REQUIRE(g.edges[0].flow == 0);
	REQUIRE(g.edges[1].flow == 0);
}

TEST_CASE("cost_flow_graph: forward + reverse with opposite cost", "[flow_graph]") {
	cost_flow_graph<int, int> g(2);
	g.add_edge(0, 1, 5, 0, 7);
	REQUIRE(g.edges.size() == 2);
	REQUIRE(g.edges[0].cost == 7);
	REQUIRE(g.edges[1].cost == - 7);
}

TEST_CASE("cost_flow_graph: lower bound stored", "[flow_graph]") {
	cost_flow_graph<int, int> g(2);
	g.add_edge(0, 1, 10, 3, 5);
	REQUIRE(g.edges.size() == 2);
	REQUIRE(g.edges[0].cap == 10);
}
