#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "graph/graph.hpp"

namespace {
struct edge_t { int from, to; };
} // namespace

TEST_CASE("undigraph: add_edge and adjacency", "[graph]") {
	undigraph<edge_t> g(4);
	int e0 = g.add_edge(0, 1);
	int e1 = g.add_edge(1, 2);
	int e2 = g.add_edge(2, 0);
	REQUIRE(e0 == 0);
	REQUIRE(e1 == 1);
	REQUIRE(e2 == 2);
	REQUIRE(g.adj[0].size() == 2);
	REQUIRE(g.adj[1].size() == 2);
	REQUIRE(g.adj[2].size() == 2);
	REQUIRE(g.adj[3].empty());
	REQUIRE(g(0, e0) == 1);
	REQUIRE(g(1, e0) == 0);
	REQUIRE(g(2, e1) == 1);
}

TEST_CASE("undigraph: parallel edges & self loop", "[graph]") {
	undigraph<edge_t> g(2);
	g.add_edge(0, 1);
	g.add_edge(0, 1);
	g.add_edge(0, 0);
	REQUIRE(g.edges.size() == 3);
	REQUIRE(g.adj[0].size() == 4);
	REQUIRE(g.adj[1].size() == 2);
}

TEST_CASE("digraph: add_edge directed adjacency", "[graph]") {
	digraph<edge_t> g(3);
	g.add_edge(0, 1);
	g.add_edge(0, 2);
	REQUIRE(g.adj[0].size() == 2);
	REQUIRE(g.adj[1].empty());
	REQUIRE(g.adj[2].empty());
}

TEST_CASE("digraph: edges record from/to", "[graph]") {
	digraph<edge_t> g(3);
	g.add_edge(0, 1);
	g.add_edge(2, 1);
	REQUIRE(g.edges[0].from == 0);
	REQUIRE(g.edges[0].to == 1);
	REQUIRE(g.edges[1].from == 2);
	REQUIRE(g.edges[1].to == 1);
}

TEST_CASE("condi_undigraph: ignore filter", "[graph]") {
	condi_undigraph<edge_t> g(3);
	g.add_edge(0, 1);
	g.add_edge(1, 2);
	REQUIRE_FALSE(g.is_ignore(0));
	g.set_ignore([](const int& id) { return id == 0; });
	REQUIRE(g.is_ignore(0));
	REQUIRE_FALSE(g.is_ignore(1));
}

TEST_CASE("condi_digraph: clear ignore", "[graph]") {
	condi_digraph<edge_t> g(2);
	g.add_edge(0, 1);
	g.set_ignore([](const int&) { return true; });
	REQUIRE(g.is_ignore(0));
	g.clear_ignore<bool>();
	REQUIRE_FALSE(g.is_ignore(0));
}

TEST_CASE("graph: reserve constructor", "[graph]") {
	undigraph<edge_t> g(5, 10);
	for (int i = 0; i < 10; ++ i) g.add_edge(0, 1);
	REQUIRE(g.edges.size() == 10);
}
