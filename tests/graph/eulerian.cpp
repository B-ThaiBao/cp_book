#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "graph/graph.hpp"
#include "graph/eulerian.hpp"

namespace {
struct edge_t { int from, to; };

bool validate_euler(const digraph<edge_t>& g, const int& start, const std::vector<int>& path) {
	if ((int)path.size() != (int)g.edges.size()) return false;
	std::vector<bool> used(g.edges.size(), false);
	int cur = start;
	for (const int& id : path) {
		if (used[id]) return false;
		used[id] = true;
		if (g.edges[id].from != cur) return false;
		cur = g.edges[id].to;
	}
	return true;
}
} // namespace

TEST_CASE("find_eulerian_path: simple cycle", "[eulerian]") {
	digraph<edge_t> g(3);
	g.add_edge(0, 1);
	g.add_edge(1, 2);
	g.add_edge(2, 0);
	auto [rt, path] = find_eulerian_path(g);
	REQUIRE(rt >= 0);
	REQUIRE((int)path.size() == 3);
	REQUIRE(validate_euler(g, rt, path));
}

TEST_CASE("find_eulerian_path: open path", "[eulerian]") {
	digraph<edge_t> g(4);
	g.add_edge(0, 1);
	g.add_edge(1, 2);
	g.add_edge(2, 3);
	auto [rt, path] = find_eulerian_path(g);
	REQUIRE(rt == 0);
	REQUIRE((int)path.size() == 3);
	REQUIRE(validate_euler(g, rt, path));
}

TEST_CASE("find_eulerian_path: no edges", "[eulerian]") {
	digraph<edge_t> g(3);
	auto [rt, path] = find_eulerian_path(g);
	REQUIRE(path.empty());
}

TEST_CASE("find_eulerian_path: double cycle figure-8", "[eulerian]") {
	digraph<edge_t> g(5);
	// 0->1->2->0 and 0->3->4->0
	g.add_edge(0, 1); g.add_edge(1, 2); g.add_edge(2, 0);
	g.add_edge(0, 3); g.add_edge(3, 4); g.add_edge(4, 0);
	auto [rt, path] = find_eulerian_path(g);
	REQUIRE(rt >= 0);
	REQUIRE((int)path.size() == 6);
	REQUIRE(validate_euler(g, rt, path));
}

TEST_CASE("find_eulerian_path: infeasible too many odd", "[eulerian]") {
	digraph<edge_t> g(4);
	g.add_edge(0, 1);
	g.add_edge(0, 2);
	g.add_edge(0, 3);
	auto [rt, path] = find_eulerian_path(g);
	REQUIRE(rt == - 1);
}

TEST_CASE("find_eulerian_path: single self-loop", "[eulerian]") {
	digraph<edge_t> g(2);
	g.add_edge(0, 0);
	auto [rt, path] = find_eulerian_path(g);
	REQUIRE(rt == 0);
	REQUIRE(path.size() == 1);
}
