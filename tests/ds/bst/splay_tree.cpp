#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "ds/bst/splay_tree.hpp"

namespace {
std::vector<splay_tree_node*> collect_inorder(splay_tree_node* v) {
	std::vector<splay_tree_node*> out;
	std::function<void(splay_tree_node*)> dfs = [&](splay_tree_node* u) {
		if (u == nullptr) return;
		u -> downdate();
		dfs(u -> c[0]);
		out.push_back(u);
		dfs(u -> c[1]);
	};
	dfs(v);
	return out;
}
std::vector<int> linearize(splay_tree_node* v, const std::vector<splay_tree_node>& nodes) {
	std::vector<int> ids;
	auto seq = collect_inorder(v);
	for (auto p : seq) ids.push_back(int(p - nodes.data()));
	return ids;
}
} // namespace

TEST_CASE("splay_tree: build + find_implicit / find_pos", "[splay]") {
	int N = 8;
	std::vector<splay_tree_node> nodes(N);
	build(nodes);
	splay_tree_node* root = &nodes[0];
	for (int k = 0; k < N; ++ k) {
		auto p = find_implicit(root, k);
		REQUIRE(p == &nodes[k]);
		REQUIRE(find_pos(p) == k);
		root = find_root(p);
	}
}

TEST_CASE("splay_tree: split + merge round-trip", "[splay]") {
	int N = 6;
	std::vector<splay_tree_node> nodes(N);
	build(nodes);
	splay_tree_node* root = &nodes[0];
	auto [l, r] = split_implicit(root, 3);
	REQUIRE(size(l) == 3);
	REQUIRE(size(r) == 3);
	auto merged = merge(l, r);
	REQUIRE(linearize(merged, nodes) == std::vector<int>({0, 1, 2, 3, 4, 5}));
}

TEST_CASE("splay_tree: reverse via flip", "[splay]") {
	int N = 6;
	std::vector<splay_tree_node> nodes(N);
	build(nodes);
	splay_tree_node* root = &nodes[0];
	auto [l, mid_r] = split_implicit(root, 1);
	auto [mid, r] = split_implicit(mid_r, 4);
	mid -> reverse();
	auto merged = merge(merge(l, mid), r);
	REQUIRE(linearize(merged, nodes) == std::vector<int>({0, 4, 3, 2, 1, 5}));
}

TEST_CASE("splay_tree: erase by position", "[splay]") {
	int N = 5;
	std::vector<splay_tree_node> nodes(N);
	build(nodes);
	splay_tree_node* root = &nodes[0];
	auto p = find_implicit(root, 2);
	root = erase(p);
	REQUIRE(size(root) == 4);
	REQUIRE(linearize(root, nodes) == std::vector<int>({0, 1, 3, 4}));
}

TEST_CASE("splay_tree: single node", "[splay]") {
	std::vector<splay_tree_node> nodes(1);
	build(nodes);
	REQUIRE(size(&nodes[0]) == 1);
	REQUIRE(find_implicit(&nodes[0], 0) == &nodes[0]);
}
