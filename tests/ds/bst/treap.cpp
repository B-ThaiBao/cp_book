#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "ds/bst/treap.hpp"

namespace {
std::vector<treap_node*> collect_inorder(treap_node* v) {
	std::vector<treap_node*> out;
	std::function<void(treap_node*)> dfs = [&](treap_node* u) {
		if (u == nullptr) return;
		u -> downdate();
		dfs(u -> c[0]);
		out.push_back(u);
		dfs(u -> c[1]);
	};
	dfs(v);
	return out;
}
std::vector<int> linearize(treap_node* v, const std::vector<treap_node>& nodes) {
	std::vector<int> ids;
	auto seq = collect_inorder(v);
	for (auto p : seq) ids.push_back(int(p - nodes.data()));
	return ids;
}
} // namespace

TEST_CASE("treap: build_by_heap + find_implicit", "[treap]") {
	int N = 8;
	std::vector<treap_node> nodes(N);
	int rid = build_by_heap(nodes);
	treap_node* root = &nodes[rid];
	for (int k = 0; k < N; ++ k) {
		auto p = find_implicit(root, k);
		REQUIRE(p == &nodes[k]);
	}
	REQUIRE(linearize(root, nodes) == std::vector<int>({0, 1, 2, 3, 4, 5, 6, 7}));
}

TEST_CASE("treap: split + merge round-trip", "[treap]") {
	int N = 7;
	std::vector<treap_node> nodes(N);
	int rid = build_by_heap(nodes);
	treap_node* root = &nodes[rid];
	auto [l, r] = split_implicit(root, 3);
	REQUIRE(size(l) == 3);
	REQUIRE(size(r) == 4);
	auto merged = merge(l, r);
	REQUIRE(linearize(merged, nodes) == std::vector<int>({0, 1, 2, 3, 4, 5, 6}));
}

TEST_CASE("treap: reverse via flip", "[treap]") {
	int N = 6;
	std::vector<treap_node> nodes(N);
	int rid = build_by_heap(nodes);
	treap_node* root = &nodes[rid];
	auto [l, mid_r] = split_implicit(root, 1);
	auto [mid, r] = split_implicit(mid_r, 4);
	mid -> reverse();
	auto merged = merge(merge(l, mid), r);
	REQUIRE(linearize(merged, nodes) == std::vector<int>({0, 4, 3, 2, 1, 5}));
}

TEST_CASE("treap: erase node", "[treap]") {
	int N = 5;
	std::vector<treap_node> nodes(N);
	int rid = build_by_heap(nodes);
	treap_node* root = &nodes[rid];
	auto p = find_implicit(root, 2);
	root = erase(p);
	REQUIRE(size(root) == 4);
	REQUIRE(linearize(root, nodes) == std::vector<int>({0, 1, 3, 4}));
}

TEST_CASE("treap: next/previous traversal", "[treap]") {
	int N = 5;
	std::vector<treap_node> nodes(N);
	build_by_heap(nodes);
	treap_node* cur = &nodes[0];
	while (cur -> par != nullptr) cur = cur -> par;
	auto first = collect_inorder(cur).front();
	std::vector<treap_node*> seq;
	for (auto x = first; x != nullptr; x = next(x)) seq.push_back(x);
	REQUIRE((int)seq.size() == N);
	for (int i = 0; i < N; ++ i) REQUIRE(seq[i] == &nodes[i]);
}
