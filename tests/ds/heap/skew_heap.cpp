#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "ds/heap/skew_heap.hpp"

namespace {
struct min_node : public skew_heap_node_base<min_node> {
	int key = 0;
	min_node() = default;
	min_node(const int& k) : key(k) {}
	void do_downdate() {}
};

auto cmp = [](const min_node* a, const min_node* b) { return a -> key <= b -> key; };
} // namespace

TEST_CASE("skew_heap: merge two heaps", "[skew_heap]") {
	min_node a(5), b(3);
	auto* root = merge(static_cast<min_node*>(&a), static_cast<min_node*>(&b), cmp);
	REQUIRE(root -> key == 3);
}

TEST_CASE("skew_heap: insert sequence pop sorted", "[skew_heap]") {
	std::vector<min_node> nodes;
	std::vector<int> keys = {7, 3, 5, 1, 9, 2, 4, 6, 8};
	nodes.reserve(keys.size());
	for (int k : keys) nodes.emplace_back(k);
	min_node* root = nullptr;
	for (auto& n : nodes) root = merge(root, &n, cmp);
	std::vector<int> got;
	while (root != nullptr) {
		got.push_back(root -> key);
		root = pop(root, cmp);
	}
	auto sorted = keys;
	std::sort(sorted.begin(), sorted.end());
	REQUIRE(got == sorted);
}

TEST_CASE("skew_heap: empty merge identity", "[skew_heap]") {
	min_node a(42);
	REQUIRE(merge(static_cast<min_node*>(nullptr), &a, cmp) == &a);
	REQUIRE(merge(&a, static_cast<min_node*>(nullptr), cmp) == &a);
}

TEST_CASE("skew_heap: stress vs priority_queue", "[skew_heap][stress]") {
	std::mt19937 rng(7);
	int N = 200;
	std::vector<min_node> nodes;
	nodes.reserve(N);
	std::priority_queue<int, std::vector<int>, std::greater<int>> ref;
	min_node* root = nullptr;
	for (int i = 0; i < N; ++ i) {
		int k = int(rng() % 10000);
		nodes.emplace_back(k);
		root = merge(root, &nodes.back(), cmp);
		ref.push(k);
	}
	std::vector<int> got;
	while (root != nullptr) {
		got.push_back(root -> key);
		root = pop(root, cmp);
	}
	std::vector<int> expect;
	while (!ref.empty()) { expect.push_back(ref.top()); ref.pop(); }
	REQUIRE(got == expect);
}
