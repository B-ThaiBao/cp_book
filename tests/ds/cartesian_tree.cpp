#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "ds/cartesian_tree.hpp"

namespace {
std::vector<int> naive_cartesian(const std::vector<int>& A) {
	int N = (int)A.size();
	std::vector<int> par(N, - 1);
	std::function<void(int, int, int)> rec = [&](int l, int r, int p) {
		if (l > r) return;
		int idx = l;
		for (int i = l + 1; i <= r; ++ i) if (A[i] < A[idx]) idx = i;
		par[idx] = p;
		rec(l, idx - 1, idx);
		rec(idx + 1, r, idx);
	};
	rec(0, N - 1, - 1);
	return par;
}
} // namespace

TEST_CASE("cartesian_tree: leftmost min basic", "[cartesian]") {
	std::vector<int> A = {3, 1, 4, 1, 5, 9, 2, 6};
	auto par = build_cartesian_tree(A);
	REQUIRE(par == naive_cartesian(A));
}

TEST_CASE("cartesian_tree: single element", "[cartesian]") {
	std::vector<int> A = {7};
	auto par = build_cartesian_tree(A);
	REQUIRE(par == std::vector<int>({- 1}));
}

TEST_CASE("cartesian_tree: sorted ascending", "[cartesian]") {
	std::vector<int> A = {1, 2, 3, 4, 5};
	auto par = build_cartesian_tree(A);
	REQUIRE(par == naive_cartesian(A));
}

TEST_CASE("cartesian_tree: sorted descending", "[cartesian]") {
	std::vector<int> A = {5, 4, 3, 2, 1};
	auto par = build_cartesian_tree(A);
	REQUIRE(par == naive_cartesian(A));
}

TEST_CASE("cartesian_tree: random stress", "[cartesian][stress]") {
	std::mt19937 rng(42);
	for (int t = 0; t < 100; ++ t) {
		int N = 1 + int(rng() % 50);
		std::vector<int> A(N);
		for (auto& x : A) x = int(rng() % 30);
		REQUIRE(build_cartesian_tree(A) == naive_cartesian(A));
	}
}

TEST_CASE("cartesian_tree: leftmost vs rightmost min", "[cartesian]") {
	std::vector<int> A = {1, 1, 1};
	auto lmt = build_cartesian_tree(A, std::less<>());
	auto rmt = build_cartesian_tree(A, std::less_equal<>());
	int lroot = - 1, rroot = - 1;
	for (int i = 0; i < (int)A.size(); ++ i) {
		if (lmt[i] == - 1) lroot = i;
		if (rmt[i] == - 1) rroot = i;
	}
	REQUIRE(lroot == 0);
	REQUIRE(rroot == 2);
}

TEST_CASE("cartesian_tree: empty input", "[cartesian]") {
	std::vector<int> A;
	REQUIRE(build_cartesian_tree(A).empty());
}
