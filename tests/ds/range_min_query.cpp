#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "ds/range_min_query.hpp"

TEST_CASE("range_min_query: leftmost min basic", "[rmq]") {
	std::vector<int> A = {3, 1, 4, 1, 5, 9, 2, 6, 5, 3};
	range_min_query<int> rmq(A);
	for (int l = 0; l < (int)A.size(); ++ l) {
		for (int r = l + 1; r <= (int)A.size(); ++ r) {
			auto [idx, val] = rmq.range_query(l, r);
			int expect_val = *std::min_element(A.begin() + l, A.begin() + r);
			REQUIRE(val == expect_val);
			REQUIRE(A[idx] == expect_val);
			int expect_idx = int(std::find(A.begin() + l, A.begin() + r, expect_val) - A.begin());
			REQUIRE(idx == expect_idx);
		}
	}
}

TEST_CASE("range_min_query: rightmost via less_equal", "[rmq]") {
	std::vector<int> A = {3, 1, 4, 1, 5, 9, 2, 6, 5, 3};
	range_min_query<int, std::less_equal<int>> rmq(A);
	for (int l = 0; l < (int)A.size(); ++ l) {
		for (int r = l + 1; r <= (int)A.size(); ++ r) {
			auto [idx, val] = rmq.range_query(l, r);
			int expect_val = *std::min_element(A.begin() + l, A.begin() + r);
			REQUIRE(val == expect_val);
			int last = - 1;
			for (int i = l; i < r; ++ i) if (A[i] == expect_val) last = i;
			REQUIRE(idx == last);
		}
	}
}

TEST_CASE("range_min_query: max via greater", "[rmq]") {
	std::vector<int> A = {3, 1, 4, 1, 5, 9, 2, 6};
	range_min_query<int, std::greater<int>> rmax(A);
	for (int l = 0; l < (int)A.size(); ++ l) {
		for (int r = l + 1; r <= (int)A.size(); ++ r) {
			auto [idx, val] = rmax.range_query(l, r);
			int expect = *std::max_element(A.begin() + l, A.begin() + r);
			REQUIRE(val == expect);
			REQUIRE(A[idx] == expect);
		}
	}
}

TEST_CASE("range_min_query: single element", "[rmq]") {
	std::vector<int> A = {7};
	range_min_query<int> rmq(A);
	auto [idx, val] = rmq.range_query(0, 1);
	REQUIRE(idx == 0);
	REQUIRE(val == 7);
}

TEST_CASE("range_min_query: large random stress", "[rmq][stress]") {
	std::mt19937 rng(2024);
	int N = 500;
	std::vector<int> A(N);
	for (auto& x : A) x = int(rng() % 100000);
	range_min_query<int> rmq(A);
	for (int t = 0; t < 3000; ++ t) {
		int l = int(rng() % N), r = int(rng() % N);
		if (l > r) std::swap(l, r);
		if (l == r) ++ r;
		auto [idx, val] = rmq.range_query(l, r);
		int expect = *std::min_element(A.begin() + l, A.begin() + r);
		REQUIRE(val == expect);
	}
}
