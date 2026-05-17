#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "ds/binary_indexed_tree.hpp"

namespace {
auto add_fn = [](auto& a, const auto& b) { a += b; };
} // namespace

TEST_CASE("BIT: prefix tree point update + prefix query", "[bit]") {
	int N = 8;
	binary_indexed_tree<long long> bit(N, 0);
	auto point_update = [&](const int& i, const long long& v) {
		for (auto& x : bit.suffix(i)) x += v;
	};
	auto prefix_query = [&](const int& r) -> long long {
		long long s = 0;
		for (auto& x : bit.prefix(r)) s += x;
		return s;
	};
	std::vector<long long> ref(N, 0);
	std::mt19937 rng(7);
	for (int t = 0; t < 400; ++ t) {
		int op = rng() % 2;
		if (op == 0) {
			int i = rng() % N;
			long long v = (long long)(rng() % 200) - 100;
			point_update(i, v);
			ref[i] += v;
		} else {
			int r = rng() % (N + 1);
			long long expect = 0;
			for (int j = 0; j < r; ++ j) expect += ref[j];
			REQUIRE(prefix_query(r) == expect);
		}
	}
}

TEST_CASE("BIT: build_prefix matches incremental", "[bit]") {
	std::mt19937 rng(11);
	int N = 16;
	std::vector<long long> A(N);
	for (auto& x : A) x = (long long)(rng() % 50);
	binary_indexed_tree<long long> bit1(N, 0);
	for (int i = 0; i < N; ++ i) for (auto& x : bit1.suffix(i)) x += A[i];
	binary_indexed_tree<long long> bit2;
	bit2.build_prefix(A, add_fn);
	for (int r = 0; r <= N; ++ r) {
		long long s1 = 0, s2 = 0;
		for (auto& x : bit1.prefix(r)) s1 += x;
		for (auto& x : bit2.prefix(r)) s2 += x;
		REQUIRE(s1 == s2);
	}
}

TEST_CASE("BIT: suffix tree point + suffix query", "[bit]") {
	int N = 8;
	binary_indexed_tree<long long> bit(N, 0);
	auto point_update = [&](const int& i, const long long& v) {
		for (auto& x : bit.prefix(i + 1)) x += v;
	};
	auto suffix_query = [&](const int& l) -> long long {
		long long s = 0;
		for (auto& x : bit.suffix(l)) s += x;
		return s;
	};
	std::vector<long long> ref(N, 0);
	std::mt19937 rng(17);
	for (int t = 0; t < 400; ++ t) {
		int op = rng() % 2;
		if (op == 0) {
			int i = rng() % N;
			long long v = (long long)(rng() % 100);
			point_update(i, v);
			ref[i] += v;
		} else {
			int l = rng() % (N + 1);
			long long expect = 0;
			for (int j = l; j < N; ++ j) expect += ref[j];
			REQUIRE(suffix_query(l) == expect);
		}
	}
}

TEST_CASE("BIT: lower_bound / upper_bound", "[bit]") {
	int N = 16;
	binary_indexed_tree<long long> bit(N, 0);
	std::vector<long long> A = {1, 0, 2, 0, 3, 0, 1, 4, 0, 0, 5, 0, 1, 1, 1, 1};
	bit.build_prefix(A, add_fn);
	std::vector<long long> pref(N + 1, 0);
	for (int i = 0; i < N; ++ i) pref[i + 1] = pref[i] + A[i];
	for (long long v = - 1; v <= pref[N] + 2; ++ v) {
		int got_lb = bit.lower_bound(v);
		int expect_lb = int(std::lower_bound(pref.begin() + 1, pref.end(), v) - pref.begin());
		if (expect_lb > N) expect_lb = N + 1;
		REQUIRE(got_lb == expect_lb);
		int got_ub = bit.upper_bound(v);
		int expect_ub = int(std::upper_bound(pref.begin() + 1, pref.end(), v) - pref.begin());
		if (expect_ub > N) expect_ub = N + 1;
		REQUIRE(got_ub == expect_ub);
	}
}

TEST_CASE("BIT: empty tree zero queries", "[bit]") {
	binary_indexed_tree<long long> bit(5, 0);
	long long s = 0;
	for (auto& x : bit.prefix(5)) s += x;
	REQUIRE(s == 0);
}

TEST_CASE("BIT: single element", "[bit]") {
	binary_indexed_tree<long long> bit(1, 0);
	for (auto& x : bit.suffix(0)) x += 42;
	long long s = 0;
	for (auto& x : bit.prefix(1)) s += x;
	REQUIRE(s == 42);
}

TEST_CASE("range_update_sum_query: range update + range query", "[bit][rusq]") {
	int N = 12;
	range_update_sum_query<long long> rusq(N);
	std::vector<long long> ref(N, 0);
	std::mt19937 rng(31);
	for (int t = 0; t < 500; ++ t) {
		int op = rng() % 2;
		if (op == 0) {
			int l = int(rng() % unsigned(N + 1)), r = int(rng() % unsigned(N + 1));
			if (l > r) std::swap(l, r);
			long long v = (long long)(rng() % 50) - 25;
			rusq.range_update(l, r, v);
			for (int i = l; i < r; ++ i) ref[i] += v;
		} else {
			int l = int(rng() % unsigned(N + 1)), r = int(rng() % unsigned(N + 1));
			if (l > r) std::swap(l, r);
			long long expect = 0;
			for (int i = l; i < r; ++ i) expect += ref[i];
			REQUIRE(rusq.range_query(l, r) == expect);
		}
	}
}

TEST_CASE("range_update_sum_query: build from vector", "[bit][rusq]") {
	std::vector<long long> A = {3, 1, 4, 1, 5, 9, 2, 6};
	range_update_sum_query<long long> rusq(A);
	long long sum = 0;
	for (int i = 0; i < (int)A.size(); ++ i) {
		sum += A[i];
		REQUIRE(rusq.prefix_sum(i + 1) == sum);
	}
	REQUIRE(rusq.range_query(2, 5) == A[2] + A[3] + A[4]);
}
