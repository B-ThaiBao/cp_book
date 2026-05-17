#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "ds/sparse_table.hpp"

TEST_CASE("sparse_table: range min stress", "[sparse_table][stress]") {
	std::mt19937 rng(2024);
	int N = 200;
	std::vector<int> A(N);
	for (auto& x : A) x = int(rng() % 10000);
	sparse_table<int> st(N);
	for (int i = 0; i < N; ++ i) st(0, i) = A[i];
	st.build([&](auto p) { st[p] = std::min(st[p.c(0)], st[p.c(1)]); });
	for (int t = 0; t < 1500; ++ t) {
		int l = int(rng() % unsigned(N)), r = int(rng() % unsigned(N));
		if (l > r) std::swap(l, r);
		if (l == r) ++ r;
		auto rng_pair = st.range(l, r);
		int got = std::min(st[rng_pair[0]], st[rng_pair[1]]);
		int expect = *std::min_element(A.begin() + l, A.begin() + r);
		REQUIRE(got == expect);
	}
}

TEST_CASE("sparse_table: for_range sum", "[sparse_table]") {
	int N = 32;
	std::vector<long long> A(N);
	for (int i = 0; i < N; ++ i) A[i] = i + 1;
	sparse_table<long long> st(N);
	for (int i = 0; i < N; ++ i) st(0, i) = A[i];
	st.build([&](auto p) { st[p] = st[p.c(0)] + st[p.c(1)]; });
	for (int l = 0; l < N; ++ l) {
		for (int r = l + 1; r <= N; ++ r) {
			long long got = 0;
			st.for_range(l, r, [&](auto p) { got += st[p]; });
			long long expect = 0;
			for (int j = l; j < r; ++ j) expect += A[j];
			REQUIRE(got == expect);
		}
	}
}

TEST_CASE("sparse_table: for_reverse_range concatenation", "[sparse_table]") {
	int N = 16;
	std::vector<std::string> A(N);
	for (int i = 0; i < N; ++ i) A[i] = std::string(1, char('a' + i));
	sparse_table<std::string> st(N);
	for (int i = 0; i < N; ++ i) st(0, i) = A[i];
	st.build([&](auto p) { st[p] = st[p.c(0)] + st[p.c(1)]; });
	for (int l = 0; l < N; ++ l) {
		for (int r = l + 1; r <= N; ++ r) {
			std::string got;
			st.for_reverse_range(l, r, [&](auto p) { got = st[p] + got; });
			std::string expect;
			for (int j = l; j < r; ++ j) expect += A[j];
			REQUIRE(got == expect);
		}
	}
}

TEST_CASE("sparse_table: single element", "[sparse_table]") {
	sparse_table<int> st(1);
	st(0, 0) = 42;
	st.build([&](auto) {});
	auto rng_pair = st.range(0, 1);
	REQUIRE(st[rng_pair[0]] == 42);
}
