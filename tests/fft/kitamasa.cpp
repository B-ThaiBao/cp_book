#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "mod/modnum.hpp"
#include "fft/kitamasa.hpp"

namespace {
using num = modnum<constant<int, 1000000007>, naive_multiplier<int>>;

static num naive_kth(const std::vector<num>& A0, const std::vector<num>& C, long long k) {
	int n = int(A0.size());
	if (k < n) return A0[k];
	std::vector<num> A = A0;
	for (long long i = n; i <= k; ++ i) {
		num s(0);
		for (int j = 0; j < n; ++ j) s += C[j] * A[A.size() - 1 - j];
		A.push_back(s);
	}
	return A[k];
}
} // namespace

TEST_CASE("kitamasa: nth fibonacci", "[kitamasa]") {
	std::vector<num> A = {num(1), num(1)};
	std::vector<num> C = {num(1), num(1)};
	REQUIRE(kitamasa(A, C, 0) == num(1));
	REQUIRE(kitamasa(A, C, 1) == num(1));
	REQUIRE(kitamasa(A, C, 2) == num(2));
	REQUIRE(kitamasa(A, C, 5) == num(8));
	REQUIRE(kitamasa(A, C, 9) == num(55));
	REQUIRE(kitamasa(A, C, 19) == num(6765));
}

TEST_CASE("kitamasa: geometric", "[kitamasa]") {
	std::vector<num> A = {num(3)};
	std::vector<num> C = {num(2)};
	REQUIRE(kitamasa(A, C, 0) == num(3));
	REQUIRE(kitamasa(A, C, 4) == num(48));
	REQUIRE(kitamasa(A, C, 10) == num(3 * 1024));
	REQUIRE(kitamasa(A, C, 20) == num(3) * num(1 << 20));
}

TEST_CASE("kitamasa: constant recurrence", "[kitamasa]") {
	std::vector<num> A = {num(7)};
	std::vector<num> C = {num(1)};
	for (int k : {0, 1, 5, 100, 10000}) REQUIRE(kitamasa(A, C, k) == num(7));
}

TEST_CASE("kitamasa: tribonacci first few", "[kitamasa]") {
	std::vector<num> A = {num(1), num(1), num(2)};
	std::vector<num> C = {num(1), num(1), num(1)};
	std::vector<long long> expected = {1, 1, 2, 4, 7, 13, 24, 44, 81};
	for (int k = 0; k < int(expected.size()); ++ k) {
		REQUIRE(kitamasa(A, C, k) == num(expected[k]));
	}
}

TEST_CASE("kitamasa: matches naive for small k (stress)", "[kitamasa][stress]") {
	std::mt19937 rng(2024);
	for (int t = 0; t < 8; ++ t) {
		int order = 1 + int(rng() % 5u);
		std::vector<num> C(order), A(order);
		for (auto& c : C) c = num(int(rng() % 100u));
		for (auto& a : A) a = num(int(rng() % 100u));
		for (long long k = 0; k < 15; ++ k) {
			REQUIRE(kitamasa(A, C, k) == naive_kth(A, C, k));
		}
	}
}

TEST_CASE("kitamasa: large k", "[kitamasa]") {
	std::vector<num> A = {num(1), num(1)};
	std::vector<num> C = {num(1), num(1)};
	REQUIRE(kitamasa(A, C, 1000000000LL) == kitamasa(A, C, 1000000000LL));
}

TEST_CASE("kitamasa: a[i]=a[i-1] (identity)", "[kitamasa]") {
	std::vector<num> A = {num(42)};
	std::vector<num> C = {num(1)};
	for (long long k : {0LL, 1LL, 100LL, 1000000LL}) REQUIRE(kitamasa(A, C, k) == num(42));
}

TEST_CASE("kitamasa: alternating sign recurrence", "[kitamasa]") {
	// a[i] = -a[i-1] starting from 5: 5, -5, 5, -5, ...
	std::vector<num> A = {num(5)};
	std::vector<num> C = {num(- 1)};
	REQUIRE(kitamasa(A, C, 0) == num(5));
	REQUIRE(kitamasa(A, C, 1) == num(- 5));
	REQUIRE(kitamasa(A, C, 100) == num(5));
	REQUIRE(kitamasa(A, C, 101) == num(- 5));
}
