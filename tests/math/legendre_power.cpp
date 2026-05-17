#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "math/legendre_power.hpp"

namespace {
static long long ref_legendre(const long long& N, const long long& P) {
	long long res = 0;
	long long pk = P;
	while (pk <= N) {
		res += N / pk;
		if (pk > N / P) break;
		pk *= P;
	}
	return res;
}
} // namespace

TEST_CASE("legendre_power: documented small cases", "[legendre]") {
	REQUIRE(legendre_power<long long>(10, 2) == 8);
	REQUIRE(legendre_power<long long>(10, 5) == 2);
	REQUIRE(legendre_power<long long>(0, 2) == 0);
	REQUIRE(legendre_power<long long>(1, 2) == 0);
	REQUIRE(legendre_power<long long>(25, 5) == 6);
	REQUIRE(legendre_power<long long>(100, 5) == 24);
}

TEST_CASE("legendre_power: P > N gives 0", "[legendre]") {
	REQUIRE(legendre_power<long long>(5, 7) == 0);
	REQUIRE(legendre_power<long long>(100, 101) == 0);
}

TEST_CASE("legendre_power: factorial vs naive prime exponent", "[legendre][stress]") {
	std::mt19937 rng(42);
	for (int t = 0; t < 500; ++ t) {
		long long N = (long long)(rng() % 1000000u);
		long long P = 2 + (long long)(rng() % 30u);
		REQUIRE(legendre_power<long long>(N, P) == ref_legendre(N, P));
	}
}

TEST_CASE("legendre_power: stress small P=2", "[legendre][stress]") {
	for (long long N = 0; N <= 200; ++ N)
		REQUIRE(legendre_power<long long>(N, 2) == ref_legendre(N, 2));
}

TEST_CASE("legendre_power: large N", "[legendre]") {
	REQUIRE(legendre_power<long long>(1000000000000LL, 2) == ref_legendre(1000000000000LL, 2));
	REQUIRE(legendre_power<long long>(1000000000000LL, 7) == ref_legendre(1000000000000LL, 7));
}
