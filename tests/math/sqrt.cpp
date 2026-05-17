#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "math/sqrt.hpp"

TEST_CASE("floor_sqrt: small values", "[sqrt]") {
	REQUIRE(floor_sqrt<int64_t>(0) == 0);
	REQUIRE(floor_sqrt<int64_t>(1) == 1);
	REQUIRE(floor_sqrt<int64_t>(2) == 1);
	REQUIRE(floor_sqrt<int64_t>(3) == 1);
	REQUIRE(floor_sqrt<int64_t>(4) == 2);
	REQUIRE(floor_sqrt<int64_t>(99) == 9);
	REQUIRE(floor_sqrt<int64_t>(100) == 10);
	REQUIRE(floor_sqrt<int64_t>(101) == 10);
}

TEST_CASE("ceil_sqrt: small values", "[sqrt]") {
	REQUIRE(ceil_sqrt<int64_t>(0) == 0);
	REQUIRE(ceil_sqrt<int64_t>(1) == 1);
	REQUIRE(ceil_sqrt<int64_t>(2) == 2);
	REQUIRE(ceil_sqrt<int64_t>(3) == 2);
	REQUIRE(ceil_sqrt<int64_t>(4) == 2);
	REQUIRE(ceil_sqrt<int64_t>(99) == 10);
	REQUIRE(ceil_sqrt<int64_t>(100) == 10);
	REQUIRE(ceil_sqrt<int64_t>(101) == 11);
}

TEST_CASE("floor_sqrt: dense small range", "[sqrt][stress]") {
	for (int64_t n = 0; n <= 1000; ++ n) {
		int64_t f = floor_sqrt<int64_t>(n);
		REQUIRE(f * f <= n);
		REQUIRE((f + 1) * (f + 1) > n);
	}
}

TEST_CASE("ceil_sqrt: dense small range", "[sqrt][stress]") {
	for (int64_t n = 0; n <= 1000; ++ n) {
		int64_t c = ceil_sqrt<int64_t>(n);
		REQUIRE(c * c >= n);
		if (c > 0) REQUIRE((c - 1) * (c - 1) < n);
	}
}

TEST_CASE("floor_sqrt: large random stress", "[sqrt][stress]") {
	std::mt19937_64 rng(1);
	for (int t = 0; t < 500; ++ t) {
		int64_t n = int64_t(rng() % uint64_t(1000000000000LL));
		int64_t f = floor_sqrt<int64_t>(n);
		REQUIRE(f * f <= n);
		REQUIRE((f + 1) * (f + 1) > n);
	}
}

TEST_CASE("ceil_sqrt: large random stress", "[sqrt][stress]") {
	std::mt19937_64 rng(2);
	for (int t = 0; t < 500; ++ t) {
		int64_t n = int64_t(rng() % uint64_t(1000000000000LL));
		int64_t c = ceil_sqrt<int64_t>(n);
		REQUIRE(c * c >= n);
		if (c > 0) REQUIRE((c - 1) * (c - 1) < n);
	}
}

TEST_CASE("floor_sqrt: perfect squares", "[sqrt][stress]") {
	for (int64_t k = 0; k <= 100000; k += 7) {
		REQUIRE(floor_sqrt<int64_t>(k * k) == k);
		REQUIRE(ceil_sqrt<int64_t>(k * k) == k);
	}
}
