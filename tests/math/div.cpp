#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "math/div.hpp"

TEST_CASE("floor_div: basic mixed signs", "[div]") {
	REQUIRE(floor_div(10, 3) == 3);
	REQUIRE(floor_div(- 10, 3) == - 4);
	REQUIRE(floor_div(10, - 3) == - 4);
	REQUIRE(floor_div(- 10, - 3) == 3);
	REQUIRE(floor_div(9, 3) == 3);
	REQUIRE(floor_div(- 9, 3) == - 3);
}

TEST_CASE("ceil_div: basic mixed signs", "[div]") {
	REQUIRE(ceil_div(10, 3) == 4);
	REQUIRE(ceil_div(- 10, 3) == - 3);
	REQUIRE(ceil_div(10, - 3) == - 3);
	REQUIRE(ceil_div(- 10, - 3) == 4);
	REQUIRE(ceil_div(9, 3) == 3);
	REQUIRE(ceil_div(- 9, 3) == - 3);
}

TEST_CASE("floor_div: exact division", "[div]") {
	for (int a = - 20; a <= 20; ++ a)
		for (int b : {- 5, - 2, - 1, 1, 2, 5}) {
			if (a % b == 0) REQUIRE(floor_div(a, b) == a / b);
		}
}

TEST_CASE("floor/ceil_div: stress vs std::floor/ceil", "[div][stress]") {
	std::mt19937 rng(7);
	for (int t = 0; t < 2000; ++ t) {
		int a = int(rng() % 2001u) - 1000;
		int b = int(rng() % 2001u) - 1000;
		if (b == 0) continue;
		long long f = floor_div<long long, long long>(a, b);
		long long c = ceil_div<long long, long long>(a, b);
		long long ef = (long long)std::floor((double)a / b);
		long long ec = (long long)std::ceil((double)a / b);
		REQUIRE(f == ef);
		REQUIRE(c == ec);
		if (a % b == 0) REQUIRE(f == c);
		else REQUIRE(c == f + 1);
	}
}

TEST_CASE("floor_div: large values", "[div]") {
	REQUIRE(floor_div<long long, long long>(1000000000000LL, 7LL) == 142857142857LL);
	REQUIRE(floor_div<long long, long long>(- 1000000000000LL, 7LL) == - 142857142858LL);
}

TEST_CASE("ceil_div: large values", "[div]") {
	REQUIRE(ceil_div<long long, long long>(1000000000001LL, 7LL) == 142857142858LL);
	REQUIRE(ceil_div<long long, long long>(- 1000000000001LL, 7LL) == - 142857142857LL);
}
