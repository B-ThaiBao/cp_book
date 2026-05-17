#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "math/derangement.hpp"

TEST_CASE("count_derangement: first 10 values", "[derangement]") {
	auto d = count_derangement<long long>(10);
	std::vector<long long> ref = {1, 0, 1, 2, 9, 44, 265, 1854, 14833, 133496};
	REQUIRE(d == ref);
}

TEST_CASE("count_derangement: empty input", "[derangement]") {
	REQUIRE(count_derangement<long long>(0).empty());
}

TEST_CASE("count_derangement: n=1", "[derangement]") {
	REQUIRE(count_derangement<long long>(1) == std::vector<long long>({0}));
}

TEST_CASE("count_derangement: recurrence d[n] = (n-1)*(d[n-1]+d[n-2])", "[derangement][stress]") {
	auto d = count_derangement<long long>(15);
	for (int n = 2; n < 15; ++ n)
		REQUIRE(d[n] == (long long)(n - 1) * (d[n - 1] + d[n - 2]));
}

TEST_CASE("count_derangement: alternative recurrence d[n] = n*d[n-1] + (-1)^n", "[derangement][stress]") {
	auto d = count_derangement<long long>(15);
	for (int n = 1; n < 15; ++ n) {
		long long sign = (n % 2 == 0) ? 1 : - 1;
		REQUIRE(d[n] == (long long)n * d[n - 1] + sign);
	}
}

TEST_CASE("count_derangement: D(3)=2, D(4)=9", "[derangement]") {
	auto d = count_derangement<long long>(5);
	REQUIRE(d[3] == 2);
	REQUIRE(d[4] == 9);
}
