#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "math/diophantine.hpp"

TEST_CASE("extended_gcd: bezout identity", "[diophantine][stress]") {
	std::mt19937_64 rng(11);
	for (int t = 0; t < 1000; ++ t) {
		long long a = (long long)(rng() % 2001) - 1000;
		long long b = (long long)(rng() % 2001) - 1000;
		if (a == 0 && b == 0) continue;
		long long x, y;
		long long g = extended_gcd<long long>(a, b, x, y);
		REQUIRE(a * x + b * y == g);
		REQUIRE(std::abs(g) == (long long)std::gcd(std::abs(a), std::abs(b)));
	}
}

TEST_CASE("extended_gcd: zero cases", "[diophantine]") {
	long long x, y;
	REQUIRE(extended_gcd<long long>(0, 5, x, y) == 5);
	REQUIRE(0 * x + 5 * y == 5);
	REQUIRE(extended_gcd<long long>(7, 0, x, y) == 7);
	REQUIRE(7 * x + 0 * y == 7);
}

TEST_CASE("diophantine: solvable case", "[diophantine]") {
	long long x, y, g;
	REQUIRE(diophantine<long long>(6, 9, 15, x, y, g));
	REQUIRE(6 * x + 9 * y == 15);
	REQUIRE(g == 3);
}

TEST_CASE("diophantine: unsolvable case", "[diophantine]") {
	long long x, y, g;
	REQUIRE_FALSE(diophantine<long long>(6, 9, 1, x, y, g));
}

TEST_CASE("diophantine: trivial zero=zero", "[diophantine]") {
	long long x, y, g;
	REQUIRE(diophantine<long long>(0, 0, 0, x, y, g));
	REQUIRE(x == 0);
	REQUIRE(y == 0);
	REQUIRE(g == 0);
}

TEST_CASE("diophantine: zero coefficients with nonzero target unsolvable", "[diophantine]") {
	long long x, y, g;
	REQUIRE_FALSE(diophantine<long long>(0, 0, 5, x, y, g));
}

TEST_CASE("diophantine: one zero coefficient", "[diophantine]") {
	long long x, y, g;
	REQUIRE(diophantine<long long>(0, 5, 10, x, y, g));
	REQUIRE(y == 2);
	REQUIRE(g == 5);
}

TEST_CASE("diophantine: stress", "[diophantine][stress]") {
	std::mt19937_64 rng(13);
	for (int t = 0; t < 800; ++ t) {
		long long a = (long long)(rng() % 201) - 100;
		long long b = (long long)(rng() % 201) - 100;
		long long c = (long long)(rng() % 4001) - 2000;
		long long x, y, g;
		bool ok = diophantine<long long>(a, b, c, x, y, g);
		if (ok) {
			REQUIRE(a * x + b * y == c);
		} else {
			if (a == 0 && b == 0) {
				REQUIRE(c != 0);
			} else {
				long long gg = (long long)std::gcd(std::abs(a), std::abs(b));
				REQUIRE(c % gg != 0);
			}
		}
	}
}
