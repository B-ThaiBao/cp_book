#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "geo/pair_point.hpp"

namespace {
using pt = pair_point<long long>;
} // namespace

TEST_CASE("pair_point: arithmetic", "[pair_point]") {
	pt a(1, 2), b(3, 4);
	REQUIRE(a + b == pt(4, 6));
	REQUIRE(b - a == pt(2, 2));
	REQUIRE(a * 3 == pt(3, 6));
	REQUIRE(2 * a == pt(2, 4));
}

TEST_CASE("pair_point: comparisons", "[pair_point]") {
	pt a(1, 2), b(3, 4);
	REQUIRE(a == pt(1, 2));
	REQUIRE_FALSE(a == b);
	REQUIRE(a < b);
	REQUIRE_FALSE(b < a);
}

TEST_CASE("pair_point: dot and cross", "[pair_point]") {
	pt a(1, 0), b(0, 1);
	REQUIRE(a.dot(b) == 0);
	REQUIRE(a.cross(b) == 1);
	REQUIRE(b.cross(a) == - 1);
	pt c(2, 3), d(4, 5);
	REQUIRE(c.dot(d) == 8 + 15);
	REQUIRE(c.cross(d) == 10 - 12);
}

TEST_CASE("pair_point: norm", "[pair_point]") {
	pt a(3, 4);
	REQUIRE(a.norm() == 25);
	REQUIRE(a.norm(pt(0, 0)) == 25);
	REQUIRE(a.norm(pt(3, 0)) == 16);
}

TEST_CASE("pair_point: perp", "[pair_point]") {
	pt a(2, 3);
	REQUIRE(a.perp_cw() == pt(3, - 2));
	REQUIRE(a.perp_ccw() == pt(- 3, 2));
}

TEST_CASE("pair_point: cross3 orientation", "[pair_point]") {
	pt a(0, 0), b(1, 0), c(0, 1), d(2, 0);
	REQUIRE(cross3(a, b, c) > 0);
	REQUIRE(cross3(a, c, b) < 0);
	REQUIRE(cross3(a, b, d) == 0);
}

TEST_CASE("pair_point: stress dot+cross properties", "[pair_point][stress]") {
	std::mt19937 rng(11);
	for (int t = 0; t < 200; ++ t) {
		long long ax = (long long)(rng() % 201) - 100;
		long long ay = (long long)(rng() % 201) - 100;
		long long bx = (long long)(rng() % 201) - 100;
		long long by = (long long)(rng() % 201) - 100;
		pt a(ax, ay), b(bx, by);
		REQUIRE(a.dot(b) == b.dot(a));
		REQUIRE(a.cross(b) == - b.cross(a));
		REQUIRE((a + b).dot(a - b) == a.dot(a) - b.dot(b));
	}
}
