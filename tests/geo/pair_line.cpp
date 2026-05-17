#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "geo/pair_point.hpp"
#include "geo/pair_line.hpp"

namespace {
using pt = pair_point<long long>;
using ln = pair_line<long long>;
} // namespace

TEST_CASE("pair_line: on_line basic", "[pair_line]") {
	ln l(pt(0, 0), pt(4, 4));
	REQUIRE(l.on_line(pt(2, 2)));
	REQUIRE(l.on_line(pt(5, 5)));
	REQUIRE(l.on_line(pt(- 1, - 1)));
	REQUIRE_FALSE(l.on_line(pt(1, 2)));
}

TEST_CASE("pair_line: on_segment basic", "[pair_line]") {
	ln l(pt(0, 0), pt(4, 4));
	REQUIRE(l.on_segment(pt(2, 2)));
	REQUIRE(l.on_segment(pt(0, 0)));
	REQUIRE(l.on_segment(pt(4, 4)));
	REQUIRE_FALSE(l.on_segment(pt(5, 5)));
	REQUIRE_FALSE(l.on_segment(pt(- 1, - 1)));
}

TEST_CASE("pair_line: horizontal segment on_line", "[pair_line]") {
	ln l(pt(0, 5), pt(10, 5));
	REQUIRE(l.on_line(pt(- 100, 5)));
	REQUIRE_FALSE(l.on_line(pt(0, 6)));
}

TEST_CASE("pair_line: intersect crossing", "[pair_line]") {
	ln a(pt(0, 0), pt(4, 0));
	ln b(pt(2, - 1), pt(2, 1));
	auto r = intersect(a, b);
	REQUIRE(r[2] != 0);
	REQUIRE(r[0] * 2 == r[2]);
}

TEST_CASE("pair_line: intersect parallel non-coincident", "[pair_line]") {
	ln a(pt(0, 0), pt(4, 0));
	ln b(pt(0, 1), pt(4, 1));
	auto r = intersect(a, b);
	REQUIRE(r[2] == 0);
	REQUIRE(r[0] != 0);
}

TEST_CASE("pair_line: intersect coincident", "[pair_line]") {
	ln a(pt(0, 0), pt(4, 0));
	ln b(pt(1, 0), pt(5, 0));
	auto r = intersect(a, b);
	REQUIRE(r[2] == 0);
	REQUIRE(r[0] == 0);
}

TEST_CASE("pair_line: intersect general crossing", "[pair_line]") {
	ln a(pt(0, 0), pt(10, 10));
	ln b(pt(0, 10), pt(10, 0));
	auto r = intersect(a, b);
	REQUIRE(r[2] != 0);
	REQUIRE(r[0] * 2 == r[2]);
	REQUIRE(r[1] * 2 == r[2]);
}
