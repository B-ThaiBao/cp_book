#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "geo/pair_point.hpp"
#include "geo/convex_inclusion.hpp"

namespace {
using pt = pair_point<long long>;
} // namespace

TEST_CASE("convex_inclusion: inside square", "[convex_inclusion]") {
	std::vector<pt> sq = {{0, 0}, {4, 0}, {4, 4}, {0, 4}};
	REQUIRE(convex_inclusion(pt(2, 2), 0, sq)[0] == - 1);
	REQUIRE(convex_inclusion(pt(1, 3), 0, sq)[0] == - 1);
	REQUIRE(convex_inclusion(pt(3, 1), 0, sq)[0] == - 1);
}

TEST_CASE("convex_inclusion: outside square", "[convex_inclusion]") {
	std::vector<pt> sq = {{0, 0}, {4, 0}, {4, 4}, {0, 4}};
	REQUIRE(convex_inclusion(pt(5, 5), 0, sq)[0] == 1);
	REQUIRE(convex_inclusion(pt(- 1, 2), 0, sq)[0] == 1);
	REQUIRE(convex_inclusion(pt(2, 10), 0, sq)[0] == 1);
}

TEST_CASE("convex_inclusion: on edge", "[convex_inclusion]") {
	std::vector<pt> sq = {{0, 0}, {4, 0}, {4, 4}, {0, 4}};
	REQUIRE(convex_inclusion(pt(2, 0), 0, sq)[0] == 0);
	REQUIRE(convex_inclusion(pt(4, 2), 0, sq)[0] == 0);
}

TEST_CASE("convex_inclusion: triangle inside/outside", "[convex_inclusion]") {
	std::vector<pt> tri = {{0, 0}, {6, 0}, {0, 6}};
	REQUIRE(convex_inclusion(pt(1, 1), 0, tri)[0] == - 1);
	REQUIRE(convex_inclusion(pt(10, 10), 0, tri)[0] == 1);
	REQUIRE(convex_inclusion(pt(- 1, 3), 0, tri)[0] == 1);
}

TEST_CASE("convex_inclusion: pentagon center inside", "[convex_inclusion]") {
	std::vector<pt> pen = {{0, 0}, {10, 0}, {12, 5}, {5, 10}, {- 2, 5}};
	REQUIRE(convex_inclusion(pt(5, 4), 0, pen)[0] == - 1);
	REQUIRE(convex_inclusion(pt(100, 100), 0, pen)[0] == 1);
}
