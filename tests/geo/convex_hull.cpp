#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "geo/pair_point.hpp"
#include "geo/convex_hull.hpp"

namespace {
using pt = pair_point<long long>;

static std::vector<int> sort_order(const std::vector<pt>& pts) {
	std::vector<int> order(pts.size());
	std::iota(order.begin(), order.end(), 0);
	std::sort(order.begin(), order.end(), [&](const int& a, const int& b) { return pts[a] < pts[b]; });
	return order;
}
} // namespace

TEST_CASE("convex_hull: square with interior point", "[convex_hull]") {
	std::vector<pt> pts = {{0, 0}, {2, 0}, {2, 2}, {0, 2}, {1, 1}};
	auto order = sort_order(pts);
	auto ch = convex_hull(pts, order);
	REQUIRE(ch[0].size() == 4);
	REQUIRE(ch[1][4] == 6);
}

TEST_CASE("convex_hull: collinear on edge marked N", "[convex_hull]") {
	std::vector<pt> pts = {{0, 0}, {1, 0}, {2, 0}, {2, 2}, {0, 2}};
	auto order = sort_order(pts);
	auto ch = convex_hull(pts, order);
	REQUIRE(ch[1][1] == 5);
}

TEST_CASE("convex_hull: triangle equals hull", "[convex_hull]") {
	std::vector<pt> pts = {{0, 0}, {3, 0}, {0, 4}};
	auto order = sort_order(pts);
	auto ch = convex_hull(pts, order);
	REQUIRE(ch[0].size() == 3);
}

TEST_CASE("convex_hull: many points hull size", "[convex_hull]") {
	std::vector<pt> pts = {{0, 0}, {10, 0}, {10, 10}, {0, 10}, {5, 5}, {3, 3}, {7, 2}, {2, 7}};
	auto order = sort_order(pts);
	auto ch = convex_hull(pts, order);
	REQUIRE(ch[0].size() == 4);
}

TEST_CASE("convex_hull: single point", "[convex_hull]") {
	std::vector<pt> pts = {{5, 5}};
	auto order = sort_order(pts);
	auto ch = convex_hull(pts, order);
	REQUIRE(ch[0].size() == 1);
}

TEST_CASE("convex_hull: two points", "[convex_hull]") {
	std::vector<pt> pts = {{0, 0}, {5, 5}};
	auto order = sort_order(pts);
	auto ch = convex_hull(pts, order);
	REQUIRE(ch[0].size() == 2);
}
