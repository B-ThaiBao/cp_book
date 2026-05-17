#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "geo/pair_point.hpp"
#include "geo/pair_line.hpp"
#include "geo/half_plane_intersect.hpp"

namespace {
using pt = pair_point<long long>;
using ln = pair_line<long long>;
} // namespace

TEST_CASE("half_plane_intersect: bounded square", "[hpi]") {
	std::vector<ln> pls = {
		ln(pt(0, 0), pt(4, 0)),
		ln(pt(4, 0), pt(4, 4)),
		ln(pt(4, 4), pt(0, 4)),
		ln(pt(0, 4), pt(0, 0)),
	};
	std::vector<int> order = {0, 1, 2, 3};
	auto res = half_plane_intersect(pls, order);
	REQUIRE(res.size() == 4);
}

TEST_CASE("half_plane_intersect: triangle", "[hpi]") {
	std::vector<ln> pls = {
		ln(pt(0, 0), pt(6, 0)),
		ln(pt(6, 0), pt(0, 6)),
		ln(pt(0, 6), pt(0, 0)),
	};
	std::vector<int> order = {0, 1, 2};
	auto res = half_plane_intersect(pls, order);
	REQUIRE(res.size() == 3);
}

TEST_CASE("half_plane_intersect: unbounded after dropping one edge", "[hpi]") {
	std::vector<ln> pls = {
		ln(pt(0, 0), pt(4, 0)),
		ln(pt(4, 0), pt(4, 4)),
		ln(pt(4, 4), pt(0, 4)),
	};
	std::vector<int> order = {0, 1, 2};
	auto res = half_plane_intersect(pls, order);
	REQUIRE(res.size() <= 3);
}

TEST_CASE("half_plane_intersect: larger pentagon", "[hpi]") {
	std::vector<ln> pls = {
		ln(pt(0, 0), pt(10, 0)),
		ln(pt(10, 0), pt(12, 5)),
		ln(pt(12, 5), pt(5, 10)),
		ln(pt(5, 10), pt(- 2, 5)),
		ln(pt(- 2, 5), pt(0, 0)),
	};
	std::vector<int> order = {0, 1, 2, 3, 4};
	auto res = half_plane_intersect(pls, order);
	REQUIRE(res.size() == 5);
}
