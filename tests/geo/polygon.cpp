#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "geo/pair_point.hpp"
#include "geo/polygon.hpp"

namespace {
using ptl = pair_point<long long>;
using ptd = pair_point<double>;
} // namespace

TEST_CASE("polygon: area square (doubled)", "[polygon]") {
	std::vector<ptl> sq = {{0, 0}, {2, 0}, {2, 2}, {0, 2}};
	REQUIRE(area<long long>(sq) == 8);
}

TEST_CASE("polygon: area clockwise gives negative", "[polygon]") {
	std::vector<ptl> sq = {{0, 0}, {0, 2}, {2, 2}, {2, 0}};
	REQUIRE(area<long long>(sq) == - 8);
}

TEST_CASE("polygon: area of triangle", "[polygon]") {
	std::vector<ptl> tri = {{0, 0}, {4, 0}, {0, 3}};
	REQUIRE(area<long long>(tri) == 12);
}

TEST_CASE("polygon: area of large rectangle", "[polygon]") {
	std::vector<ptl> r = {{0, 0}, {100, 0}, {100, 50}, {0, 50}};
	REQUIRE(area<long long>(r) == 100 * 50 * 2);
}

TEST_CASE("polygon: pentagon area", "[polygon]") {
	std::vector<ptl> p = {{0, 0}, {4, 0}, {5, 3}, {2, 5}, {- 1, 3}};
	long long got = area<long long>(p);
	REQUIRE(got > 0);
}

TEST_CASE("polygon: centroid square", "[polygon]") {
	std::vector<ptd> sq = {{0.0, 0.0}, {2.0, 0.0}, {2.0, 2.0}, {0.0, 2.0}};
	auto c = centroid<double>(sq);
	REQUIRE(std::abs(c.x - 1.0) < 1e-9);
	REQUIRE(std::abs(c.y - 1.0) < 1e-9);
}

TEST_CASE("polygon: centroid triangle", "[polygon]") {
	std::vector<ptd> tri = {{0.0, 0.0}, {6.0, 0.0}, {0.0, 6.0}};
	auto c = centroid<double>(tri);
	REQUIRE(std::abs(c.x - 2.0) < 1e-9);
	REQUIRE(std::abs(c.y - 2.0) < 1e-9);
}
