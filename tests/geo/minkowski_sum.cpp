#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "geo/pair_point.hpp"
#include "geo/minkowski_sum.hpp"

namespace {
using pt = pair_point<long long>;
} // namespace

TEST_CASE("minkowski_sum: square + square", "[minkowski]") {
	std::vector<pt> A = {{0, 0}, {2, 0}, {2, 2}, {0, 2}};
	std::vector<pt> B = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
	auto C = minkowski_sum(A, B, 0, 0);
	REQUIRE(C.size() == 4);
	std::set<std::pair<long long, long long>> got;
	for (auto& p : C) got.insert({p.x, p.y});
	std::set<std::pair<long long, long long>> want = {{0, 0}, {3, 0}, {3, 3}, {0, 3}};
	REQUIRE(got == want);
}

TEST_CASE("minkowski_sum: triangle + triangle", "[minkowski]") {
	std::vector<pt> A = {{0, 0}, {2, 0}, {0, 2}};
	std::vector<pt> B = {{0, 0}, {1, 0}, {0, 1}};
	auto C = minkowski_sum(A, B, 0, 0);
	std::set<std::pair<long long, long long>> got;
	for (auto& p : C) got.insert({p.x, p.y});
	REQUIRE(got.count({0, 0}) == 1);
	REQUIRE(got.count({3, 0}) == 1);
	REQUIRE(got.count({0, 3}) == 1);
}

TEST_CASE("minkowski_sum: square + triangle", "[minkowski]") {
	std::vector<pt> A = {{0, 0}, {2, 0}, {2, 2}, {0, 2}};
	std::vector<pt> B = {{0, 0}, {1, 0}, {0, 1}};
	auto C = minkowski_sum(A, B, 0, 0);
	REQUIRE(C.size() >= 4);
	std::set<std::pair<long long, long long>> got;
	for (auto& p : C) got.insert({p.x, p.y});
	REQUIRE(got.count({0, 0}) == 1);
	REQUIRE(got.count({3, 0}) == 1);
}

TEST_CASE("minkowski_sum: single point + square", "[minkowski]") {
	std::vector<pt> A = {{5, 5}};
	std::vector<pt> B = {{0, 0}, {2, 0}, {2, 2}, {0, 2}};
	auto C = minkowski_sum(A, B, 0, 0);
	REQUIRE(C.size() == 4);
	std::set<std::pair<long long, long long>> got;
	for (auto& p : C) got.insert({p.x, p.y});
	REQUIRE(got.count({5, 5}) == 1);
	REQUIRE(got.count({7, 5}) == 1);
	REQUIRE(got.count({7, 7}) == 1);
	REQUIRE(got.count({5, 7}) == 1);
}
