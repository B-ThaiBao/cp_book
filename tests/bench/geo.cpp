// Benchmarks for geo/convex_hull.hpp.
#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>

#include "geo/pair_point.hpp"
#include "geo/convex_hull.hpp"

using pt = pair_point<long long>;

TEST_CASE("convex_hull: 1e6 random points in a square", "[!benchmark][convex_hull]") {
	constexpr int N = 1'000'000;
	std::mt19937_64 rng(31337);
	std::vector<pt> pts(N);
	for (auto& p : pts) p = pt((long long)(rng() & 0xfffff), (long long)(rng() & 0xfffff));

	BENCHMARK("convex_hull(N=1e6)") {
		std::vector<int> order(N); std::iota(order.begin(), order.end(), 0);
		std::sort(order.begin(), order.end(), [&](int a, int b) { return pts[a] < pts[b]; });
		auto ch = convex_hull(pts, order);
		return ch[0].size();
	};
}

TEST_CASE("convex_hull: 1e5 points on a circle (worst case = all on hull)", "[!benchmark][convex_hull]") {
	constexpr int N = 100'000;
	std::vector<pt> pts(N);
	const long long R = 1'000'000'000LL;
	std::mt19937_64 rng(1);
	std::vector<double> angles(N);
	for (auto& a : angles) a = (double)rng() / (double)UINT64_MAX * 2 * acos(-1.0);
	std::sort(angles.begin(), angles.end());
	for (int i = 0; i < N; ++i) {
		pts[i] = pt((long long)(R * cos(angles[i])), (long long)(R * sin(angles[i])));
	}

	BENCHMARK("convex_hull(N=1e5, all on hull)") {
		std::vector<int> order(N); std::iota(order.begin(), order.end(), 0);
		std::sort(order.begin(), order.end(), [&](int a, int b) { return pts[a] < pts[b]; });
		auto ch = convex_hull(pts, order);
		return ch[0].size();
	};
}
