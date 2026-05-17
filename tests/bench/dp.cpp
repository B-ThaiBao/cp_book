// Benchmarks for dp/line_multiset.hpp (Convex Hull Trick).
#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>

#include "dp/line_multiset.hpp"

TEST_CASE("line_multiset (CHT): 1e5 add + 1e5 max-query", "[!benchmark][cht]") {
	constexpr int N = 100'000;
	constexpr int Q = 100'000;
	std::mt19937_64 rng(0xCBA7ULL);
	std::vector<std::pair<long long, long long>> lines(N);
	for (auto& [m, b] : lines) {
		m = (long long)((rng() & 0x1ffff) - 0x10000);
		b = (long long)((rng() & 0x1ffff) - 0x10000);
	}
	std::vector<long long> qs(Q);
	for (auto& x : qs) x = (long long)((rng() & 0xfffff) - 0x80000);

	BENCHMARK("CHT add 1e5 lines + query 1e5 max") {
		line_multiset<line<long long>> hull;
		for (auto& [m, b] : lines) hull.add_line(m, b);
		long long acc = 0;
		for (auto& x : qs) {
			auto it = hull.lower_bound(x);
			if (it != hull.end()) acc += (*it)[0] * x + (*it)[1];
		}
		return acc;
	};
}
