// Benchmarks for heap/* and BST/* heavy structures (radix_heap vs priority_queue,
// treap, splay, link_cut_tree).
#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>

#include "ds/heap/radix_heap.hpp"

TEST_CASE("radix_heap vs std::priority_queue (1e6 monotone pushes + 1e6 pops)", "[!benchmark][heap]") {
	constexpr int N = 1'000'000;
	std::mt19937_64 rng(2025);
	std::vector<int32_t> vals(N);
	for (auto& v : vals) v = int32_t(rng() & 0x0fffffffu);

	BENCHMARK("radix_heap<int32_t,int>: 1e6 push + 1e6 pop") {
		radix_heap<int32_t, int> h;
		for (int i = 0; i < N; ++i) h.emplace(vals[i], i);
		uint64_t acc = 0;
		while (!h.empty()) {
			auto p = h.pop();
			acc += uint32_t(p.first);
		}
		return acc;
	};

	BENCHMARK("std::priority_queue min-heap: 1e6 push + 1e6 pop") {
		std::priority_queue<std::pair<int32_t, int>,
			std::vector<std::pair<int32_t, int>>, std::greater<>> pq;
		for (int i = 0; i < N; ++i) pq.push({vals[i], i});
		uint64_t acc = 0;
		while (!pq.empty()) {
			acc += pq.top().first;
			pq.pop();
		}
		return acc;
	};
}
