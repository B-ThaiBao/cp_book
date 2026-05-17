#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "ds/heap/radix_heap.hpp"

TEST_CASE("radix_heap: push + pop ascending", "[radix_heap]") {
	radix_heap<int32_t, int> h;
	h.emplace(0, 100);
	h.emplace(3, 13);
	h.emplace(5, 25);
	h.emplace(10, 50);
	auto p = h.pop(); REQUIRE(p.first == 0); REQUIRE(p.second == 100);
	p = h.pop(); REQUIRE(p.first == 3);
	p = h.pop(); REQUIRE(p.first == 5);
	p = h.pop(); REQUIRE(p.first == 10);
	REQUIRE(h.empty());
}

TEST_CASE("radix_heap: top returns min", "[radix_heap]") {
	radix_heap<int32_t, int> h;
	h.emplace(7, 0);
	h.emplace(8, 1);
	REQUIRE(h.top().first == 7);
	h.pop();
	REQUIRE(h.top().first == 8);
}

TEST_CASE("radix_heap: int64 keys", "[radix_heap]") {
	radix_heap<int64_t, int> h;
	for (int64_t k : {(int64_t)1, (int64_t)1000, (int64_t)1000000LL, (int64_t)10000000000LL})
		h.emplace(k, int(k % 10));
	int64_t last = - 1;
	while (!h.empty()) {
		auto p = h.pop();
		REQUIRE(p.first >= last);
		last = p.first;
	}
}

TEST_CASE("radix_heap: monotone stress vs priority_queue", "[radix_heap][stress]") {
	std::mt19937 rng(2024);
	radix_heap<int32_t, int> h;
	std::priority_queue<int, std::vector<int>, std::greater<int>> ref;
	int32_t last_popped = 0;
	for (int t = 0; t < 1500; ++ t) {
		int op = int(rng() % 2);
		if (op == 0 || ref.empty()) {
			int32_t k = last_popped + int32_t(rng() % 100);
			h.emplace(k, 0);
			ref.push(k);
		} else {
			int expect = ref.top(); ref.pop();
			auto p = h.pop();
			REQUIRE(p.first == expect);
			last_popped = p.first;
		}
	}
	while (!ref.empty()) {
		int expect = ref.top(); ref.pop();
		REQUIRE(h.pop().first == expect);
	}
}
