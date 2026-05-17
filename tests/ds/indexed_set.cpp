#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "ds/indexed_set.hpp"

TEST_CASE("indexed_set: insert/find/erase", "[indexed_set]") {
	indexed_set s(100);
	REQUIRE_FALSE(s.find(5));
	s.insert(5);
	REQUIRE(s.find(5));
	s.insert(20);
	s.insert(42);
	REQUIRE(s.find(20));
	REQUIRE(s.find(42));
	REQUIRE_FALSE(s.find(43));
	s.erase(20);
	REQUIRE_FALSE(s.find(20));
}

TEST_CASE("indexed_set: find_next / find_prev", "[indexed_set]") {
	indexed_set s(1000);
	for (int x : {1, 5, 17, 100, 200, 500}) s.insert(x);
	REQUIRE(s.find_next(0) == 1);
	REQUIRE(s.find_next(1) == 1);
	REQUIRE(s.find_next(2) == 5);
	REQUIRE(s.find_next(18) == 100);
	REQUIRE(s.find_next(501) == 1000);
	REQUIRE(s.find_prev(999) == 500);
	REQUIRE(s.find_prev(150) == 100);
	REQUIRE(s.find_prev(5) == 5);
	REQUIRE(s.find_prev(0) == - 1);
}

TEST_CASE("indexed_set: empty boundary", "[indexed_set]") {
	indexed_set s(10);
	REQUIRE(s.find_next(0) == 10);
	REQUIRE(s.find_prev(9) == - 1);
}

TEST_CASE("indexed_set: for_each in range", "[indexed_set]") {
	indexed_set s(100);
	for (int x : {2, 10, 25, 26, 30, 70}) s.insert(x);
	std::vector<int> got;
	s.for_each(10, 30, [&](int x) { got.push_back(x); });
	REQUIRE(got == std::vector<int>({10, 25, 26}));
}

TEST_CASE("indexed_set: build from predicate", "[indexed_set]") {
	std::vector<int> vals = {1, 0, 1, 1, 0, 0, 1, 0, 0, 1};
	indexed_set s(int(vals.size()), [&](int i) { return vals[i]; });
	for (int i = 0; i < (int)vals.size(); ++ i) REQUIRE(s.find(i) == bool(vals[i]));
}

TEST_CASE("indexed_set: stress vs std::set", "[indexed_set][stress]") {
	std::mt19937 rng(99);
	int N = 2000;
	indexed_set is(N);
	std::set<int> ref;
	for (int t = 0; t < 8000; ++ t) {
		int op = int(rng() % 4);
		int x = int(rng() % N);
		if (op == 0) { is.insert(x); ref.insert(x); }
		else if (op == 1) { is.erase(x); ref.erase(x); }
		else if (op == 2) {
			auto it = ref.lower_bound(x);
			int expect = (it == ref.end()) ? N : *it;
			REQUIRE(is.find_next(x) == expect);
		} else {
			auto it = ref.upper_bound(x);
			int expect = (it == ref.begin()) ? - 1 : *std::prev(it);
			REQUIRE(is.find_prev(x) == expect);
		}
	}
}
