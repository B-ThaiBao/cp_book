#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "ds/van_emde_boas_tree.hpp"

TEST_CASE("van_emde_boas_tree: insert/find_next/find_prev", "[vEB]") {
	van_emde_boas_tree<8> s;
	s.insert(10);
	s.insert(50);
	s.insert(100);
	REQUIRE(s.find_next(0) == 10);
	REQUIRE(s.find_next(10) == 10);
	REQUIRE(s.find_next(11) == 50);
	REQUIRE(s.find_next(101) == 256);
	REQUIRE(s.find_prev(255) == 100);
	REQUIRE(s.find_prev(50) == 50);
	REQUIRE(s.find_prev(49) == 10);
	REQUIRE(s.find_prev(9) == - 1);
}

TEST_CASE("van_emde_boas_tree: erase", "[vEB]") {
	van_emde_boas_tree<8> s;
	for (int x : {5, 10, 50, 100, 200}) s.insert(x);
	s.erase(50);
	REQUIRE(s.find_next(11) == 100);
	s.erase(5);
	REQUIRE(s.find_next(0) == 10);
	s.erase(200);
	REQUIRE(s.find_prev(255) == 100);
}

TEST_CASE("van_emde_boas_tree: empty/min/max", "[vEB]") {
	van_emde_boas_tree<8> s;
	REQUIRE(s.empty());
	s.insert(42);
	REQUIRE_FALSE(s.empty());
	REQUIRE(s.min() == 42);
	REQUIRE(s.max() == 42);
	s.erase(42);
	REQUIRE(s.empty());
}

TEST_CASE("van_emde_boas_tree: stress vs std::set", "[vEB][stress]") {
	std::mt19937 rng(33);
	van_emde_boas_tree<10> s;
	std::set<int> ref;
	int LIM = 1 << 10;
	for (int t = 0; t < 5000; ++ t) {
		int op = int(rng() % 4);
		int x = int(rng() % unsigned(LIM));
		if (op == 0) { s.insert(x); ref.insert(x); }
		else if (op == 1) { s.erase(x); ref.erase(x); }
		else if (op == 2) {
			auto it = ref.lower_bound(x);
			int expect = (it == ref.end()) ? LIM : *it;
			REQUIRE(s.find_next(x) == expect);
		} else {
			auto it = ref.upper_bound(x);
			int expect = (it == ref.begin()) ? - 1 : *std::prev(it);
			REQUIRE(s.find_prev(x) == expect);
		}
	}
}
