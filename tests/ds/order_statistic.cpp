#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "ds/order_statistic.hpp"

TEST_CASE("order_statistic_set: find_by_order / order_of_key", "[order_statistic]") {
	order_statistic_set<int> s;
	for (int x : {10, 5, 7, 3, 20, 15}) s.insert(x);
	REQUIRE(*s.find_by_order(0) == 3);
	REQUIRE(*s.find_by_order(2) == 7);
	REQUIRE(*s.find_by_order(5) == 20);
	REQUIRE(s.find_by_order(6) == s.end());
	REQUIRE(s.order_of_key(3) == 0);
	REQUIRE(s.order_of_key(11) == 4);
	REQUIRE(s.order_of_key(0) == 0);
	REQUIRE(s.order_of_key(100) == 6);
}

TEST_CASE("order_statistic_set: empty", "[order_statistic]") {
	order_statistic_set<int> s;
	REQUIRE(s.size() == 0);
	REQUIRE(s.order_of_key(5) == 0);
}

TEST_CASE("order_statistic_set: single element", "[order_statistic]") {
	order_statistic_set<int> s;
	s.insert(42);
	REQUIRE(*s.find_by_order(0) == 42);
	REQUIRE(s.order_of_key(42) == 0);
	REQUIRE(s.order_of_key(43) == 1);
}

TEST_CASE("order_statistic_set: stress vs sorted vector", "[order_statistic][stress]") {
	std::mt19937 rng(33);
	order_statistic_set<int> s;
	std::vector<int> ref;
	for (int t = 0; t < 500; ++ t) {
		int v = int(rng() % 1000);
		if (s.find(v) == s.end()) {
			s.insert(v);
			ref.insert(std::lower_bound(ref.begin(), ref.end(), v), v);
		}
	}
	for (int k = 0; k < (int)ref.size(); ++ k)
		REQUIRE(*s.find_by_order(k) == ref[k]);
	for (int v : {0, 50, 500, 999}) {
		int expect = int(std::lower_bound(ref.begin(), ref.end(), v) - ref.begin());
		REQUIRE((int)s.order_of_key(v) == expect);
	}
}

TEST_CASE("order_statistic_map: basic", "[order_statistic]") {
	order_statistic_map<int, std::string> m;
	m[10] = "ten";
	m[3] = "three";
	m[7] = "seven";
	REQUIRE(m.find_by_order(0)->first == 3);
	REQUIRE(m.find_by_order(1)->first == 7);
	REQUIRE(m.find_by_order(2)->first == 10);
	REQUIRE(m.order_of_key(7) == 1);
}
