#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "ds/hash_map.hpp"

TEST_CASE("hash_map: insert/find/erase basic", "[hash_map]") {
	hash_map<int, std::string> m;
	m[1] = "one";
	m[2] = "two";
	m[42] = "answer";
	REQUIRE(m.size() == 3);
	REQUIRE(m.find(1) != m.end());
	REQUIRE(m.find(99) == m.end());
	REQUIRE(m[42] == "answer");
	m.erase(2);
	REQUIRE(m.find(2) == m.end());
	REQUIRE(m.size() == 2);
}

TEST_CASE("hash_map: overwrite value", "[hash_map]") {
	hash_map<int, int> m;
	m[1] = 10;
	m[1] = 20;
	REQUIRE(m.size() == 1);
	REQUIRE(m[1] == 20);
}

TEST_CASE("hash_set: basic", "[hash_set]") {
	hash_set<int> s;
	for (int i = 0; i < 100; ++ i) s.insert(i * 7);
	REQUIRE(s.size() == 100);
	for (int i = 0; i < 100; ++ i) REQUIRE(s.find(i * 7) != s.end());
	REQUIRE(s.find(1) == s.end());
}

TEST_CASE("hash_set: duplicate insertion", "[hash_set]") {
	hash_set<int> s;
	s.insert(5);
	s.insert(5);
	s.insert(5);
	REQUIRE(s.size() == 1);
}

TEST_CASE("hash_map: stress vs std::map", "[hash_map][stress]") {
	std::mt19937 rng(2024);
	hash_map<long long, int> hm;
	std::map<long long, int> ref;
	for (int t = 0; t < 3000; ++ t) {
		long long k = (long long)(rng() % 200);
		int op = rng() % 3;
		if (op == 0) {
			int v = int(rng() % 1000);
			hm[k] = v;
			ref[k] = v;
		} else if (op == 1) {
			hm.erase(k);
			ref.erase(k);
		} else {
			auto it = hm.find(k);
			auto rit = ref.find(k);
			REQUIRE((it == hm.end()) == (rit == ref.end()));
			if (rit != ref.end()) REQUIRE(it->second == rit->second);
		}
	}
	REQUIRE(hm.size() == ref.size());
}
