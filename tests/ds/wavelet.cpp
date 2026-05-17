#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "ds/wavelet.hpp"

TEST_CASE("bit_vector: set/get basic", "[wavelet][bit_vector]") {
	bit_vector bv(10);
	for (int i = 0; i < 10; ++ i) REQUIRE_FALSE(bv.get(i));
	bv.set(0);
	bv.set(3);
	bv.set(7);
	REQUIRE(bv.get(0));
	REQUIRE_FALSE(bv.get(1));
	REQUIRE(bv.get(3));
	REQUIRE(bv.get(7));
	bv.set(3, false);
	REQUIRE_FALSE(bv.get(3));
}

TEST_CASE("bit_vector: spans multiple words", "[wavelet][bit_vector]") {
	bit_vector bv(100);
	for (int i = 0; i < 100; i += 3) bv.set(i);
	for (int i = 0; i < 100; ++ i) {
		REQUIRE(bv[i] == ((i % 3) == 0));
	}
}

TEST_CASE("bit_vector: large size", "[wavelet][bit_vector]") {
	bit_vector bv(1000);
	std::mt19937 rng(1);
	std::vector<bool> ref(1000, false);
	for (int t = 0; t < 500; ++ t) {
		int i = int(rng() % 1000);
		bool v = bool(rng() & 1);
		bv.set(i, v);
		ref[i] = v;
	}
	for (int i = 0; i < 1000; ++ i) REQUIRE(bv.get(i) == ref[i]);
}

TEST_CASE("bit_vector: resize", "[wavelet][bit_vector]") {
	bit_vector bv;
	bv.resize(50);
	REQUIRE(bv.size() == 50);
	bv.set(25);
	REQUIRE(bv.get(25));
}
