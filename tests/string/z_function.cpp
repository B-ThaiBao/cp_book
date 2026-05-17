#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "string/z_function.hpp"

namespace {
static std::vector<int> naive_z(const std::string& s) {
	int N = int(s.size());
	std::vector<int> z(N, 0);
	for (int i = 1; i < N; ++ i) {
		while (i + z[i] < N && s[z[i]] == s[i + z[i]]) ++ z[i];
	}
	return z;
}
} // namespace

TEST_CASE("z_function: all same", "[z_function]") {
	z_function z(std::string("aaaaa"));
	REQUIRE(std::vector<int>(z.begin(), z.end()) == std::vector<int>({0, 4, 3, 2, 1}));
}

TEST_CASE("z_function: abacaba", "[z_function]") {
	z_function z(std::string("abacaba"));
	REQUIRE(std::vector<int>(z.begin(), z.end()) == std::vector<int>({0, 0, 1, 0, 3, 0, 1}));
}

TEST_CASE("z_function: single char", "[z_function]") {
	z_function z(std::string("a"));
	REQUIRE(std::vector<int>(z.begin(), z.end()) == std::vector<int>({0}));
}

TEST_CASE("z_function: distinct chars", "[z_function]") {
	z_function z(std::string("abcde"));
	REQUIRE(std::vector<int>(z.begin(), z.end()) == std::vector<int>({0, 0, 0, 0, 0}));
}

TEST_CASE("z_function: stress vs naive (binary)", "[z_function][stress]") {
	std::mt19937 rng(3);
	for (int t = 0; t < 60; ++ t) {
		int N = 1 + int(rng() % 40u);
		std::string s;
		for (int i = 0; i < N; ++ i) s += char('a' + int(rng() % 2u));
		z_function z(s);
		REQUIRE(std::vector<int>(z.begin(), z.end()) == naive_z(s));
	}
}

TEST_CASE("z_function: stress vs naive (alpha 4)", "[z_function][stress]") {
	std::mt19937 rng(13);
	for (int t = 0; t < 50; ++ t) {
		int N = 1 + int(rng() % 50u);
		std::string s;
		for (int i = 0; i < N; ++ i) s += char('a' + int(rng() % 4u));
		z_function z(s);
		REQUIRE(std::vector<int>(z.begin(), z.end()) == naive_z(s));
	}
}

TEST_CASE("z_function: compress periodicity", "[z_function]") {
	REQUIRE(z_function(std::string("abcabcabc")).compress() == 3);
	REQUIRE(z_function(std::string("aaaaa")).compress() == 1);
	REQUIRE(z_function(std::string("abcde")).compress() == 5);
	REQUIRE(z_function(std::string("abab")).compress() == 2);
}
