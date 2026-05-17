#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "string/manacher.hpp"

namespace {
static bool is_pal_naive(const std::string& s, const int& l, const int& r) {
	if (l >= r) return false;
	int a = l, b = r;
	while (a < b - 1) { if (s[a] != s[b - 1]) return false; ++ a; -- b; }
	return true;
}

static std::string rand_str(std::mt19937& rng, const int& N, const int& alpha) {
	std::string s; s.reserve(N);
	for (int i = 0; i < N; ++ i) s += char('a' + int(rng() % uint32_t(alpha)));
	return s;
}
} // namespace

TEST_CASE("manacher: documented example", "[manacher]") {
	manacher m;
	m.build(std::string("abaa"));
	REQUIRE(int(m.size()) == 7);
	REQUIRE(std::vector<int>(m.begin(), m.end()) == std::vector<int>({0, 0, 1, 0, 0, 1, 0}));
}

TEST_CASE("manacher: palindrome of full length", "[manacher]") {
	std::string s = "abcba";
	manacher m(s);
	REQUIRE(m.is_palindrome(0, 5));
	REQUIRE_FALSE(m.is_palindrome(0, 4));
}

TEST_CASE("manacher: single char palindromes", "[manacher]") {
	manacher m(std::string("xyz"));
	REQUIRE(m.is_palindrome(0, 1));
	REQUIRE(m.is_palindrome(1, 2));
	REQUIRE(m.is_palindrome(2, 3));
	REQUIRE_FALSE(m.is_palindrome(0, 2));
}

TEST_CASE("manacher: all-same string", "[manacher]") {
	manacher m(std::string("aaaaaa"));
	for (int l = 0; l < 6; ++ l)
		for (int r = l + 1; r <= 6; ++ r)
			REQUIRE(m.is_palindrome(l, r));
}

TEST_CASE("manacher: stress vs naive (binary)", "[manacher][stress]") {
	std::mt19937 rng(7);
	for (int trial = 0; trial < 40; ++ trial) {
		int N = 1 + int(rng() % 30u);
		std::string s = rand_str(rng, N, 2);
		manacher m(s);
		for (int l = 0; l <= N; ++ l)
			for (int r = l + 1; r <= N; ++ r)
				REQUIRE(m.is_palindrome(l, r) == is_pal_naive(s, l, r));
	}
}

TEST_CASE("manacher: stress vs naive (alphabet 3)", "[manacher][stress]") {
	std::mt19937 rng(17);
	for (int trial = 0; trial < 30; ++ trial) {
		int N = 1 + int(rng() % 40u);
		std::string s = rand_str(rng, N, 3);
		manacher m(s);
		for (int l = 0; l <= N; ++ l)
			for (int r = l + 1; r <= N; ++ r)
				REQUIRE(m.is_palindrome(l, r) == is_pal_naive(s, l, r));
	}
}

TEST_CASE("manacher: even-length palindrome", "[manacher]") {
	manacher m(std::string("abba"));
	REQUIRE(m.is_palindrome(0, 4));
	REQUIRE(m.is_palindrome(1, 3));
}

TEST_CASE("manacher: longest palindrome substring via scan", "[manacher][stress]") {
	std::mt19937 rng(33);
	for (int trial = 0; trial < 20; ++ trial) {
		int N = 1 + int(rng() % 25u);
		std::string s = rand_str(rng, N, 3);
		manacher m(s);
		int best_m = 0, best_n = 0;
		for (int l = 0; l <= N; ++ l)
			for (int r = l + 1; r <= N; ++ r) {
				if (m.is_palindrome(l, r)) best_m = std::max(best_m, r - l);
				if (is_pal_naive(s, l, r)) best_n = std::max(best_n, r - l);
			}
		REQUIRE(best_m == best_n);
	}
}
