#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "string/prefix_function.hpp"

namespace {
static std::vector<int> naive_pi(const std::string& s) {
	int N = int(s.size());
	std::vector<int> p(N, 0);
	for (int i = 1; i < N; ++ i) {
		int j = p[i - 1];
		while (j > 0 && s[i] != s[j]) j = p[j - 1];
		if (s[i] == s[j]) ++ j;
		p[i] = j;
	}
	return p;
}

static std::vector<int> naive_find(const std::string& w, const std::string& s) {
	std::vector<int> r;
	if (s.empty() || s.size() > w.size()) return r;
	for (size_t i = 0; i + s.size() <= w.size(); ++ i)
		if (w.compare(i, s.size(), s) == 0) r.push_back(int(i));
	return r;
}

static std::string rand_str(std::mt19937& rng, const int& N, const int& alpha) {
	std::string s; s.reserve(N);
	for (int i = 0; i < N; ++ i) s += char('a' + int(rng() % uint32_t(alpha)));
	return s;
}
} // namespace

TEST_CASE("prefix_function: basic example", "[kmp]") {
	prefix_function p(std::string("abacaba"));
	REQUIRE(std::vector<int>(p.begin(), p.end()) == std::vector<int>({0, 0, 1, 0, 1, 2, 3}));
}

TEST_CASE("prefix_function: single char", "[kmp]") {
	prefix_function p(std::string("a"));
	REQUIRE(std::vector<int>(p.begin(), p.end()) == std::vector<int>({0}));
}

TEST_CASE("prefix_function: all same", "[kmp]") {
	prefix_function p(std::string("aaaaa"));
	REQUIRE(std::vector<int>(p.begin(), p.end()) == std::vector<int>({0, 1, 2, 3, 4}));
}

TEST_CASE("prefix_function: stress vs naive binary", "[kmp][stress]") {
	std::mt19937 rng(9);
	for (int t = 0; t < 50; ++ t) {
		int N = 1 + int(rng() % 30u);
		std::string s = rand_str(rng, N, 2);
		prefix_function p(s);
		REQUIRE(std::vector<int>(p.begin(), p.end()) == naive_pi(s));
	}
}

TEST_CASE("prefix_function: stress vs naive alpha 3", "[kmp][stress]") {
	std::mt19937 rng(21);
	for (int t = 0; t < 50; ++ t) {
		int N = 1 + int(rng() % 50u);
		std::string s = rand_str(rng, N, 3);
		prefix_function p(s);
		REQUIRE(std::vector<int>(p.begin(), p.end()) == naive_pi(s));
	}
}

TEST_CASE("find_pos: finds all occurrences", "[kmp]") {
	std::string w = "ababcababab", s = "ab";
	prefix_function pi(s);
	REQUIRE(find_pos(w, s, pi) == naive_find(w, s));
}

TEST_CASE("find_pos: no occurrence", "[kmp]") {
	std::string w = "abcdef", s = "xyz";
	prefix_function pi(s);
	REQUIRE(find_pos(w, s, pi).empty());
}

TEST_CASE("find_pos: pattern equals text", "[kmp]") {
	std::string w = "hello", s = "hello";
	prefix_function pi(s);
	REQUIRE(find_pos(w, s, pi) == std::vector<int>({0}));
}

TEST_CASE("find_pos: stress vs naive", "[kmp][stress]") {
	std::mt19937 rng(11);
	for (int t = 0; t < 50; ++ t) {
		int N = 5 + int(rng() % 40u);
		int M = 1 + int(rng() % 5u);
		std::string w = rand_str(rng, N, 2);
		std::string s = rand_str(rng, M, 2);
		prefix_function pi(s);
		REQUIRE(find_pos(w, s, pi) == naive_find(w, s));
	}
}

TEST_CASE("compress_prefix periodicity", "[kmp]") {
	REQUIRE(compress_prefix(prefix_function(std::string("abcabcabc"))) == 3);
	REQUIRE(compress_prefix(prefix_function(std::string("abcd"))) == 4);
	REQUIRE(compress_prefix(prefix_function(std::string("aaaaa"))) == 1);
	REQUIRE(compress_prefix(prefix_function(std::string("abab"))) == 2);
	REQUIRE(compress_prefix(prefix_function(std::string("aabaab"))) == 3);
}
