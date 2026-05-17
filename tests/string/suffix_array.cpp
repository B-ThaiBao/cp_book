#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "ds/range_min_query.hpp"
#include "string/suffix_array.hpp"

namespace {
static std::vector<int> naive_sa(const std::vector<int>& s) {
	int N = int(s.size());
	std::vector<int> sa(N);
	std::iota(sa.begin(), sa.end(), 0);
	std::sort(sa.begin(), sa.end(), [&](const int& a, const int& b) {
		int x = a, y = b;
		while (x < N && y < N) {
			if (s[x] != s[y]) return s[x] < s[y];
			++ x; ++ y;
		}
		return x == N && y < N;
	});
	return sa;
}

static std::vector<int> naive_lcp(const std::vector<int>& s, const std::vector<int>& sa) {
	int N = int(s.size());
	std::vector<int> lcp(N - 1, 0);
	for (int i = 0; i + 1 < N; ++ i) {
		int a = sa[i], b = sa[i + 1], k = 0;
		while (a + k < N && b + k < N && s[a + k] == s[b + k]) ++ k;
		lcp[i] = k;
	}
	return lcp;
}
} // namespace

TEST_CASE("suffix_array: example 'banana'", "[sa]") {
	std::vector<int> s = {'b', 'a', 'n', 'a', 'n', 'a'};
	suffix_array sa(s);
	REQUIRE(std::vector<int>(sa.begin(), sa.end()) == std::vector<int>({5, 3, 1, 0, 4, 2}));
}

TEST_CASE("suffix_array: single char", "[sa]") {
	std::vector<int> s = {7};
	suffix_array sa(s);
	REQUIRE(std::vector<int>(sa.begin(), sa.end()) == std::vector<int>({0}));
}

TEST_CASE("suffix_array: all same", "[sa]") {
	std::vector<int> s = {5, 5, 5, 5};
	suffix_array sa(s);
	REQUIRE(std::vector<int>(sa.begin(), sa.end()) == std::vector<int>({3, 2, 1, 0}));
}

TEST_CASE("suffix_array: stress vs naive (alpha 4)", "[sa][stress]") {
	std::mt19937 rng(123);
	for (int t = 0; t < 40; ++ t) {
		int N = 1 + int(rng() % 30u);
		std::vector<int> s(N);
		for (int i = 0; i < N; ++ i) s[i] = 1 + int(rng() % 4u);
		suffix_array sa(s);
		REQUIRE(std::vector<int>(sa.begin(), sa.end()) == naive_sa(s));
	}
}

TEST_CASE("suffix_array: stress vs naive (binary)", "[sa][stress]") {
	std::mt19937 rng(321);
	for (int t = 0; t < 40; ++ t) {
		int N = 1 + int(rng() % 50u);
		std::vector<int> s(N);
		for (int i = 0; i < N; ++ i) s[i] = 1 + int(rng() % 2u);
		suffix_array sa(s);
		REQUIRE(std::vector<int>(sa.begin(), sa.end()) == naive_sa(s));
	}
}

TEST_CASE("suffix_lcp_array: stress vs naive", "[sa][lcp][stress]") {
	std::mt19937 rng(456);
	for (int t = 0; t < 25; ++ t) {
		int N = 2 + int(rng() % 20u);
		std::vector<int> s(N);
		for (int i = 0; i < N; ++ i) s[i] = 1 + int(rng() % 3u);
		suffix_array sa(s);
		suffix_lcp_array lcp;
		lcp.build(s, sa, true);
		auto ref = naive_lcp(s, std::vector<int>(sa.begin(), sa.end()));
		REQUIRE(lcp.lcp == std::vector<int32_t>(ref.begin(), ref.end()));
		for (int a = 0; a < N; ++ a)
			if (lcp.rank[a] != N - 1)
				REQUIRE(lcp.find_lcp(a) == lcp.lcp[lcp.rank[a]]);
	}
}

TEST_CASE("suffix_lcp_array: count_diff small", "[sa][lcp]") {
	{
		std::vector<int> s = {'a', 'b', 'a', 'b'};
		suffix_array sa(s);
		suffix_lcp_array lcp;
		lcp.build(s, sa, false);
		REQUIRE(lcp.count_diff() == 7);
	}
	{
		std::vector<int> s = {'a'};
		suffix_array sa(s);
		suffix_lcp_array lcp;
		lcp.build(s, sa, false);
		REQUIRE(lcp.count_diff() == 1);
	}
	{
		std::vector<int> s = {'a', 'a', 'a'};
		suffix_array sa(s);
		suffix_lcp_array lcp;
		lcp.build(s, sa, false);
		REQUIRE(lcp.count_diff() == 3);
	}
}

TEST_CASE("suffix_lcp_array: count_diff stress", "[sa][lcp][stress]") {
	std::mt19937 rng(2024);
	for (int t = 0; t < 15; ++ t) {
		int N = 2 + int(rng() % 12u);
		std::vector<int> s(N);
		for (int i = 0; i < N; ++ i) s[i] = 1 + int(rng() % 3u);
		suffix_array sa(s);
		suffix_lcp_array lcp;
		lcp.build(s, sa, false);
		std::set<std::vector<int>> distinct;
		for (int i = 0; i < N; ++ i)
			for (int j = i + 1; j <= N; ++ j)
				distinct.insert(std::vector<int>(s.begin() + i, s.begin() + j));
		REQUIRE(lcp.count_diff() == (long long)distinct.size());
	}
}
