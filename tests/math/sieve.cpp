#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "math/sieve.hpp"

namespace {
static std::vector<bool> trial_sieve(const int& N) {
	std::vector<bool> p(N, true);
	if (N > 0) p[0] = false;
	if (N > 1) p[1] = false;
	for (int i = 2; i < N; ++ i) {
		if (!p[i]) continue;
		for (int j = 2 * i; j < N; j += i) p[j] = false;
	}
	return p;
}
} // namespace

TEST_CASE("eratos_sieve: basic", "[sieve]") {
	int N = 500;
	REQUIRE(eratos_sieve(N) == trial_sieve(N));
}

TEST_CASE("eratos_sieve: tiny", "[sieve]") {
	REQUIRE(eratos_sieve(2) == std::vector<bool>({false, false}));
	REQUIRE(eratos_sieve(3) == std::vector<bool>({false, false, true}));
}

TEST_CASE("block_sieve: prime list", "[sieve]") {
	int K = 500;
	auto primes = block_sieve(K);
	auto ref = trial_sieve(K);
	std::vector<int> ref_list;
	for (int i = 2; i < K; ++ i) if (ref[i]) ref_list.push_back(i);
	REQUIRE(primes == ref_list);
}

TEST_CASE("block_sieve: small K", "[sieve]") {
	REQUIRE(block_sieve(10) == std::vector<int>({2, 3, 5, 7}));
	REQUIRE(block_sieve(30) == std::vector<int>({2, 3, 5, 7, 11, 13, 17, 19, 23, 29}));
}

TEST_CASE("linear_sieve: smallest prime factor", "[sieve][linear]") {
	int N = 500;
	linear_sieve ls(N);
	for (int x = 2; x < N; ++ x) {
		int spf = - 1;
		for (int p = 2; p * p <= x; ++ p) if (x % p == 0) { spf = p; break; }
		if (spf == - 1) spf = x;
		REQUIRE(ls[x] == spf);
	}
}

TEST_CASE("linear_sieve: primes list matches eratos", "[sieve][linear]") {
	int N = 500;
	linear_sieve ls(N);
	std::set<int> got(ls.primes.begin(), ls.primes.end());
	auto ref = trial_sieve(N);
	for (int i = 2; i < N; ++ i)
		REQUIRE(got.count(i) == size_t(ref[i] ? 1 : 0));
}

TEST_CASE("linear_sieve: prime count", "[sieve][linear]") {
	linear_sieve ls(100);
	REQUIRE(ls.primes.size() == 25);
}
