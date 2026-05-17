#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "mod/modnum.hpp"
#include "math/miller_rabin.hpp"
#include "math/factorizer.hpp"

namespace {
static std::vector<std::pair<uint64_t, int>> trial_factor(uint64_t x) {
	std::vector<std::pair<uint64_t, int>> out;
	for (uint64_t p = 2; p * p <= x; ++ p) {
		if (x % p == 0) {
			int c = 0;
			while (x % p == 0) { x /= p; ++ c; }
			out.push_back({p, c});
		}
	}
	if (x > 1) out.push_back({x, 1});
	return out;
}

struct sieve_init {
	sieve_init() { factorizer::linear_sieve(10000); }
};
static sieve_init _init;
} // namespace

TEST_CASE("factorizer: is_prime cache", "[factorizer]") {
	REQUIRE(factorizer::is_prime<int>(2));
	REQUIRE(factorizer::is_prime<int>(7));
	REQUIRE(factorizer::is_prime<int>(997));
	REQUIRE(factorizer::is_prime<int>(7919));
	REQUIRE_FALSE(factorizer::is_prime<int>(1));
	REQUIRE_FALSE(factorizer::is_prime<int>(100));
	REQUIRE_FALSE(factorizer::is_prime<int>(9999));
}

TEST_CASE("factorizer: factorize small numbers", "[factorizer]") {
	for (int x : {2, 3, 12, 60, 720, 1024, 999, 9999, 8191}) {
		auto got = factorizer::factorize<uint64_t>(uint64_t(x));
		REQUIRE(got == trial_factor(uint64_t(x)));
	}
}

TEST_CASE("factorizer: factorize stress small", "[factorizer][stress]") {
	std::mt19937 rng(7);
	for (int t = 0; t < 300; ++ t) {
		uint64_t x = 1 + uint64_t(rng() % 9999u);
		REQUIRE(factorizer::factorize<uint64_t>(x) == trial_factor(x));
	}
}

TEST_CASE("factorizer: pollard_rho for large semiprimes", "[factorizer][pollard]") {
	std::vector<uint64_t> xs = {
		uint64_t(1000003ULL) * 1000033ULL,
		uint64_t(999983ULL) * 999979ULL,
		uint64_t(1) << 40,
		uint64_t(3) * 1000000007ULL,
	};
	for (auto x : xs) {
		auto got = factorizer::factorize<uint64_t>(x);
		uint64_t prod = 1;
		for (auto [p, e] : got) {
			REQUIRE(miller_rabin::is_prime<uint64_t>(p));
			for (int i = 0; i < e; ++ i) prod *= p;
		}
		REQUIRE(prod == x);
	}
}

TEST_CASE("factorizer: find_divisor enumerates all", "[factorizer]") {
	for (int n : {1, 6, 12, 60, 100, 360, 720}) {
		auto fact = factorizer::factorize<uint64_t>(uint64_t(n));
		auto divs = factorizer::find_divisor(fact, true);
		std::vector<uint64_t> ref;
		for (int d = 1; d <= n; ++ d) if (n % d == 0) ref.push_back(uint64_t(d));
		REQUIRE(divs == ref);
	}
}

TEST_CASE("factorizer: prime factorization unique", "[factorizer][stress]") {
	std::mt19937 rng(99);
	for (int t = 0; t < 100; ++ t) {
		uint64_t x = 2 + uint64_t(rng() % 9998u);
		auto fact = factorizer::factorize<uint64_t>(x);
		for (auto [p, e] : fact) REQUIRE(miller_rabin::is_prime<uint64_t>(p));
		for (size_t i = 1; i < fact.size(); ++ i) REQUIRE(fact[i - 1].first < fact[i].first);
	}
}
