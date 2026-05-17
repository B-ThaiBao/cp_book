// Benchmarks for math/* (miller_rabin, factorizer, combinatorics, sieve).
#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>

#include "mod/modnum.hpp"
#include "math/miller_rabin.hpp"
#include "math/factorizer.hpp"
#include "math/combinatorics.hpp"
#include "math/sieve.hpp"

using mint = modnum<constant<int, 1000000007>, naive_multiplier<int>>;

TEST_CASE("miller_rabin: 5e4 primality tests on ~62-bit numbers", "[!benchmark][miller_rabin]") {
	constexpr int N = 50'000;
	std::mt19937_64 rng(101);
	std::vector<uint64_t> nums(N);
	for (auto& x : nums) x = (rng() | 1ULL) & ((1ULL << 62) - 1);

	BENCHMARK("miller_rabin::is_prime x 5e4") {
		int acc = 0;
		for (auto& x : nums) acc += miller_rabin::is_prime<uint64_t>(x);
		return acc;
	};
}

TEST_CASE("factorizer: 1e4 Pollard-rho factorizations (~60-bit)", "[!benchmark][factorizer]") {
	constexpr int N = 10'000;
	std::mt19937_64 rng(103);
	std::vector<uint64_t> nums(N);
	for (auto& x : nums) x = (rng() | 1ULL) & ((1ULL << 60) - 1);

	BENCHMARK("factorize x 1e4") {
		long long acc = 0;
		for (auto& x : nums) {
			auto f = factorizer::factorize<uint64_t>(x);
			acc += (long long)f.size();
		}
		return acc;
	};
}

TEST_CASE("sieve: build prime sieves up to 1e7", "[!benchmark][sieve]") {
	BENCHMARK("eratos_sieve(1e7)") {
		auto p = eratos_sieve(10'000'000);
		return p.size();
	};
	BENCHMARK("block_sieve(1e7)") {
		auto primes = block_sieve(10'000'000);
		return primes.size();
	};
	BENCHMARK("linear_sieve(1e7)") {
		linear_sieve ls(10'000'000);
		return ls.primes.size();
	};
}

TEST_CASE("combinatorics: build (N=1e6) + 1e6 binomial queries", "[!benchmark][combinatorics]") {
	constexpr int N = 1'000'000;
	std::mt19937_64 rng(107);
	std::vector<std::pair<int, int>> q(N);
	for (auto& [n, k] : q) {
		n = 1 + int(rng() % unsigned(N));
		k = int(rng() % unsigned(n + 1));
	}

	BENCHMARK("combinatorics<mint>(N+1) + 1e6 choose(n,k)") {
		combinatorics<mint> C(N + 1);
		uint64_t acc = 0;
		for (auto& [n, k] : q) acc += (uint64_t)int(C.choose(n, k)());
		return acc;
	};
}
