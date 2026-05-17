#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "mod/modnum.hpp"
#include "math/miller_rabin.hpp"

namespace {
static std::vector<bool> trial_sieve(const int& N) {
	std::vector<bool> p(N, true);
	if (N > 0) p[0] = false;
	if (N > 1) p[1] = false;
	for (int i = 2; i < N; ++ i)
		if (p[i]) for (long long j = (long long)i * i; j < N; j += i) p[j] = false;
	return p;
}
} // namespace

TEST_CASE("miller_rabin: small range vs sieve", "[miller_rabin]") {
	int N = 5000;
	auto ref = trial_sieve(N);
	for (int x = 0; x < N; ++ x)
		REQUIRE(miller_rabin::is_prime<uint64_t>(uint64_t(x)) == ref[x]);
}

TEST_CASE("miller_rabin: known large primes", "[miller_rabin]") {
	REQUIRE(miller_rabin::is_prime<uint64_t>(1000003ULL));
	REQUIRE(miller_rabin::is_prime<uint64_t>(1000000007ULL));
	REQUIRE(miller_rabin::is_prime<uint64_t>(1000000009ULL));
	REQUIRE(miller_rabin::is_prime<uint64_t>(998244353ULL));
	REQUIRE(miller_rabin::is_prime<uint64_t>(9223372036854775783ULL));
}

TEST_CASE("miller_rabin: Carmichael numbers", "[miller_rabin]") {
	REQUIRE_FALSE(miller_rabin::is_prime<uint64_t>(561ULL));
	REQUIRE_FALSE(miller_rabin::is_prime<uint64_t>(1105ULL));
	REQUIRE_FALSE(miller_rabin::is_prime<uint64_t>(1729ULL));
	REQUIRE_FALSE(miller_rabin::is_prime<uint64_t>(2821ULL));
	REQUIRE_FALSE(miller_rabin::is_prime<uint64_t>(6601ULL));
}

TEST_CASE("miller_rabin: large semiprimes are composite", "[miller_rabin]") {
	REQUIRE_FALSE(miller_rabin::is_prime<uint64_t>(1000003ULL * 1000033ULL));
	REQUIRE_FALSE(miller_rabin::is_prime<uint64_t>(999983ULL * 999979ULL));
}

TEST_CASE("miller_rabin: 0 and 1 are not prime", "[miller_rabin]") {
	REQUIRE_FALSE(miller_rabin::is_prime<uint64_t>(0));
	REQUIRE_FALSE(miller_rabin::is_prime<uint64_t>(1));
}

TEST_CASE("miller_rabin: 2 and 3 are prime", "[miller_rabin]") {
	REQUIRE(miller_rabin::is_prime<uint64_t>(2));
	REQUIRE(miller_rabin::is_prime<uint64_t>(3));
}
