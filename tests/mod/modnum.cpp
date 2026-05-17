#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "mod/modnum.hpp"

namespace {
using i64 = long long;

constexpr int MOD_PRIME = 1'000'000'007;
constexpr int MOD_SMALL = 17;
constexpr int MOD_998   = 998'244'353;

using mint      = modnum<constant<int, MOD_PRIME>, naive_multiplier<int>>;
using msm       = modnum<constant<int, MOD_SMALL>, naive_multiplier<int>>;
using mint_b    = modnum<constant<int, MOD_PRIME>, barrett_multiplier<int>>;
using m998      = modnum<constant<int, MOD_998>,   naive_multiplier<int>>;
using m998_b    = modnum<constant<int, MOD_998>,   barrett_multiplier<int>>;
using mint_dyn  = modnum<inconstant<int>, naive_multiplier<int>>;
using mint64    = modnum<constant<int64_t, (int64_t(1) << 61) - 1>, naive_multiplier<int64_t>>;

static i64 ref_pow(i64 a, i64 b, i64 m) {
	a %= m; if (a < 0) a += m;
	i64 r = 1 % m;
	while (b > 0) {
		if (b & 1) r = r * a % m;
		a = a * a % m;
		b >>= 1;
	}
	return r;
}
} // namespace

TEST_CASE("modnum: default ctor is zero", "[modnum]") {
	mint x;
	REQUIRE(int(x) == 0);
	REQUIRE(is_zero(x));
	REQUIRE(mint::mod() == MOD_PRIME);
}

TEST_CASE("modnum: normalize negatives and large", "[modnum]") {
	REQUIRE(int(mint(-1)) == MOD_PRIME - 1);
	REQUIRE(int(mint(-MOD_PRIME)) == 0);
	REQUIRE(int(mint(MOD_PRIME)) == 0);
	REQUIRE(int(mint(i64(MOD_PRIME) * 5 + 3)) == 3);
	REQUIRE(int(mint(i64(- MOD_PRIME) * 5 - 3)) == MOD_PRIME - 3);
	REQUIRE(int(mint(0)) == 0);
	REQUIRE(int(mint(MOD_PRIME - 1)) == MOD_PRIME - 1);
	REQUIRE(int(mint(i64(- 1) * MOD_PRIME * 1000)) == 0);
}

TEST_CASE("modnum: basic arithmetic", "[modnum]") {
	mint a = 10, b = 7;
	REQUIRE(int(a + b) == 17);
	REQUIRE(int(a - b) == 3);
	REQUIRE(int(b - a) == MOD_PRIME - 3);
	REQUIRE(int(a * b) == 70);
	REQUIRE(int(- a) == MOD_PRIME - 10);
	REQUIRE(int(+ a) == 10);
	REQUIRE(int(- mint(0)) == 0);
	REQUIRE(int(mint(0) - mint(0)) == 0);
	REQUIRE(int(mint(MOD_PRIME - 1) + mint(1)) == 0);
	REQUIRE(int(mint(0) - mint(1)) == MOD_PRIME - 1);
}

TEST_CASE("modnum: comparison operators", "[modnum]") {
	mint a = 3, b = 5, c = 3;
	REQUIRE(a == c);
	REQUIRE(a != b);
	REQUIRE(a < b);
	REQUIRE(b > a);
	REQUIRE(a <= c);
	REQUIRE(a >= c);
	REQUIRE_FALSE(a < c);
	REQUIRE_FALSE(a > c);
}

TEST_CASE("modnum: increment / decrement", "[modnum]") {
	mint a = MOD_PRIME - 1;
	REQUIRE(int(++ a) == 0);
	REQUIRE(int(a ++) == 0);
	REQUIRE(int(a) == 1);
	REQUIRE(int(-- a) == 0);
	REQUIRE(int(a --) == 0);
	REQUIRE(int(a) == MOD_PRIME - 1);

	mint b = 0;
	REQUIRE(int(-- b) == MOD_PRIME - 1);
	REQUIRE(int(++ b) == 0);
}

TEST_CASE("modnum: division and inverse", "[modnum]") {
	std::mt19937 rng(42);
	for (int t = 0; t < 500; ++ t) {
		int av = int(rng() % unsigned(MOD_PRIME - 1)) + 1;
		int bv = int(rng() % unsigned(MOD_PRIME - 1)) + 1;
		mint a = av, b = bv;
		mint q = a / b;
		REQUIRE(int(q * b) == av);
	}
	for (int b = 1; b <= 200; ++ b) {
		mint inv = mint(1) / mint(b);
		REQUIRE(int(mint(b) * inv) == 1);
	}
}

TEST_CASE("modnum: bin_pow vs naive", "[modnum]") {
	std::mt19937_64 rng(123);
	for (int t = 0; t < 500; ++ t) {
		i64 a = i64(rng() % unsigned(MOD_PRIME));
		i64 b = i64(rng() % 1000u);
		REQUIRE(int(bin_pow(mint(a), b)) == int(ref_pow(a, b, MOD_PRIME)));
	}
	REQUIRE(int(bin_pow(mint(0), 0)) == 1);
	REQUIRE(int(bin_pow(mint(0), 5)) == 0);
	REQUIRE(int(bin_pow(mint(123), 0)) == 1);
	REQUIRE(int(bin_pow(mint(1), 1000000)) == 1);
	REQUIRE(int(bin_pow(mint(2), MOD_PRIME - 1)) == 1);
}

TEST_CASE("modnum: small modulus exhaustive", "[modnum]") {
	for (int a = 0; a < MOD_SMALL; ++ a) {
		for (int b = 0; b < MOD_SMALL; ++ b) {
			REQUIRE(int(msm(a) + msm(b)) == (a + b) % MOD_SMALL);
			REQUIRE(int(msm(a) * msm(b)) == (a * b) % MOD_SMALL);
			int diff = ((a - b) % MOD_SMALL + MOD_SMALL) % MOD_SMALL;
			REQUIRE(int(msm(a) - msm(b)) == diff);
		}
	}
	for (int a = 1; a < MOD_SMALL; ++ a) {
		msm inv = msm(1) / msm(a);
		REQUIRE(int(msm(a) * inv) == 1);
	}
}

TEST_CASE("modnum: barrett multiplier agrees with naive", "[modnum]") {
	std::mt19937 rng(7);
	for (int t = 0; t < 1000; ++ t) {
		int a = int(rng() % unsigned(MOD_PRIME));
		int b = int(rng() % unsigned(MOD_PRIME));
		REQUIRE(int(mint_b(a) * mint_b(b)) == int(mint(a) * mint(b)));
		REQUIRE(int(m998_b(a) * m998_b(b)) == int(m998(a) * m998(b)));
	}
}

TEST_CASE("modnum: barrett exhaustive small range", "[modnum]") {
	for (int a = 0; a < 100; ++ a) {
		for (int b = 0; b < 100; ++ b) {
			REQUIRE(int(mint_b(a) * mint_b(b)) == int(mint(a) * mint(b)));
		}
	}
}

TEST_CASE("modnum: dynamic modulus (inconstant)", "[modnum]") {
	inconstant<int>::set_value(13);
	REQUIRE(mint_dyn::mod() == 13);
	for (int a = 0; a < 13; ++ a) {
		for (int b = 0; b < 13; ++ b) {
			REQUIRE(int(mint_dyn(a) + mint_dyn(b)) == (a + b) % 13);
			REQUIRE(int(mint_dyn(a) * mint_dyn(b)) == (a * b) % 13);
		}
	}
	inconstant<int>::set_value(MOD_PRIME);
	REQUIRE(mint_dyn::mod() == MOD_PRIME);
	REQUIRE(int(mint_dyn(MOD_PRIME - 1) + mint_dyn(2)) == 1);
}

TEST_CASE("modnum: to_string / abs / call operator", "[modnum]") {
	mint a = 42;
	REQUIRE(to_string(a) == "42");
	REQUIRE(abs(a) == 42);
	REQUIRE(a() == 42);
	mint b = 0;
	REQUIRE(to_string(b) == "0");
	REQUIRE(is_zero(b));
}

TEST_CASE("modnum: compound assignments", "[modnum]") {
	mint a = 10;
	a += mint(20); REQUIRE(int(a) == 30);
	a -= mint(5);  REQUIRE(int(a) == 25);
	a *= mint(4);  REQUIRE(int(a) == 100);
	a /= mint(4);  REQUIRE(int(a) == 25);
	a += mint(- 100); REQUIRE(int(a) == (MOD_PRIME - 75) % MOD_PRIME);
}

TEST_CASE("modnum: 64-bit modulus arithmetic", "[modnum]") {
	// naive_multiplier for 64-bit uses long-double trick; only safe for
	// mod * mod fitting comfortably. Keep operands below 2^30 for safety.
	constexpr int64_t M = (int64_t(1) << 61) - 1;
	std::mt19937_64 rng(2024);
	for (int t = 0; t < 500; ++ t) {
		int64_t a = int64_t(rng() % (uint64_t(1) << 30));
		int64_t b = int64_t(rng() % (uint64_t(1) << 30));
		int64_t got = int64_t(mint64(a) * mint64(b));
		int64_t ref = int64_t((__int128(a) * b) % M);
		REQUIRE(got == ref);
	}
}

TEST_CASE("modnum: distributive and identity laws (stress)", "[modnum][stress]") {
	std::mt19937 rng(2025);
	for (int t = 0; t < 1000; ++ t) {
		int av = int(rng() % unsigned(MOD_PRIME));
		int bv = int(rng() % unsigned(MOD_PRIME));
		int cv = int(rng() % unsigned(MOD_PRIME));
		mint a = av, b = bv, c = cv;
		REQUIRE(a * (b + c) == a * b + a * c);
		REQUIRE((a + b) * c == a * c + b * c);
		REQUIRE(a + mint(0) == a);
		REQUIRE(a * mint(1) == a);
		REQUIRE(a * mint(0) == mint(0));
	}
}

TEST_CASE("modnum: stream output", "[modnum]") {
	std::ostringstream oss;
	std::ostream& os = oss;
	os << mint(123);
	REQUIRE(oss.str() == "123");
}

TEST_CASE("modnum: copy and assign", "[modnum]") {
	mint a = 42;
	mint b = a;
	REQUIRE(int(b) == 42);
	mint c;
	c = a;
	REQUIRE(int(c) == 42);
	a = mint(99);
	REQUIRE(int(b) == 42);
	REQUIRE(int(c) == 42);
}

TEST_CASE("modnum: division by inverse path", "[modnum]") {
	std::mt19937 rng(31415);
	for (int t = 0; t < 200; ++ t) {
		int b = int(rng() % unsigned(MOD_PRIME - 1)) + 1;
		mint inv = mint(1) / mint(b);
		i64 expected = ref_pow(b, MOD_PRIME - 2, MOD_PRIME);
		REQUIRE(int(inv) == int(expected));
	}
}

TEST_CASE("modnum: barrett + naive mixed via constants", "[modnum]") {
	for (int a = 0; a < 50; ++ a) {
		REQUIRE(int(mint(a)) == int(mint_b(a)));
		REQUIRE(int(- mint(a)) == int(- mint_b(a)));
		REQUIRE(int(mint(a) + mint(a)) == int(mint_b(a) + mint_b(a)));
	}
}
