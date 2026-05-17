#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "fft/fft.hpp"
#include "mod/modnum.hpp"

namespace {
using ll = long long;
using mint = modnum<constant<int, 998244353>, naive_multiplier<int>>;

static std::vector<ll> naive_ll(const std::vector<ll>& a, const std::vector<ll>& b) {
	std::vector<ll> r(a.size() + b.size() - 1, 0);
	for (size_t i = 0; i < a.size(); ++ i)
		for (size_t j = 0; j < b.size(); ++ j)
			r[i + j] += a[i] * b[j];
	return r;
}

static std::vector<mint> naive_mint(const std::vector<mint>& a, const std::vector<mint>& b) {
	std::vector<mint> r(a.size() + b.size() - 1, mint(0));
	for (size_t i = 0; i < a.size(); ++ i)
		for (size_t j = 0; j < b.size(); ++ j)
			r[i + j] += a[i] * b[j];
	return r;
}

static std::vector<ll> rand_ll(std::mt19937& rng, int N, ll lo, ll hi) {
	std::vector<ll> a(N);
	ll span = hi - lo + 1;
	for (auto& x : a) x = lo + ll(rng() % uint32_t(span));
	return a;
}
} // namespace

TEST_CASE("fft::naive_multiply: small example", "[fft]") {
	std::vector<ll> a = {1, 2, 3}, b = {4, 5};
	auto r = fft::naive_multiply<ll>(a, b);
	REQUIRE(r == std::vector<ll>({4, 13, 22, 15}));
}

TEST_CASE("fft::naive_multiply: identity", "[fft]") {
	std::vector<ll> a = {3, 0, 5, -1};
	std::vector<ll> id = {1};
	auto r = fft::naive_multiply<ll>(a, id);
	REQUIRE(r == a);
}

TEST_CASE("fft::naive_multiply: random stress", "[fft][stress]") {
	std::mt19937 rng(31);
	for (int t = 0; t < 50; ++ t) {
		int N = 1 + int(rng() % 30u), M = 1 + int(rng() % 30u);
		auto a = rand_ll(rng, N, -100, 100);
		auto b = rand_ll(rng, M, -100, 100);
		REQUIRE(fft::naive_multiply<ll>(a, b) == naive_ll(a, b));
	}
}

TEST_CASE("fft::fft_complex_multiply: matches naive small ints", "[fft][stress]") {
	std::mt19937 rng(7);
	for (int t = 0; t < 30; ++ t) {
		int N = 1 + int(rng() % 30u), M = 1 + int(rng() % 30u);
		auto a = rand_ll(rng, N, 0, 100);
		auto b = rand_ll(rng, M, 0, 100);
		auto ref = naive_ll(a, b);
		auto got = fft::fft_complex_multiply<ll>(a, b);
		REQUIRE(got == ref);
	}
}

TEST_CASE("fft::fft_complex_multiply: signed values", "[fft][stress]") {
	std::mt19937 rng(99);
	for (int t = 0; t < 30; ++ t) {
		int N = 1 + int(rng() % 30u), M = 1 + int(rng() % 30u);
		auto a = rand_ll(rng, N, -50, 50);
		auto b = rand_ll(rng, M, -50, 50);
		REQUIRE(fft::fft_complex_multiply<ll>(a, b) == naive_ll(a, b));
	}
}

TEST_CASE("fft::fft_complex_double_multiply: matches naive larger values", "[fft][stress]") {
	std::mt19937 rng(11);
	for (int t = 0; t < 20; ++ t) {
		int N = 1 + int(rng() % 50u), M = 1 + int(rng() % 50u);
		auto a = rand_ll(rng, N, 0, 1000);
		auto b = rand_ll(rng, M, 0, 1000);
		auto ref = naive_ll(a, b);
		auto got = fft::fft_complex_double_multiply<ll>(a, b);
		REQUIRE(got.size() == ref.size());
		for (size_t i = 0; i < ref.size(); ++ i) REQUIRE(std::abs(got[i] - ref[i]) <= 1);
	}
}

TEST_CASE("fft::fft_numeric_multiply: matches naive (NTT)", "[fft][stress]") {
	std::mt19937 rng(2024);
	for (int t = 0; t < 30; ++ t) {
		int N = 1 + int(rng() % 30u), M = 1 + int(rng() % 30u);
		std::vector<mint> a(N), b(M);
		for (auto& x : a) x = mint(int(rng() % 1000u));
		for (auto& x : b) x = mint(int(rng() % 1000u));
		auto ref = naive_mint(a, b);
		auto got = fft::fft_numeric_multiply<mint>(a, b);
		REQUIRE(got.size() == ref.size());
		for (size_t i = 0; i < got.size(); ++ i) REQUIRE(got[i] == ref[i]);
	}
}

TEST_CASE("fft::fft_numeric_multiply: identity", "[fft]") {
	std::vector<mint> a = {mint(3), mint(0), mint(5), mint(7)};
	std::vector<mint> id = {mint(1)};
	auto r = fft::fft_numeric_multiply<mint>(a, id);
	REQUIRE(r.size() == a.size());
	for (size_t i = 0; i < a.size(); ++ i) REQUIRE(r[i] == a[i]);
}

TEST_CASE("fft::fft_numeric_multiply: associativity", "[fft][stress]") {
	std::mt19937 rng(57);
	for (int t = 0; t < 10; ++ t) {
		int N = 1 + int(rng() % 10u);
		int M = 1 + int(rng() % 10u);
		int K = 1 + int(rng() % 10u);
		std::vector<mint> a(N), b(M), c(K);
		for (auto& x : a) x = mint(int(rng() % 100u));
		for (auto& x : b) x = mint(int(rng() % 100u));
		for (auto& x : c) x = mint(int(rng() % 100u));
		auto ab  = fft::fft_numeric_multiply<mint>(a, b);
		auto bc  = fft::fft_numeric_multiply<mint>(b, c);
		auto abc1 = fft::fft_numeric_multiply<mint>(ab, c);
		auto abc2 = fft::fft_numeric_multiply<mint>(a, bc);
		REQUIRE(abc1.size() == abc2.size());
		for (size_t i = 0; i < abc1.size(); ++ i) REQUIRE(abc1[i] == abc2[i]);
	}
}

TEST_CASE("fft::fft_complex_mod_multiply: matches naive", "[fft][stress]") {
	std::mt19937 rng(123);
	for (int t = 0; t < 15; ++ t) {
		int N = 1 + int(rng() % 30u), M = 1 + int(rng() % 30u);
		std::vector<mint> a(N), b(M);
		for (auto& x : a) x = mint(int(rng() % 1000000u));
		for (auto& x : b) x = mint(int(rng() % 1000000u));
		auto ref = naive_mint(a, b);
		auto got = fft::fft_complex_mod_multiply<mint>(a, b);
		REQUIRE(got.size() == ref.size());
		for (size_t i = 0; i < got.size(); ++ i) REQUIRE(got[i] == ref[i]);
	}
}

TEST_CASE("fft::fft_numeric_multiply: larger sizes", "[fft][stress]") {
	std::mt19937 rng(2025);
	for (int sz : {64, 128, 256}) {
		std::vector<mint> a(sz), b(sz);
		for (auto& x : a) x = mint(int(rng() % 1000u));
		for (auto& x : b) x = mint(int(rng() % 1000u));
		auto ref = naive_mint(a, b);
		auto got = fft::fft_numeric_multiply<mint>(a, b);
		REQUIRE(got.size() == ref.size());
		for (size_t i = 0; i < got.size(); ++ i) REQUIRE(got[i] == ref[i]);
	}
}

TEST_CASE("fft::fft_complex_multiply: one-sided lengths", "[fft]") {
	std::vector<ll> a = {1, 2, 3, 4, 5};
	std::vector<ll> b = {10};
	auto got = fft::fft_complex_multiply<ll>(a, b);
	REQUIRE(got == std::vector<ll>({10, 20, 30, 40, 50}));
}

TEST_CASE("fft::naive_multiply: empty-ish (size 1 elements)", "[fft]") {
	std::vector<ll> a = {0}, b = {0};
	auto r = fft::naive_multiply<ll>(a, b);
	REQUIRE(r == std::vector<ll>({0}));
}
