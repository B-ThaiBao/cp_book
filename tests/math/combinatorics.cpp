#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "mod/modnum.hpp"
#include "math/combinatorics.hpp"

namespace {
using mint = modnum<constant<int, 1000000007>, naive_multiplier<int>>;
} // namespace

TEST_CASE("combinatorics: factorial values", "[combinatorics]") {
	combinatorics<mint> C(21);
	long long f = 1;
	for (int i = 0; i <= 20; ++ i) {
		REQUIRE(int(C.fact(i)()) == int(f % 1000000007));
		f = f * (i + 1) % 1000000007;
	}
}

TEST_CASE("combinatorics: inv_fact is inverse of fact", "[combinatorics]") {
	combinatorics<mint> C(50);
	for (int i = 0; i < 50; ++ i) REQUIRE(int((C.fact(i) * C.inv_fact(i))()) == 1);
}

TEST_CASE("combinatorics: choose basic", "[combinatorics]") {
	combinatorics<mint> C(30);
	REQUIRE(int(C.choose(5, 2)()) == 10);
	REQUIRE(int(C.choose(10, 3)()) == 120);
	REQUIRE(int(C.choose(10, 0)()) == 1);
	REQUIRE(int(C.choose(10, 10)()) == 1);
	REQUIRE(int(C.choose(20, 10)()) == 184756);
}

TEST_CASE("combinatorics: choose out of range", "[combinatorics]") {
	combinatorics<mint> C(30);
	REQUIRE(int(C.choose(5, - 1)()) == 0);
	REQUIRE(int(C.choose(5, 6)()) == 0);
	REQUIRE(int(C.choose(- 1, 0)()) == 0);
}

TEST_CASE("combinatorics: permute basic", "[combinatorics]") {
	combinatorics<mint> C(20);
	REQUIRE(int(C.permute(5, 2)()) == 20);
	REQUIRE(int(C.permute(5, 5)()) == 120);
	REQUIRE(int(C.permute(10, 3)()) == 720);
	REQUIRE(int(C.permute(7, 0)()) == 1);
}

TEST_CASE("combinatorics: catalan numbers", "[combinatorics][catalan]") {
	combinatorics<mint> C(40);
	std::vector<long long> cat = {1, 1, 2, 5, 14, 42, 132, 429, 1430, 4862, 16796, 58786, 208012};
	for (int i = 0; i < int(cat.size()); ++ i) REQUIRE(int(C.catalan(i)()) == int(cat[i] % 1000000007));
}

TEST_CASE("combinatorics: pascal triangle stress", "[combinatorics][stress]") {
	int N = 30;
	combinatorics<mint> C(N + 1);
	std::vector<std::vector<long long>> p(N + 1, std::vector<long long>(N + 1, 0));
	for (int i = 0; i <= N; ++ i) {
		p[i][0] = 1;
		for (int j = 1; j <= i; ++ j)
			p[i][j] = (p[i - 1][j - 1] + (j <= i - 1 ? p[i - 1][j] : 0)) % 1000000007;
	}
	for (int n = 0; n <= N; ++ n)
		for (int k = 0; k <= n; ++ k)
			REQUIRE(int(C.choose(n, k)()) == int(p[n][k]));
}

TEST_CASE("combinatorics: symmetric choose(n,k)==choose(n,n-k)", "[combinatorics][stress]") {
	combinatorics<mint> C(40);
	for (int n = 0; n < 40; ++ n)
		for (int k = 0; k <= n; ++ k)
			REQUIRE(int(C.choose(n, k)()) == int(C.choose(n, n - k)()));
}

TEST_CASE("combinatorics: row sums to 2^n", "[combinatorics][stress]") {
	combinatorics<mint> C(20);
	for (int n = 0; n < 20; ++ n) {
		mint s(0);
		for (int k = 0; k <= n; ++ k) s += C.choose(n, k);
		REQUIRE(int(s()) == (1 << n));
	}
}
