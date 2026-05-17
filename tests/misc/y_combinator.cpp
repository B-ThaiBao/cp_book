#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "misc/y_combinator.hpp"

TEST_CASE("y_combinator: factorial", "[y_combinator]") {
	auto fact = std::y_combinator([](auto self, int n) -> long long {
		return n <= 1 ? 1LL : n * self(n - 1);
	});
	REQUIRE(fact(0)  == 1);
	REQUIRE(fact(1)  == 1);
	REQUIRE(fact(2)  == 2);
	REQUIRE(fact(5)  == 120);
	REQUIRE(fact(10) == 3628800LL);
	REQUIRE(fact(15) == 1307674368000LL);
	REQUIRE(fact(20) == 2432902008176640000LL);
}

TEST_CASE("y_combinator: fibonacci with memo", "[y_combinator]") {
	std::vector<long long> memo(60, -1);
	auto fib = std::y_combinator([&](auto self, int n) -> long long {
		if (n < 2) return n;
		if (memo[n] != -1) return memo[n];
		return memo[n] = self(n - 1) + self(n - 2);
	});
	REQUIRE(fib(0)  == 0);
	REQUIRE(fib(1)  == 1);
	REQUIRE(fib(2)  == 1);
	REQUIRE(fib(10) == 55);
	REQUIRE(fib(20) == 6765);
	REQUIRE(fib(30) == 832040);
	REQUIRE(fib(50) == 12586269025LL);
}

TEST_CASE("y_combinator: tree depth via dfs", "[y_combinator]") {
	std::vector<std::vector<int>> g(7);
	auto add = [&](int u, int v) { g[u].push_back(v); g[v].push_back(u); };
	add(0, 1); add(0, 2); add(1, 3); add(1, 4); add(2, 5); add(5, 6);

	int max_depth = 0;
	std::y_combinator([&](auto self, int u, int p, int d) -> void {
		max_depth = std::max(max_depth, d);
		for (int v : g[u]) if (v != p) self(v, u, d + 1);
	})(0, -1, 0);
	REQUIRE(max_depth == 3);
}

TEST_CASE("y_combinator: gcd", "[y_combinator]") {
	auto gcd = std::y_combinator([](auto self, long long a, long long b) -> long long {
		return b == 0 ? a : self(b, a % b);
	});
	REQUIRE(gcd(48, 18) == 6);
	REQUIRE(gcd(17, 5)  == 1);
	REQUIRE(gcd(100, 0) == 100);
	REQUIRE(gcd(0, 7)   == 7);
	REQUIRE(gcd(1000000, 1000) == 1000);
}

TEST_CASE("y_combinator: subtree sum via dfs", "[y_combinator]") {
	std::vector<std::vector<int>> g(6);
	auto add = [&](int u, int v) { g[u].push_back(v); g[v].push_back(u); };
	add(0, 1); add(0, 2); add(1, 3); add(1, 4); add(2, 5);
	std::vector<int> val = {1, 2, 3, 4, 5, 6};
	std::vector<int> sub(6, 0);
	std::y_combinator([&](auto self, int u, int p) -> int {
		sub[u] = val[u];
		for (int v : g[u]) if (v != p) sub[u] += self(v, u);
		return sub[u];
	})(0, -1);
	REQUIRE(sub[0] == 21);
	REQUIRE(sub[1] == 11);
	REQUIRE(sub[2] == 9);
	REQUIRE(sub[3] == 4);
	REQUIRE(sub[4] == 5);
	REQUIRE(sub[5] == 6);
}

TEST_CASE("y_combinator: void return permutations", "[y_combinator]") {
	std::vector<int> a = {1, 2, 3};
	int cnt = 0;
	std::y_combinator([&](auto self, int i) -> void {
		if (i == int(a.size())) { ++ cnt; return; }
		for (int j = i; j < int(a.size()); ++ j) {
			std::swap(a[i], a[j]);
			self(i + 1);
			std::swap(a[i], a[j]);
		}
	})(0);
	REQUIRE(cnt == 6);
}

TEST_CASE("y_combinator: mutable capture", "[y_combinator]") {
	int sum = 0;
	auto rec = std::y_combinator([&](auto self, int n) -> void {
		if (n <= 0) return;
		sum += n;
		self(n - 1);
	});
	rec(10);
	REQUIRE(sum == 55);
	rec(5);
	REQUIRE(sum == 55 + 15);
}

TEST_CASE("y_combinator: deeper recursion", "[y_combinator]") {
	auto rec = std::y_combinator([](auto self, int n) -> int {
		return n == 0 ? 0 : 1 + self(n - 1);
	});
	REQUIRE(rec(1000) == 1000);
	REQUIRE(rec(0) == 0);
}
