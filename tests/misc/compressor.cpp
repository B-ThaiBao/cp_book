#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "misc/compressor.hpp"

TEST_CASE("compressor: sort + compress + lookup", "[compressor]") {
	compressor<int> c{5, 1, 3, 1, 5, 2};
	std::sort(c.begin(), c.end());
	c.compress();
	REQUIRE(c.size() == 4);
	REQUIRE(c(1) == 0);
	REQUIRE(c(2) == 1);
	REQUIRE(c(3) == 2);
	REQUIRE(c(5) == 3);
	REQUIRE(c(0) == 0);
	REQUIRE(c(4) == 3);
	REQUIRE(c(100) == 4);
}

TEST_CASE("compressor: empty", "[compressor]") {
	compressor<int> c;
	c.compress();
	REQUIRE(c.empty());
	REQUIRE(c(42) == 0);
	REQUIRE(c(0)  == 0);
	REQUIRE(c(-100) == 0);
}

TEST_CASE("compressor: single value", "[compressor]") {
	compressor<int> c{7};
	c.compress();
	REQUIRE(c.size() == 1);
	REQUIRE(c(6) == 0);
	REQUIRE(c(7) == 0);
	REQUIRE(c(8) == 1);
}

TEST_CASE("compressor: strings", "[compressor]") {
	compressor<std::string> c{"banana", "apple", "cherry", "apple"};
	std::sort(c.begin(), c.end());
	c.compress();
	REQUIRE(c.size() == 3);
	REQUIRE(c(std::string("apple"))  == 0);
	REQUIRE(c(std::string("banana")) == 1);
	REQUIRE(c(std::string("cherry")) == 2);
	REQUIRE(c(std::string("date"))   == 3);
}

TEST_CASE("compressor: long long with negatives", "[compressor]") {
	compressor<long long> c{-5, 100, 0, -5, 100, 7};
	std::sort(c.begin(), c.end());
	c.compress();
	REQUIRE(c.size() == 4);
	REQUIRE(c(-5)  == 0);
	REQUIRE(c(0)   == 1);
	REQUIRE(c(7)   == 2);
	REQUIRE(c(100) == 3);
}

TEST_CASE("compressor: random stress", "[compressor][stress]") {
	std::mt19937 rng(1);
	for (int t = 0; t < 50; ++ t) {
		int N = 1 + int(rng() % 100u);
		compressor<int> c;
		std::vector<int> src;
		for (int i = 0; i < N; ++ i) {
			int v = int(rng() % 30u) - 15;
			c.push_back(v);
			src.push_back(v);
		}
		std::sort(c.begin(), c.end());
		c.compress();

		std::sort(src.begin(), src.end());
		src.erase(std::unique(src.begin(), src.end()), src.end());
		REQUIRE(c.size() == src.size());
		for (int q = 0; q < 50; ++ q) {
			int x = int(rng() % 40u) - 20;
			int expected = int(std::lower_bound(src.begin(), src.end(), x) - src.begin());
			REQUIRE(c(x) == expected);
		}
	}
}

TEST_CASE("compress(): preserves order, assigns ranks", "[compress]") {
	std::vector<int> a = {10, 5, 5, 100, 7};
	auto r = compress(a);
	std::vector<int> expect = {2, 0, 0, 3, 1};
	REQUIRE(r == expect);
}

TEST_CASE("compress(): empty input", "[compress]") {
	std::vector<int> a;
	auto r = compress(a);
	REQUIRE(r.empty());
}

TEST_CASE("compress(): single element", "[compress]") {
	std::vector<int> a = {42};
	auto r = compress(a);
	REQUIRE(r == std::vector<int>({0}));
	auto r2 = compress(a, 5);
	REQUIRE(r2 == std::vector<int>({5}));
}

TEST_CASE("compress(): all same", "[compress]") {
	std::vector<int> a(10, 7);
	auto r = compress(a);
	REQUIRE(r == std::vector<int>(10, 0));
}

TEST_CASE("compress(): already sorted ascending", "[compress]") {
	std::vector<int> a = {1, 2, 3, 4, 5};
	auto r = compress(a);
	REQUIRE(r == std::vector<int>({0, 1, 2, 3, 4}));
}

TEST_CASE("compress(): already sorted descending", "[compress]") {
	std::vector<int> a = {5, 4, 3, 2, 1};
	auto r = compress(a);
	REQUIRE(r == std::vector<int>({4, 3, 2, 1, 0}));
}

TEST_CASE("compress(): nxt offset", "[compress]") {
	std::vector<int> a = {1, 2, 1};
	auto r = compress(a, 10);
	REQUIRE(r == std::vector<int>({10, 11, 10}));
}

TEST_CASE("compress(): random stress vs naive", "[compress][stress]") {
	std::mt19937 rng(2024);
	for (int trial = 0; trial < 100; ++ trial) {
		int n = 1 + int(rng() % 100u);
		std::vector<int> a(n);
		for (auto& x : a) x = int(rng() % 40u) - 20;
		auto got = compress(a);

		auto sorted_unique = a;
		std::sort(sorted_unique.begin(), sorted_unique.end());
		sorted_unique.erase(std::unique(sorted_unique.begin(), sorted_unique.end()), sorted_unique.end());
		std::vector<int> expect(n);
		for (int i = 0; i < n; ++ i) {
			expect[i] = int(std::lower_bound(sorted_unique.begin(), sorted_unique.end(), a[i]) - sorted_unique.begin());
		}
		REQUIRE(got == expect);
	}
}

TEST_CASE("compress(): strings", "[compress]") {
	std::vector<std::string> a = {"cat", "apple", "banana", "apple", "cat"};
	auto r = compress(a);
	REQUIRE(r == std::vector<int>({2, 0, 1, 0, 2}));
}
