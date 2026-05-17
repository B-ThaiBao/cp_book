#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "misc/ska_sort.hpp"

template <typename T>
static void check_sorted_equal(std::vector<T> a) {
	auto b = a;
	std::sort(b.begin(), b.end());
	std::ska_sort(a.begin(), a.end());
	REQUIRE(a == b);
}

TEST_CASE("ska_sort: empty / single", "[ska_sort]") {
	std::vector<int> a;
	std::ska_sort(a.begin(), a.end());
	REQUIRE(a.empty());
	std::vector<int> b = {42};
	std::ska_sort(b.begin(), b.end());
	REQUIRE(b == std::vector<int>{42});
}

TEST_CASE("ska_sort: two elements", "[ska_sort]") {
	check_sorted_equal<int>({2, 1});
	check_sorted_equal<int>({1, 2});
	check_sorted_equal<int>({-1, 1});
	check_sorted_equal<int>({0, 0});
}

TEST_CASE("ska_sort: unsigned int", "[ska_sort]") {
	std::vector<unsigned int> a = {5u, 3u, 9u, 0u, 7u, 1u, 5u};
	std::vector<unsigned int> b = a;
	std::sort(b.begin(), b.end());
	std::ska_sort(a.begin(), a.end());
	REQUIRE(a == b);
}

TEST_CASE("ska_sort: signed int with negatives", "[ska_sort]") {
	check_sorted_equal<int>({-3, 7, -1, 0, 5, -10, 4});
	check_sorted_equal<int>({INT_MIN, INT_MAX, 0});
	check_sorted_equal<int>({-1, -1, -1, 0, 0, 1, 1});
}

TEST_CASE("ska_sort: long long", "[ska_sort]") {
	check_sorted_equal<long long>({(1LL << 40), -(1LL << 35), 0, 12345, -1});
	check_sorted_equal<long long>({LLONG_MIN, 0, LLONG_MAX});
}

TEST_CASE("ska_sort: doubles", "[ska_sort]") {
	std::vector<double> a = {3.14, -1.5, 0.0, 2.71, -0.0, 1e9, -1e9};
	std::vector<double> b = a;
	std::sort(b.begin(), b.end());
	std::ska_sort(a.begin(), a.end());
	REQUIRE(a.size() == b.size());
	for (size_t i = 0; i < a.size(); ++ i) REQUIRE(a[i] == b[i]);
}

TEST_CASE("ska_sort: random large stress vs std::sort (int)", "[ska_sort][stress]") {
	std::mt19937 rng(98765);
	for (int trial = 0; trial < 10; ++ trial) {
		int n = 2000 + int(rng() % 3000u);
		std::vector<int> a(n);
		for (auto& x : a) x = int(rng()) - int(rng());
		auto b = a;
		std::sort(b.begin(), b.end());
		std::ska_sort(a.begin(), a.end());
		REQUIRE(a == b);
	}
}

TEST_CASE("ska_sort: random small stress (signed)", "[ska_sort][stress]") {
	std::mt19937 rng(31);
	for (int trial = 0; trial < 100; ++ trial) {
		int n = int(rng() % 200u);
		std::vector<int> a(n);
		for (auto& x : a) x = int(rng() % 200u) - 100;
		auto b = a;
		std::sort(b.begin(), b.end());
		std::ska_sort(a.begin(), a.end());
		REQUIRE(a == b);
	}
}

TEST_CASE("ska_sort: random stress unsigned long long", "[ska_sort][stress]") {
	std::mt19937_64 rng(7);
	for (int trial = 0; trial < 20; ++ trial) {
		int n = 500 + int(rng() % 1500u);
		std::vector<unsigned long long> a(n);
		for (auto& x : a) x = rng();
		auto b = a;
		std::sort(b.begin(), b.end());
		std::ska_sort(a.begin(), a.end());
		REQUIRE(a == b);
	}
}

TEST_CASE("ska_sort: with custom extract_key", "[ska_sort]") {
	std::vector<std::pair<int, int>> a = {{3, 1}, {1, 2}, {2, 3}, {1, 4}};
	std::ska_sort(a.begin(), a.end(), [](const std::pair<int, int>& p) { return p.first; });
	std::vector<int> firsts;
	for (auto& p : a) firsts.push_back(p.first);
	REQUIRE(firsts == std::vector<int>({1, 1, 2, 3}));
}

TEST_CASE("ska_sort: custom extract_key on long long", "[ska_sort]") {
	std::vector<std::pair<long long, std::string>> a = {{100, "a"}, {-5, "b"}, {0, "c"}, {-5, "d"}};
	std::ska_sort(a.begin(), a.end(), [](const auto& p) { return p.first; });
	std::vector<long long> keys;
	for (auto& p : a) keys.push_back(p.first);
	REQUIRE(keys == std::vector<long long>({-5, -5, 0, 100}));
}

TEST_CASE("ska_sort_copy: unsigned ints", "[ska_sort_copy]") {
	std::mt19937 rng(11);
	int n = 4000;
	std::vector<unsigned int> a(n), buf(n);
	for (auto& x : a) x = unsigned(rng());
	auto expect = a;
	std::sort(expect.begin(), expect.end());

	bool result_in_buf = std::ska_sort_copy(a.begin(), a.end(), buf.begin());
	const auto& got = result_in_buf ? buf : a;
	REQUIRE(std::vector<unsigned int>(got.begin(), got.begin() + n) == expect);
}

TEST_CASE("ska_sort: strings", "[ska_sort]") {
	std::vector<std::string> a = {"pear", "apple", "banana", "apple", "cherry", ""};
	auto b = a;
	std::sort(b.begin(), b.end());
	std::ska_sort(a.begin(), a.end());
	REQUIRE(a == b);
}

TEST_CASE("ska_sort: strings random stress", "[ska_sort][stress]") {
	std::mt19937 rng(123);
	for (int trial = 0; trial < 30; ++ trial) {
		int n = int(rng() % 100u) + 1;
		std::vector<std::string> a(n);
		for (auto& s : a) {
			int len = int(rng() % 10u);
			s.resize(len);
			for (auto& c : s) c = char('a' + rng() % 5u);
		}
		auto b = a;
		std::sort(b.begin(), b.end());
		std::ska_sort(a.begin(), a.end());
		REQUIRE(a == b);
	}
}

TEST_CASE("ska_sort: already sorted", "[ska_sort]") {
	std::vector<int> a;
	for (int i = 0; i < 1000; ++ i) a.push_back(i);
	auto b = a;
	std::ska_sort(a.begin(), a.end());
	REQUIRE(a == b);
}

TEST_CASE("ska_sort: reverse sorted", "[ska_sort]") {
	std::vector<int> a;
	for (int i = 999; i >= 0; -- i) a.push_back(i);
	auto b = a;
	std::sort(b.begin(), b.end());
	std::ska_sort(a.begin(), a.end());
	REQUIRE(a == b);
}

TEST_CASE("ska_sort: all equal", "[ska_sort]") {
	std::vector<int> a(500, 7);
	std::ska_sort(a.begin(), a.end());
	REQUIRE(a == std::vector<int>(500, 7));
}
