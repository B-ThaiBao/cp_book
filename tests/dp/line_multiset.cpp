#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "dp/line_multiset.hpp"

namespace {
using ll = long long;
using lines_t = std::vector<std::pair<ll, ll>>;

static ll brute_max(const lines_t& lines, ll x) {
	ll best = LLONG_MIN;
	for (auto& [a, b] : lines) best = std::max(best, a * x + b);
	return best;
}
static ll brute_min(const lines_t& lines, ll x) {
	ll best = LLONG_MAX;
	for (auto& [a, b] : lines) best = std::min(best, a * x + b);
	return best;
}
static ll query_max(line_multiset<line<ll>>& h, ll x) {
	auto it = h.lower_bound(x);
	return (*it)[0] * x + (*it)[1];
}
} // namespace

TEST_CASE("line_multiset: empty container", "[cht]") {
	line_multiset<line<ll>> h;
	REQUIRE(h.empty());
	REQUIRE(h.size() == 0);
}

TEST_CASE("line_multiset: single line", "[cht]") {
	line_multiset<line<ll>> h;
	h.add_line(3, -5);
	REQUIRE(h.size() == 1);
	for (ll x = -10; x <= 10; ++ x) REQUIRE(query_max(h, x) == 3 * x - 5);
}

TEST_CASE("line_multiset: two intersecting lines", "[cht]") {
	line_multiset<line<ll>> h;
	h.add_line(1, 0);
	h.add_line(-1, 10);
	lines_t lines = {{1, 0}, {-1, 10}};
	for (ll x = -20; x <= 20; ++ x) REQUIRE(query_max(h, x) == brute_max(lines, x));
}

TEST_CASE("line_multiset: documented example (max)", "[cht]") {
	line_multiset<line<ll>> h;
	h.add_line(2, 0);
	h.add_line(-1, 5);
	h.add_line(0, 3);
	lines_t lines = {{2, 0}, {-1, 5}, {0, 3}};
	for (ll x = -10; x <= 10; ++ x) REQUIRE(query_max(h, x) == brute_max(lines, x));
}

TEST_CASE("line_multiset: parallel lines keep maximum", "[cht]") {
	line_multiset<line<ll>> h;
	h.add_line(2, 1);
	h.add_line(2, 5);
	h.add_line(2, -10);
	lines_t lines = {{2, 1}, {2, 5}, {2, -10}};
	for (ll x = -50; x <= 50; ++ x) REQUIRE(query_max(h, x) == brute_max(lines, x));
}

TEST_CASE("line_multiset: horizontal lines only", "[cht]") {
	line_multiset<line<ll>> h;
	h.add_line(0, 3);
	h.add_line(0, 7);
	h.add_line(0, -2);
	for (ll x = -100; x <= 100; ++ x) REQUIRE(query_max(h, x) == 7);
}

TEST_CASE("line_multiset: dominated line ignored", "[cht]") {
	line_multiset<line<ll>> h;
	h.add_line(2, 0);
	h.add_line(2, 5);
	lines_t lines = {{2, 0}, {2, 5}};
	for (ll x = -10; x <= 10; ++ x) REQUIRE(query_max(h, x) == brute_max(lines, x));
}

TEST_CASE("line_multiset: min query via negation", "[cht]") {
	line_multiset<line<ll>> h;
	lines_t lines = {{2, 0}, {-1, 5}, {0, 3}, {3, -2}, {-2, 4}};
	for (auto& [a, b] : lines) h.add_line(- a, - b);
	for (ll x = -10; x <= 10; ++ x) {
		auto it = h.lower_bound(x);
		ll got = -((*it)[0] * x + (*it)[1]);
		REQUIRE(got == brute_min(lines, x));
	}
}

TEST_CASE("line_multiset: insertion order does not matter", "[cht]") {
	lines_t lines = {{1, 0}, {-2, 4}, {3, -5}, {0, 2}, {-1, 7}, {2, -3}};
	line_multiset<line<ll>> h1, h2, h3;
	for (auto& [a, b] : lines) h1.add_line(a, b);
	for (auto it = lines.rbegin(); it != lines.rend(); ++ it) h2.add_line(it->first, it->second);
	auto shuffled = lines;
	std::mt19937 rng(31);
	std::shuffle(shuffled.begin(), shuffled.end(), rng);
	for (auto& [a, b] : shuffled) h3.add_line(a, b);

	for (ll x = -25; x <= 25; ++ x) {
		ll ref = brute_max(lines, x);
		REQUIRE(query_max(h1, x) == ref);
		REQUIRE(query_max(h2, x) == ref);
		REQUIRE(query_max(h3, x) == ref);
	}
}

TEST_CASE("line_multiset: small random stress", "[cht][stress]") {
	std::mt19937 rng(57);
	for (int t = 0; t < 100; ++ t) {
		line_multiset<line<ll>> h;
		lines_t lines;
		int N = 1 + int(rng() % 15u);
		for (int i = 0; i < N; ++ i) {
			ll a = ll(rng() % 21u) - 10;
			ll b = ll(rng() % 201u) - 100;
			h.add_line(a, b);
			lines.push_back({a, b});
		}
		for (int q = 0; q < 50; ++ q) {
			ll x = ll(rng() % 41u) - 20;
			REQUIRE(query_max(h, x) == brute_max(lines, x));
		}
	}
}

TEST_CASE("line_multiset: large random stress", "[cht][stress]") {
	std::mt19937 rng(98765);
	for (int t = 0; t < 20; ++ t) {
		line_multiset<line<ll>> h;
		lines_t lines;
		int N = 100 + int(rng() % 100u);
		for (int i = 0; i < N; ++ i) {
			ll a = ll(rng() % 201u) - 100;
			ll b = ll(rng() % 2001u) - 1000;
			h.add_line(a, b);
			lines.push_back({a, b});
		}
		for (int q = 0; q < 200; ++ q) {
			ll x = ll(rng() % 401u) - 200;
			REQUIRE(query_max(h, x) == brute_max(lines, x));
		}
	}
}

TEST_CASE("line_multiset: random min stress via negation", "[cht][stress]") {
	std::mt19937 rng(2024);
	for (int t = 0; t < 30; ++ t) {
		line_multiset<line<ll>> h;
		lines_t lines;
		int N = 1 + int(rng() % 30u);
		for (int i = 0; i < N; ++ i) {
			ll a = ll(rng() % 41u) - 20;
			ll b = ll(rng() % 401u) - 200;
			lines.push_back({a, b});
			h.add_line(- a, - b);
		}
		for (int q = 0; q < 50; ++ q) {
			ll x = ll(rng() % 41u) - 20;
			auto it = h.lower_bound(x);
			ll got = -((*it)[0] * x + (*it)[1]);
			REQUIRE(got == brute_min(lines, x));
		}
	}
}

TEST_CASE("line_multiset: zero-slope mixed with sloped", "[cht]") {
	line_multiset<line<ll>> h;
	lines_t lines = {{0, 10}, {1, 0}, {-1, 0}, {0, -5}};
	for (auto& [a, b] : lines) h.add_line(a, b);
	for (ll x = -20; x <= 20; ++ x) REQUIRE(query_max(h, x) == brute_max(lines, x));
}

TEST_CASE("line_multiset: extreme x values", "[cht]") {
	line_multiset<line<ll>> h;
	lines_t lines = {{1, 0}, {-1, 0}, {2, -5}, {-3, 7}};
	for (auto& [a, b] : lines) h.add_line(a, b);
	for (ll x : {ll(-1000000), ll(-1000), ll(0), ll(1000), ll(1000000)}) {
		REQUIRE(query_max(h, x) == brute_max(lines, x));
	}
}
