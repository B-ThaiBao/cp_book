#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "misc/reverse_args.hpp"
#include "ds/swag.hpp"

namespace {
auto sum_f = [](const long long& a, const long long& b) { return a + b; };
auto min_f = [](const int& a, const int& b) { return std::min(a, b); };
auto cat_f = [](const std::string& a, const std::string& b) { return a + b; };
} // namespace

TEST_CASE("swag_stack: aggregate sum", "[swag]") {
	auto st = make_swag_stack<long long, long long>(sum_f);
	st.push_back(3);
	st.push_back(1);
	st.push_back(4);
	REQUIRE(st.back().second == 8);
	st.push_back(1);
	st.push_back(5);
	REQUIRE(st.back().second == 14);
}

TEST_CASE("swag_stack: aggregate min", "[swag]") {
	auto st = make_swag_stack<int, int>(min_f);
	st.push_back(5);
	st.push_back(2);
	st.push_back(8);
	REQUIRE(st.back().second == 2);
}

TEST_CASE("swag_queue: push_back + pop_front sum", "[swag]") {
	auto q = make_swag_queue<long long, long long>(sum_f);
	q.push_back(1); q.push_back(2); q.push_back(3);
	REQUIRE(q.front().second == 6);
	q.pop_front();
	REQUIRE(q.front().second == 5);
	q.push_back(10);
	REQUIRE(q.front().second == 15);
	q.pop_front();
	REQUIRE(q.front().second == 13);
}

TEST_CASE("swag_queue: string concatenation order", "[swag]") {
	auto q = make_swag_queue<std::string, std::string>(cat_f);
	q.push_back("a"); q.push_back("b"); q.push_back("c"); q.push_back("d");
	REQUIRE(q.front().second == "abcd");
	q.pop_front();
	REQUIRE(q.front().second == "bcd");
}

TEST_CASE("swag_deque: push_front + push_back sum", "[swag]") {
	auto d = make_swag_deque<long long, long long>(sum_f);
	d.push_back(1); d.push_back(2);
	d.push_front(10); d.push_front(20);
	REQUIRE(d.front().second == 33);
	d.pop_back();
	REQUIRE(d.front().second == 31);
	d.pop_front();
	REQUIRE(d.front().second == 11);
}

TEST_CASE("swag_deque: random ops stress (min)", "[swag][stress]") {
	std::mt19937 rng(91);
	auto d = make_swag_deque<int, int>(min_f);
	std::deque<int> ref;
	for (int t = 0; t < 2000; ++ t) {
		int op = int(rng() % 4);
		if (op == 0) { int v = int(rng() % 1000); d.push_back(v); ref.push_back(v); }
		else if (op == 1) { int v = int(rng() % 1000); d.push_front(v); ref.push_front(v); }
		else if (op == 2 && !ref.empty()) { d.pop_back(); ref.pop_back(); }
		else if (op == 3 && !ref.empty()) { d.pop_front(); ref.pop_front(); }
		if (!ref.empty()) {
			int expect = *std::min_element(ref.begin(), ref.end());
			REQUIRE(d.front().second == expect);
			REQUIRE(d.back().second == expect);
		}
	}
}

TEST_CASE("swag_queue: empty after pop", "[swag]") {
	auto q = make_swag_queue<int, int>(sum_f);
	q.push_back(5);
	REQUIRE(q.size() == 1);
	q.pop_front();
	REQUIRE(q.empty());
}
