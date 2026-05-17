#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "misc/reverse_args.hpp"

TEST_CASE("reverse_args: subtraction order swapped", "[reverse_args]") {
	auto sub = [](int a, int b) { return a - b; };
	auto rsub = std::reverse_args(sub);
	REQUIRE(sub(10, 3) == 7);
	REQUIRE(rsub(10, 3) == -7);
	REQUIRE(rsub(3, 10) == 7);
	REQUIRE(rsub(0, 0) == 0);
}

TEST_CASE("reverse_args: string concat", "[reverse_args]") {
	auto cat = [](const std::string& a, const std::string& b) { return a + b; };
	auto rcat = std::reverse_args(cat);
	REQUIRE(cat("a", "b") == "ab");
	REQUIRE(rcat("a", "b") == "ba");
	REQUIRE(rcat("", "x") == "x");
	REQUIRE(rcat("hello", "world") == "worldhello");
}

TEST_CASE("reverse_args: division integer", "[reverse_args]") {
	auto divf = [](int a, int b) { return a / b; };
	auto rdiv = std::reverse_args(divf);
	REQUIRE(rdiv(2, 10) == 5);
	REQUIRE(rdiv(3, 100) == 33);
	REQUIRE(rdiv(7, 14) == 2);
}

TEST_CASE("reverse_args: works with stored functor through const ref", "[reverse_args]") {
	auto add = [](int a, int b) { return a * 10 + b; };
	const auto radd = std::reverse_args(add);
	REQUIRE(radd(1, 2) == 21);
	REQUIRE(radd(5, 0) == 5);
	REQUIRE(radd(9, 9) == 99);
}

TEST_CASE("reverse_args: rvalue invocation", "[reverse_args]") {
	auto sub = [](int a, int b) { return a - b; };
	REQUIRE(std::reverse_args(sub)(2, 5) == 3);
	REQUIRE(std::reverse_args([](int a, int b) { return a * 2 + b; })(1, 10) == 21);
}

TEST_CASE("reverse_args: with std functions", "[reverse_args]") {
	auto pwr = std::reverse_args([](int e, int b) { return int(std::pow(b, e)); });
	// pwr(b, e) wraps lambda(e, b) → pow(b, e)
	REQUIRE(pwr(3, 2) == 9);   // pow(3, 2)
	REQUIRE(pwr(2, 5) == 32);  // pow(2, 5)
}

TEST_CASE("reverse_args: pair construction order", "[reverse_args]") {
	auto mp = [](int a, int b) { return std::make_pair(a, b); };
	auto rmp = std::reverse_args(mp);
	auto p = rmp(1, 2);
	REQUIRE(p.first == 2);
	REQUIRE(p.second == 1);
}

TEST_CASE("reverse_args: void return type", "[reverse_args]") {
	std::vector<int> v;
	auto push = std::reverse_args([&](int b, int a) { v.push_back(a); v.push_back(b); });
	// rpush(a, b) -> push(b, a): pushes a then b
	push(1, 2);
	push(3, 4);
	REQUIRE(v == std::vector<int>({1, 2, 3, 4}));
}

TEST_CASE("reverse_args: chaining reverse twice is identity", "[reverse_args]") {
	auto sub = [](int a, int b) { return a - b; };
	auto rr = std::reverse_args(std::reverse_args(sub));
	for (int a = -5; a <= 5; ++ a) for (int b = -5; b <= 5; ++ b) {
		REQUIRE(rr(a, b) == sub(a, b));
	}
}
