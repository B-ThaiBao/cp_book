#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "math/matrix.hpp"

namespace {
using M22 = matrix<int, 2, 2>;
using M23 = matrix<int, 2, 3>;
using M32 = matrix<int, 3, 2>;
using M33 = matrix<long long, 3, 3>;
} // namespace

TEST_CASE("matrix: initializer_list and access", "[matrix]") {
	M22 a{1, 2, 3, 4};
	REQUIRE(a(0, 0) == 1);
	REQUIRE(a(0, 1) == 2);
	REQUIRE(a(1, 0) == 3);
	REQUIRE(a(1, 1) == 4);
	REQUIRE(a[(std::array<int, 2>{0, 1})] == 2);
}

TEST_CASE("matrix: add / sub / unary minus", "[matrix]") {
	M22 a{1, 2, 3, 4}, b{5, 6, 7, 8};
	M22 c = a + b;
	REQUIRE(c(0, 0) == 6);
	REQUIRE(c(1, 1) == 12);
	M22 d = b - a;
	REQUIRE(d(0, 0) == 4);
	REQUIRE(d(1, 1) == 4);
	M22 e = - a;
	REQUIRE(e(0, 0) == - 1);
	REQUIRE(e(1, 1) == - 4);
}

TEST_CASE("matrix: multiplication 2x3 * 3x2", "[matrix]") {
	M23 a{1, 2, 3, 4, 5, 6};
	M32 b{7, 8, 9, 10, 11, 12};
	matrix<int, 2, 2> c = a * b;
	REQUIRE(c(0, 0) == 58);
	REQUIRE(c(0, 1) == 64);
	REQUIRE(c(1, 0) == 139);
	REQUIRE(c(1, 1) == 154);
}

TEST_CASE("matrix: scalar add adds to diagonal", "[matrix]") {
	M22 a{1, 2, 3, 4};
	a += 10;
	REQUIRE(a(0, 0) == 11);
	REQUIRE(a(0, 1) == 2);
	REQUIRE(a(1, 0) == 3);
	REQUIRE(a(1, 1) == 14);
}

TEST_CASE("matrix: scalar sub subtracts from diagonal", "[matrix]") {
	M22 a{10, 2, 3, 20};
	a -= 5;
	REQUIRE(a(0, 0) == 5);
	REQUIRE(a(1, 1) == 15);
	REQUIRE(a(0, 1) == 2);
}

TEST_CASE("matrix: bin_pow on Fibonacci matrix", "[matrix][bin_pow]") {
	matrix<long long, 2, 2> F{1, 1, 1, 0};
	auto F10 = bin_pow(F, 10LL);
	REQUIRE(F10(0, 0) == 89);
	REQUIRE(F10(0, 1) == 55);
	REQUIRE(F10(1, 0) == 55);
	REQUIRE(F10(1, 1) == 34);
	auto F20 = bin_pow(F, 20LL);
	REQUIRE(F20(0, 1) == 6765);
}

TEST_CASE("matrix: bin_pow identity for k=0", "[matrix]") {
	matrix<long long, 2, 2> F{2, 3, 5, 7};
	auto I = bin_pow(F, 0LL);
	REQUIRE(I(0, 0) == 1);
	REQUIRE(I(0, 1) == 0);
	REQUIRE(I(1, 0) == 0);
	REQUIRE(I(1, 1) == 1);
}

TEST_CASE("matrix: transpose (square)", "[matrix]") {
	matrix<int, 3, 3> a{1, 2, 3, 4, 5, 6, 7, 8, 9};
	auto at = tranpose(a);
	REQUIRE(at(0, 0) == 1);
	REQUIRE(at(0, 1) == 4);
	REQUIRE(at(0, 2) == 7);
	REQUIRE(at(1, 0) == 2);
	REQUIRE(at(2, 1) == 6);
}

TEST_CASE("matrix: multiplication associativity stress", "[matrix][stress]") {
	std::mt19937 rng(2024);
	for (int t = 0; t < 20; ++ t) {
		matrix<int, 3, 3> a, b, c;
		for (int i = 0; i < 3; ++ i) for (int j = 0; j < 3; ++ j) {
			a(i, j) = int(rng() % 10u);
			b(i, j) = int(rng() % 10u);
			c(i, j) = int(rng() % 10u);
		}
		auto ab_c = (a * b) * c;
		auto a_bc = a * (b * c);
		for (int i = 0; i < 3; ++ i) for (int j = 0; j < 3; ++ j)
			REQUIRE(ab_c(i, j) == a_bc(i, j));
	}
}
