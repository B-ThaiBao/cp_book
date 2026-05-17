#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "mod/modnum.hpp"
#include "fft/berlekamp_massey.hpp"

namespace {
using num = modnum<constant<int, 1000000007>, naive_multiplier<int>>;

static std::vector<num> extend(const std::vector<num>& A, const std::vector<num>& C, int target_len) {
	std::vector<num> B = A;
	while (int(B.size()) < target_len) {
		num s(0);
		for (int j = 0; j < int(C.size()); ++ j) s += C[j] * B[B.size() - 1 - j];
		B.push_back(s);
	}
	return B;
}
} // namespace

TEST_CASE("berlekamp_massey: fibonacci", "[bm]") {
	std::vector<num> A = {num(1), num(1), num(2), num(3), num(5), num(8), num(13), num(21)};
	auto C = berlekamp_massey(A);
	REQUIRE(C.size() == 2);
	REQUIRE(C[0] == num(1));
	REQUIRE(C[1] == num(1));
}

TEST_CASE("berlekamp_massey: geometric a[i] = 2*a[i-1]", "[bm]") {
	std::vector<num> A;
	num x(3);
	for (int i = 0; i < 6; ++ i) { A.push_back(x); x = x * num(2); }
	auto C = berlekamp_massey(A);
	REQUIRE(C.size() == 1);
	REQUIRE(C[0] == num(2));
}

TEST_CASE("berlekamp_massey: empty recurrence for all zeros", "[bm]") {
	std::vector<num> A = {num(0), num(0), num(0), num(0), num(0)};
	auto C = berlekamp_massey(A);
	REQUIRE(C.empty());
}

TEST_CASE("berlekamp_massey: constant sequence", "[bm]") {
	std::vector<num> A(8, num(5));
	auto C = berlekamp_massey(A);
	REQUIRE(C.size() == 1);
	REQUIRE(C[0] == num(1));
}

TEST_CASE("berlekamp_massey: tribonacci", "[bm]") {
	std::vector<num> A = {num(1), num(1), num(2)};
	for (int i = 3; i < 12; ++ i) A.push_back(A[i - 1] + A[i - 2] + A[i - 3]);
	auto C = berlekamp_massey(A);
	REQUIRE(C.size() == 3);
	REQUIRE(C[0] == num(1));
	REQUIRE(C[1] == num(1));
	REQUIRE(C[2] == num(1));
}

TEST_CASE("berlekamp_massey: predict next term (fibonacci)", "[bm]") {
	std::vector<num> A;
	num a(1), b(1);
	for (int i = 0; i < 10; ++ i) { A.push_back(a); num t = a + b; a = b; b = t; }
	auto C = berlekamp_massey(A);
	auto B = extend(A, C, 20);
	std::vector<num> ref;
	num x(1), y(1);
	for (int i = 0; i < 20; ++ i) { ref.push_back(x); num t = x + y; x = y; y = t; }
	for (int i = 0; i < 20; ++ i) REQUIRE(B[i] == ref[i]);
}

TEST_CASE("berlekamp_massey: random linear recurrence", "[bm][stress]") {
	std::mt19937 rng(2024);
	for (int t = 0; t < 30; ++ t) {
		int order = 1 + int(rng() % 5u);
		std::vector<num> C(order);
		for (auto& c : C) c = num(int(rng() % 100u) + 1);
		std::vector<num> A(order);
		for (auto& a : A) a = num(int(rng() % 100u) + 1);
		for (int i = order; i < 4 * order + 5; ++ i) {
			num s(0);
			for (int j = 0; j < order; ++ j) s += C[j] * A[i - 1 - j];
			A.push_back(s);
		}
		auto C2 = berlekamp_massey(A);
		REQUIRE(int(C2.size()) <= order);
		auto B = extend(A, C2, int(A.size()) + 5);
		for (int i = 0; i < int(A.size()); ++ i) REQUIRE(B[i] == A[i]);
	}
}

TEST_CASE("berlekamp_massey: alternating signs", "[bm]") {
	std::vector<num> A;
	for (int i = 0; i < 8; ++ i) A.push_back(num(i % 2 ? -3 : 3));
	auto C = berlekamp_massey(A);
	REQUIRE(C.size() == 1);
	REQUIRE(C[0] == num(- 1));
}
