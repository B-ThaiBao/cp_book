#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "misc/tensor.hpp"

namespace {
using A1 = std::array<int, 1>;
using A2 = std::array<int, 2>;
using A3 = std::array<int, 3>;
using A4 = std::array<int, 4>;
} // namespace

TEST_CASE("tensor: 1D construction & indexing", "[tensor]") {
	std::tensor<int, 1> t(A1{5}, 7);
	REQUIRE(t.len == 5);
	REQUIRE(t.shape[0] == 5);
	for (int i = 0; i < 5; ++ i) REQUIRE(t[A1{i}] == 7);
	t[A1{2}] = 42;
	REQUIRE(t[A1{2}] == 42);
	REQUIRE(t.at(A1{0}) == 7);
	REQUIRE(t.at(A1{4}) == 7);
}

TEST_CASE("tensor: 2D row-major and chained indexing", "[tensor]") {
	std::tensor<int, 2> t(A2{3, 4}, 0);
	REQUIRE(t.len == 12);
	REQUIRE(t.shape[0] == 3);
	REQUIRE(t.shape[1] == 4);
	REQUIRE(t.strides[0] == 4);
	REQUIRE(t.strides[1] == 1);

	for (int i = 0; i < 3; ++ i) {
		for (int j = 0; j < 4; ++ j) {
			t[A2{i, j}] = i * 10 + j;
		}
	}
	for (int i = 0; i < 3; ++ i) {
		for (int j = 0; j < 4; ++ j) {
			REQUIRE(t[A2{i, j}] == i * 10 + j);
			REQUIRE(t[i][A1{j}] == i * 10 + j);
		}
	}
}

TEST_CASE("tensor: 3D and view", "[tensor]") {
	std::tensor<int, 3> t(A3{2, 3, 4}, 1);
	REQUIRE(t.len == 24);
	int cnt = 0;
	for (int i = 0; i < 2; ++ i) for (int j = 0; j < 3; ++ j) for (int k = 0; k < 4; ++ k) {
		t[A3{i, j, k}] = cnt ++;
	}
	auto v = t.view();
	for (int i = 0; i < 2; ++ i) for (int j = 0; j < 3; ++ j) for (int k = 0; k < 4; ++ k) {
		REQUIRE(v[A3{i, j, k}] == t[A3{i, j, k}]);
	}
	const auto& ct = t;
	auto cv = ct.view();
	REQUIRE(cv[A3{1, 2, 3}] == t[A3{1, 2, 3}]);
}

TEST_CASE("tensor: copy & assign (deep)", "[tensor]") {
	std::tensor<int, 2> a(A2{2, 2}, 0);
	a[A2{0, 0}] = 1; a[A2{0, 1}] = 2; a[A2{1, 0}] = 3; a[A2{1, 1}] = 4;
	std::tensor<int, 2> b(a);
	REQUIRE(b[A2{1, 1}] == 4);
	b[A2{1, 1}] = 99;
	REQUIRE(a[A2{1, 1}] == 4);
	REQUIRE(b[A2{1, 1}] == 99);

	std::tensor<int, 2> c(A2{1, 1}, 0);
	c = a;
	REQUIRE(c[A2{0, 0}] == 1);
	REQUIRE(c.shape[0] == 2);
}

TEST_CASE("tensor: move", "[tensor]") {
	std::tensor<int, 1> a(A1{5}, 9);
	a[A1{0}] = 1;
	std::tensor<int, 1> b(std::move(a));
	REQUIRE(b.len == 5);
	REQUIRE(b[A1{0}] == 1);
	REQUIRE(b[A1{4}] == 9);
}

TEST_CASE("tensor: assign reshapes", "[tensor]") {
	std::tensor<int, 2> t(A2{2, 2}, 0);
	t.assign(A2{3, 5}, 7);
	REQUIRE(t.len == 15);
	REQUIRE(t.shape[0] == 3);
	REQUIRE(t.shape[1] == 5);
	REQUIRE(t[A2{2, 4}] == 7);
}

TEST_CASE("tensor: 1D large stress", "[tensor][stress]") {
	int N = 500;
	std::tensor<int, 1> t(A1{N}, 0);
	for (int i = 0; i < N; ++ i) t[A1{i}] = i * i;
	for (int i = 0; i < N; ++ i) REQUIRE(t[A1{i}] == i * i);
}

TEST_CASE("tensor: 2D stress vs naive 2D vector", "[tensor][stress]") {
	std::mt19937 rng(1234);
	int R = 25, C = 30;
	std::tensor<int, 2> t(A2{R, C}, 0);
	std::vector<std::vector<int>> ref(R, std::vector<int>(C, 0));
	for (int op = 0; op < 1000; ++ op) {
		int i = int(rng() % unsigned(R));
		int j = int(rng() % unsigned(C));
		int v = int(rng() % 1000u);
		t[A2{i, j}] = v;
		ref[i][j] = v;
	}
	for (int i = 0; i < R; ++ i) for (int j = 0; j < C; ++ j) REQUIRE(t[A2{i, j}] == ref[i][j]);
}

TEST_CASE("tensor: 4D shape & strides", "[tensor]") {
	std::tensor<int, 4> t(A4{2, 3, 4, 5}, 0);
	REQUIRE(t.len == 2 * 3 * 4 * 5);
	REQUIRE(t.strides[0] == 3 * 4 * 5);
	REQUIRE(t.strides[1] == 4 * 5);
	REQUIRE(t.strides[2] == 5);
	REQUIRE(t.strides[3] == 1);
	t[A4{1, 2, 3, 4}] = 99;
	REQUIRE(t[A4{1, 2, 3, 4}] == 99);
}

TEST_CASE("tensor: default ctor empty", "[tensor]") {
	std::tensor<int, 2> t;
	REQUIRE(t.len == 0);
	REQUIRE(t.data == nullptr);
}

TEST_CASE("tensor: self-assignment via copy", "[tensor]") {
	std::tensor<int, 1> a(A1{3}, 0);
	a[A1{0}] = 1; a[A1{1}] = 2; a[A1{2}] = 3;
	std::tensor<int, 1> b(A1{3}, 0);
	b = a;
	REQUIRE(b[A1{0}] == 1);
	REQUIRE(b[A1{1}] == 2);
	REQUIRE(b[A1{2}] == 3);
	b[A1{0}] = 99;
	REQUIRE(a[A1{0}] == 1);
}

TEST_CASE("tensor: 2D fill with view", "[tensor]") {
	std::tensor<int, 2> t(A2{4, 4}, 0);
	for (int i = 0; i < 4; ++ i) {
		auto row = t[i];
		for (int j = 0; j < 4; ++ j) row[A1{j}] = i * 4 + j;
	}
	for (int i = 0; i < 4; ++ i) for (int j = 0; j < 4; ++ j) {
		REQUIRE(t[A2{i, j}] == i * 4 + j);
	}
}
