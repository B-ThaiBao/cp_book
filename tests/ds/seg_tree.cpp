#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "ds/seg_tree.hpp"

namespace {
struct sum_seg {
	std::vector<long long> data;
	seg_tree::in_order_tree layout;
	sum_seg(const int& N) : data(2 * N, 0), layout(N) {}
	void point_update(const int& i, const long long& v) {
		seg_tree::point_t pt = layout.point(i);
		data[pt] = v;
		pt.for_ancestor_up([&](seg_tree::point_t p) {
			data[p] = data[p.c(0)] + data[p.c(1)];
		});
	}
	long long range_sum(const int& l, const int& r) const {
		long long s = 0;
		layout.range(l, r).for_each([&](seg_tree::point_t p) { s += data[p]; });
		return s;
	}
};

struct lazy_seg {
	std::vector<long long> sum, lazy;
	std::vector<int> sz;
	seg_tree::in_order_tree layout;
	lazy_seg(const int& N) : sum(2 * N, 0), lazy(2 * N, 0), sz(2 * N, 0), layout(N) {
		for (int i = 0; i < N; ++ i) sz[layout.point(i)] = 1;
		for (int i = N - 1; i >= 1; -- i) sz[i] = sz[i * 2] + sz[i * 2 + 1];
	}
	void apply_node(seg_tree::point_t p, long long add) {
		sum[p] += add * sz[p];
		lazy[p] += add;
	}
	void push(seg_tree::point_t p) {
		if (lazy[p] != 0) {
			apply_node(p.c(0), lazy[p]);
			apply_node(p.c(1), lazy[p]);
			lazy[p] = 0;
		}
	}
	void pull(seg_tree::point_t p) { sum[p] = sum[p.c(0)] + sum[p.c(1)]; }
	void range_add(const int& l, const int& r, const long long& v) {
		auto rng = layout.range(l, r);
		rng.for_ancestor_down([&](seg_tree::point_t p) { push(p); });
		rng.for_each([&](seg_tree::point_t p) { apply_node(p, v); });
		rng.for_ancestor_up([&](seg_tree::point_t p) { pull(p); });
	}
	long long range_sum(const int& l, const int& r) {
		auto rng = layout.range(l, r);
		rng.for_ancestor_down([&](seg_tree::point_t p) { push(p); });
		long long s = 0;
		rng.for_each([&](seg_tree::point_t p) { s += sum[p]; });
		return s;
	}
};
} // namespace

TEST_CASE("seg_tree: floor/ceil/next_pow_2", "[seg_tree]") {
	REQUIRE(seg_tree::floor_log_2(1) == 0);
	REQUIRE(seg_tree::floor_log_2(7) == 2);
	REQUIRE(seg_tree::floor_log_2(8) == 3);
	REQUIRE(seg_tree::ceil_log_2(1) == 0);
	REQUIRE(seg_tree::ceil_log_2(5) == 3);
	REQUIRE(seg_tree::ceil_log_2(8) == 3);
	REQUIRE(seg_tree::next_pow_2(1) == 1);
	REQUIRE(seg_tree::next_pow_2(5) == 8);
	REQUIRE(seg_tree::next_pow_2(16) == 16);
}

TEST_CASE("seg_tree::in_order_tree: leaf round-trip", "[seg_tree]") {
	for (int N : {1, 2, 5, 8, 11, 16, 31, 32, 100}) {
		seg_tree::in_order_tree t(N);
		for (int i = 0; i < N; ++ i) {
			auto pt = t.point(i);
			REQUIRE(t.is_leaf(pt));
			REQUIRE(t.leaf_index(pt) == i);
		}
	}
}

TEST_CASE("seg_tree: sum_seg point update + range sum", "[seg_tree]") {
	std::mt19937 rng(7);
	int N = 30;
	sum_seg st(N);
	std::vector<long long> ref(N, 0);
	for (int t = 0; t < 500; ++ t) {
		int op = rng() % 2;
		if (op == 0) {
			int i = int(rng() % unsigned(N));
			long long v = (long long)(rng() % 100) - 50;
			st.point_update(i, v);
			ref[i] = v;
		} else {
			int l = int(rng() % unsigned(N + 1)), r = int(rng() % unsigned(N + 1));
			if (l > r) std::swap(l, r);
			long long expect = 0;
			for (int j = l; j < r; ++ j) expect += ref[j];
			REQUIRE(st.range_sum(l, r) == expect);
		}
	}
}

TEST_CASE("seg_tree: lazy range_add + range_sum", "[seg_tree]") {
	std::mt19937 rng(13);
	int N = 25;
	lazy_seg st(N);
	std::vector<long long> ref(N, 0);
	for (int t = 0; t < 400; ++ t) {
		int op = rng() % 2;
		int l = int(rng() % unsigned(N + 1)), r = int(rng() % unsigned(N + 1));
		if (l > r) std::swap(l, r);
		if (op == 0) {
			long long v = (long long)(rng() % 20) - 10;
			st.range_add(l, r, v);
			for (int j = l; j < r; ++ j) ref[j] += v;
		} else {
			long long expect = 0;
			for (int j = l; j < r; ++ j) expect += ref[j];
			REQUIRE(st.range_sum(l, r) == expect);
		}
	}
}

TEST_CASE("seg_tree::circular_tree: leaf round-trip", "[seg_tree]") {
	int N = 16;
	seg_tree::circular_tree ct(N);
	for (int i = 0; i < N; ++ i) {
		auto pt = ct.point(i);
		REQUIRE(ct.is_leaf(pt));
		REQUIRE(ct.leaf_index(pt) == i);
	}
}

TEST_CASE("seg_tree: empty range_sum is zero", "[seg_tree]") {
	sum_seg st(5);
	st.point_update(2, 100);
	REQUIRE(st.range_sum(0, 0) == 0);
	REQUIRE(st.range_sum(3, 3) == 0);
	REQUIRE(st.range_sum(2, 3) == 100);
}
