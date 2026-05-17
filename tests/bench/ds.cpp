// Benchmarks for ds/binary_indexed_tree.hpp, ds/sparse_table.hpp,
// ds/range_min_query.hpp, ds/seg_tree.hpp.
//
// Run with:
//   ./build/tests/bench/ds [!benchmark]
#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>

#include "ds/binary_indexed_tree.hpp"
#include "ds/sparse_table.hpp"
#include "ds/range_min_query.hpp"
#include "ds/seg_tree.hpp"

namespace {
std::vector<long long> make_data(int N, uint64_t seed = 0xC0DEC0DEULL) {
	std::mt19937_64 rng(seed);
	std::vector<long long> v(N);
	for (auto& x : v) x = (long long)(rng() & ((1ULL << 30) - 1));
	return v;
}
std::vector<std::pair<int, int>> make_ranges(int Q, int N, uint64_t seed = 0xBEEFULL) {
	std::mt19937_64 rng(seed);
	std::vector<std::pair<int, int>> q(Q);
	for (auto& [l, r] : q) {
		l = int(rng() % unsigned(N));
		r = int(rng() % unsigned(N));
		if (l > r) std::swap(l, r);
		++r;
	}
	return q;
}

struct SumSeg {
	seg_tree::in_order_tree layout;
	std::vector<long long> sum;
	SumSeg(int N) : layout(N), sum(2 * N, 0) {}
	void update(int i, long long v) {
		auto p = layout.point(i);
		sum[int(p)] = v;
		p.for_ancestor_up([&](seg_tree::point_t pt) {
			sum[int(pt)] = sum[int(pt.c(0))] + sum[int(pt.c(1))];
		});
	}
	long long query(int l, int r) {
		long long res = 0;
		layout.range(l, r).for_each([&](seg_tree::point_t pt) {
			res += sum[int(pt)];
		});
		return res;
	}
};
}

TEST_CASE("seg_tree (in_order, sum): N=1e5, 1e5 updates + 1e5 queries", "[!benchmark][seg_tree]") {
	constexpr int N = 100'000;
	constexpr int Q = 100'000;
	auto data = make_data(N);
	auto ranges = make_ranges(Q, N);

	BENCHMARK("build + 1e5 updates + 1e5 sum queries") {
		SumSeg seg(N);
		for (int i = 0; i < N; ++i) seg.update(i, data[i]);
		long long acc = 0;
		for (auto& [l, r] : ranges) acc += seg.query(l, r);
		return acc;
	};
}

TEST_CASE("binary_indexed_tree: 1e6 point_update + 1e6 prefix_query", "[!benchmark][bit]") {
	constexpr int N = 1'000'000;
	constexpr int Q = 1'000'000;
	auto data = make_data(N);
	std::mt19937_64 rng(123);
	std::vector<int> idx(Q);
	for (auto& x : idx) x = int(rng() % unsigned(N));

	BENCHMARK("BIT point_update + prefix_query (N=1e6)") {
		binary_indexed_tree<long long> bit(N, 0);
		for (int i = 0; i < N; ++i)
			for (auto& x : bit.suffix(i)) x += data[i];
		long long acc = 0;
		for (auto& i : idx) {
			long long s = 0;
			for (auto& x : bit.prefix(i)) s += x;
			acc += s;
		}
		return acc;
	};
}

TEST_CASE("sparse_table: build (N=1e6) + 1e6 range_min queries", "[!benchmark][sparse_table]") {
	constexpr int N = 1'000'000;
	constexpr int Q = 1'000'000;
	auto data = make_data(N);
	auto ranges = make_ranges(Q, N);

	BENCHMARK("sparse_table<long long> min (build + 1e6 queries)") {
		sparse_table<long long> st(N);
		for (int i = 0; i < N; ++i) st(0, i) = data[i];
		st.build([&](auto p) {
			st[p] = std::min(st[p.c(0)], st[p.c(1)]);
		});
		long long acc = 0;
		for (auto& [l, r] : ranges) {
			auto rg = st.range(l, r);
			acc += std::min(st[rg[0]], st[rg[1]]);
		}
		return acc;
	};
}

TEST_CASE("range_min_query: build + 1e6 queries (mask + sparse)", "[!benchmark][rmq]") {
	constexpr int N = 1'000'000;
	constexpr int Q = 1'000'000;
	auto raw = make_data(N);
	std::vector<int> data(N);
	for (int i = 0; i < N; ++i) data[i] = int(raw[i] & 0x7fffffff);
	auto ranges = make_ranges(Q, N);

	BENCHMARK("range_min_query<int> (N=1e6, 1e6 queries)") {
		range_min_query<int> rmq(data);
		long long acc = 0;
		for (auto& [l, r] : ranges) {
			auto [idx, val] = rmq.range_query(l, r - 1);
			acc += val;
		}
		return acc;
	};
}
