// Benchmarks for ds/hash_map.hpp, ds/indexed_set.hpp,
// ds/van_emde_boas_tree.hpp, ds/order_statistic.hpp, ds/disjoint_set.hpp.
#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>

#include "ds/hash_map.hpp"
#include "ds/order_statistic.hpp"
#include "ds/indexed_set.hpp"
#include "ds/van_emde_boas_tree.hpp"
#include "ds/disjoint_set.hpp"

TEST_CASE("hash_map vs std::unordered_map (1e6 inserts + 1e6 lookups)", "[!benchmark][hash_map]") {
	constexpr int N = 1'000'000;
	std::mt19937_64 rng(7);
	std::vector<long long> keys(N), lookups(N);
	for (auto& k : keys) k = (long long)rng();
	for (auto& k : lookups) k = keys[rng() % unsigned(N)];

	BENCHMARK("hash_map<long long, int>: 1e6 insert + 1e6 lookup") {
		hash_map<long long, int> m;
		for (int i = 0; i < N; ++i) m[keys[i]] = i;
		long long acc = 0;
		for (auto& k : lookups) acc += m[k];
		return acc;
	};

	BENCHMARK("std::unordered_map<long long, int>: 1e6 insert + 1e6 lookup") {
		std::unordered_map<long long, int> m;
		m.reserve(N * 2);
		for (int i = 0; i < N; ++i) m[keys[i]] = i;
		long long acc = 0;
		for (auto& k : lookups) acc += m[k];
		return acc;
	};
}

TEST_CASE("order_statistic_set: 5e5 inserts + 5e5 order queries", "[!benchmark][order_statistic]") {
	constexpr int N = 500'000;
	std::mt19937_64 rng(11);
	std::vector<int> ins(N), q(N);
	for (auto& x : ins) x = int(rng() & 0x7fffffff);
	for (auto& x : q)   x = int(rng() % unsigned(N));

	BENCHMARK("order_statistic_set<int> (pb_ds)") {
		order_statistic_set<int> s;
		for (auto& x : ins) s.insert(x);
		long long acc = 0;
		for (auto& i : q) acc += *s.find_by_order(i);
		return acc;
	};
}

TEST_CASE("indexed_set vs std::set (1e6 inserts + 1e6 contains)", "[!benchmark][indexed_set]") {
	constexpr int N = 1'000'000;
	std::mt19937_64 rng(13);
	std::vector<int> ins(N), q(N);
	for (auto& x : ins) x = int(rng() % unsigned(N));
	for (auto& x : q)   x = int(rng() % unsigned(N));

	BENCHMARK("indexed_set (universe = 1e6, bit-trie)") {
		indexed_set s(N);
		for (auto& x : ins) s.insert(x);
		int acc = 0;
		for (auto& x : q) acc += int(s.find(x));
		return acc;
	};

	BENCHMARK("std::set<int>") {
		std::set<int> s;
		for (auto& x : ins) s.insert(x);
		int acc = 0;
		for (auto& x : q) acc += int(s.count(x));
		return acc;
	};
}

TEST_CASE("van_emde_boas_tree: 5e5 inserts + 5e5 next-queries (U=2^20)", "[!benchmark][vEB]") {
	constexpr int N = 500'000;
	constexpr uint64_t U = 1ULL << 20;
	std::mt19937_64 rng(17);
	std::vector<int> ins(N), q(N);
	for (auto& x : ins) x = int(rng() % U);
	for (auto& x : q)   x = int(rng() % U);

	BENCHMARK("van_emde_boas_tree<20>") {
		van_emde_boas_tree<20> s;
		for (auto& x : ins) s.insert(x);
		long long acc = 0;
		for (auto& x : q) acc += int64_t(s.find_next(x));
		return acc;
	};
}

TEST_CASE("disjoint_set_size: 1e6 nodes, 1e6 random unions", "[!benchmark][dsu]") {
	constexpr int N = 1'000'000;
	constexpr int Q = 1'000'000;
	std::mt19937_64 rng(19);
	std::vector<std::pair<int, int>> ops(Q);
	for (auto& [a, b] : ops) {
		a = int(rng() % unsigned(N));
		b = int(rng() % unsigned(N));
	}
	BENCHMARK("1e6 merges over 1e6 nodes") {
		disjoint_set_size dsu(N);
		for (auto& [a, b] : ops) dsu.merge(a, b);
		return dsu.size();
	};
}
