#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "ds/hilbert_curve.hpp"

namespace {
struct query_t { int L, R, id; };
} // namespace

TEST_CASE("hilbert_order: deterministic", "[hilbert]") {
	uint64_t a = hilbert_curve::hilbert_order(0, 0);
	uint64_t b = hilbert_curve::hilbert_order(0, 0);
	REQUIRE(a == b);
	REQUIRE(hilbert_curve::hilbert_order(0, 0) != hilbert_curve::hilbert_order(1, 0));
}

TEST_CASE("hilbert_order: bijection on small grid", "[hilbert]") {
	int M = 16;
	std::set<uint64_t> seen;
	for (int x = 0; x < M; ++ x)
		for (int y = 0; y < M; ++ y)
			REQUIRE(seen.insert(hilbert_curve::hilbert_order(x, y)).second);
	REQUIRE((int)seen.size() == M * M);
}

TEST_CASE("hilbert_order: distinct across larger grid", "[hilbert]") {
	int M = 32;
	std::set<uint64_t> seen;
	for (int x = 0; x < M; ++ x)
		for (int y = 0; y < M; ++ y)
			seen.insert(hilbert_curve::hilbert_order(x, y));
	REQUIRE((int)seen.size() == M * M);
}

TEST_CASE("hilbert_curve::for_each: Mo sweep correctness", "[hilbert]") {
	int N = 50;
	std::vector<query_t> queries;
	std::mt19937 rng(1);
	int Qn = 30;
	for (int i = 0; i < Qn; ++ i) {
		int l = int(rng()) % N, r = int(rng()) % N + 1;
		if (l >= r) std::swap(l, r), ++ r;
		if (r > N) r = N;
		if (l == r) continue;
		queries.push_back({l, r, (int)queries.size()});
	}
	auto sorted = queries;
	std::sort(sorted.begin(), sorted.end(), [](const query_t& a, const query_t& b) {
		return hilbert_curve::hilbert_order(a.L, a.R) < hilbert_curve::hilbert_order(b.L, b.R);
	});
	std::vector<int> A(N);
	for (int i = 0; i < N; ++ i) A[i] = i * i;
	long long current = 0;
	std::vector<long long> answers(queries.size());
	hilbert_curve::for_each(sorted.begin(), sorted.end(),
		[&](int i) { current += A[i]; },
		[&](int i) { current += A[i]; },
		[&](int i) { current -= A[i]; },
		[&](int i) { current -= A[i]; },
		[&](const query_t& q) { answers[q.id] = current; }
	);
	for (const auto& q : queries) {
		long long expect = 0;
		for (int i = q.L; i < q.R; ++ i) expect += A[i];
		REQUIRE(answers[q.id] == expect);
	}
}
