// Benchmarks for flow/* (dinic vs hlpp on dense bipartite & random graphs).
#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>

#include "flow/flow_graph.hpp"
#include "flow/max_flow.hpp"

namespace {
struct GraphSpec {
	int N;
	std::vector<std::tuple<int, int, long long>> edges;
};

GraphSpec random_dense_dag(int N, int avg_deg, uint64_t seed) {
	GraphSpec g; g.N = N;
	std::mt19937_64 rng(seed);
	for (int u = 0; u < N; ++u) {
		int d = avg_deg;
		for (int k = 0; k < d; ++k) {
			int v = u + 1 + int(rng() % unsigned(std::max(1, N - u - 1)));
			if (v >= N) continue;
			g.edges.push_back({u, v, (long long)(1 + (rng() % 100ULL))});
		}
	}
	return g;
}

flow_graph<long long> build(const GraphSpec& spec) {
	flow_graph<long long> g(spec.N);
	for (auto& [u, v, c] : spec.edges) g.add_edge(u, v, c);
	return g;
}
}

TEST_CASE("max_flow: dinic vs hlpp on dense DAG (N=400, avg_deg=20)", "[!benchmark][maxflow]") {
	auto spec = random_dense_dag(400, 20, 2024);

	BENCHMARK("dinic_max_flow N=400") {
		auto g = build(spec);
		dinic_max_flow<long long> mf(spec.N);
		return mf.max_flow(g, 0, spec.N - 1);
	};

	BENCHMARK("hlpp_max_flow N=400") {
		auto g = build(spec);
		hlpp_max_flow<long long> mf(spec.N);
		return mf.max_flow(g, 0, spec.N - 1);
	};
}

TEST_CASE("max_flow: dinic vs hlpp on sparse random (N=2000, avg_deg=5)", "[!benchmark][maxflow]") {
	auto spec = random_dense_dag(2000, 5, 4242);

	BENCHMARK("dinic_max_flow N=2000") {
		auto g = build(spec);
		dinic_max_flow<long long> mf(spec.N);
		return mf.max_flow(g, 0, spec.N - 1);
	};

	BENCHMARK("hlpp_max_flow N=2000") {
		auto g = build(spec);
		hlpp_max_flow<long long> mf(spec.N);
		return mf.max_flow(g, 0, spec.N - 1);
	};
}
