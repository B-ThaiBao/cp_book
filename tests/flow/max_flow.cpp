#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "flow/flow_graph.hpp"
#include "flow/max_flow.hpp"

namespace {
struct brute_flow {
	int N;
	std::vector<std::vector<long long>> cap;
	brute_flow(const int& n) : N(n), cap(n, std::vector<long long>(n, 0)) {}
	void add_edge(const int& u, const int& v, const long long& c) { cap[u][v] += c; }
	long long max_flow(const int& s, const int& t) {
		long long ans = 0;
		while (true) {
			std::vector<int> par(N, - 1); par[s] = s;
			std::vector<long long> fbn(N, 0); fbn[s] = LLONG_MAX;
			std::queue<int> q; q.push(s);
			while (!q.empty() && par[t] == - 1) {
				int u = q.front(); q.pop();
				for (int v = 0; v < N; ++ v) {
					if (par[v] == - 1 && cap[u][v] > 0) {
						par[v] = u; fbn[v] = std::min(fbn[u], cap[u][v]);
						q.push(v);
					}
				}
			}
			if (par[t] == - 1) break;
			long long f = fbn[t];
			for (int v = t; v != s; v = par[v]) {
				cap[par[v]][v] -= f;
				cap[v][par[v]] += f;
			}
			ans += f;
		}
		return ans;
	}
};
} // namespace

TEST_CASE("dinic_max_flow: classic example", "[maxflow]") {
	flow_graph<long long> g(4);
	g.add_edge(0, 1, 3);
	g.add_edge(0, 2, 2);
	g.add_edge(1, 2, 2);
	g.add_edge(1, 3, 2);
	g.add_edge(2, 3, 3);
	dinic_max_flow<long long> mf(4);
	REQUIRE(mf.max_flow(g, 0, 3) == 5);
}

TEST_CASE("hlpp_max_flow: classic example", "[maxflow][hlpp]") {
	flow_graph<long long> g(4);
	g.add_edge(0, 1, 3);
	g.add_edge(0, 2, 2);
	g.add_edge(1, 2, 2);
	g.add_edge(1, 3, 2);
	g.add_edge(2, 3, 3);
	hlpp_max_flow<long long> mf(4);
	REQUIRE(mf.max_flow(g, 0, 3) == 5);
}

TEST_CASE("dinic_max_flow: direct path", "[maxflow]") {
	flow_graph<long long> g(2);
	g.add_edge(0, 1, 42);
	dinic_max_flow<long long> mf(2);
	REQUIRE(mf.max_flow(g, 0, 1) == 42);
}

TEST_CASE("dinic_max_flow: disconnected gives 0", "[maxflow]") {
	flow_graph<long long> g(4);
	g.add_edge(0, 1, 5);
	g.add_edge(2, 3, 5);
	dinic_max_flow<long long> mf(4);
	REQUIRE(mf.max_flow(g, 0, 3) == 0);
}

TEST_CASE("hlpp_max_flow: direct path", "[maxflow][hlpp]") {
	flow_graph<long long> g(2);
	g.add_edge(0, 1, 99);
	hlpp_max_flow<long long> mf(2);
	REQUIRE(mf.max_flow(g, 0, 1) == 99);
}

TEST_CASE("max_flow: dinic vs hlpp vs brute stress", "[maxflow][stress]") {
	std::mt19937 rng(31);
	for (int t = 0; t < 50; ++ t) {
		int N = 3 + int(rng() % 5u);
		flow_graph<long long> g1(N), g2(N);
		brute_flow ref(N);
		int M = int(rng() % 10u) + 1;
		for (int i = 0; i < M; ++ i) {
			int u = int(rng() % unsigned(N)), v = int(rng() % unsigned(N));
			if (u == v) continue;
			long long c = (long long)(rng() % 10u) + 1;
			g1.add_edge(u, v, c);
			g2.add_edge(u, v, c);
			ref.add_edge(u, v, c);
		}
		long long brute = ref.max_flow(0, N - 1);
		dinic_max_flow<long long> mf1(N);
		REQUIRE(mf1.max_flow(g1, 0, N - 1) == brute);
		hlpp_max_flow<long long> mf2(N);
		REQUIRE(mf2.max_flow(g2, 0, N - 1) == brute);
	}
}

TEST_CASE("max_flow: parallel edges", "[maxflow]") {
	flow_graph<long long> g(2);
	g.add_edge(0, 1, 3);
	g.add_edge(0, 1, 7);
	g.add_edge(0, 1, 2);
	dinic_max_flow<long long> mf(2);
	REQUIRE(mf.max_flow(g, 0, 1) == 12);
}
