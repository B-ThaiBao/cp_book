#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "graph/twosat_graph.hpp"

namespace {
bool eval(const std::vector<bool>& res, const int& x) {
	if (x >= 0) return res[x];
	return !res[~ x];
}
} // namespace

TEST_CASE("twosat_graph: simple satisfiable", "[twosat]") {
	twosat_graph sat(3);
	sat.either(0, 1);
	sat.either(~ 1, 2);
	sat.single(0);
	auto res = sat.solve();
	REQUIRE_FALSE(res.empty());
	REQUIRE(res[0]);
	REQUIRE((res[0] || res[1]));
	REQUIRE((!res[1] || res[2]));
}

TEST_CASE("twosat_graph: unsatisfiable", "[twosat]") {
	twosat_graph sat(1);
	sat.single(0);
	sat.single(~ 0);
	auto res = sat.solve();
	REQUIRE(res.empty());
}

TEST_CASE("twosat_graph: xor + single", "[twosat]") {
	twosat_graph sat(2);
	sat.diff(0, 1);
	sat.single(0);
	auto res = sat.solve();
	REQUIRE_FALSE(res.empty());
	REQUIRE(res[0]);
	REQUIRE_FALSE(res[1]);
}

TEST_CASE("twosat_graph: implies chain", "[twosat]") {
	twosat_graph sat(3);
	sat.implies(0, 1);
	sat.implies(1, 2);
	sat.single(0);
	auto res = sat.solve();
	REQUIRE_FALSE(res.empty());
	REQUIRE(res[0]);
	REQUIRE(res[1]);
	REQUIRE(res[2]);
}

TEST_CASE("twosat_graph: same equivalence", "[twosat]") {
	twosat_graph sat(2);
	sat.same(0, 1);
	sat.single(0);
	auto res = sat.solve();
	REQUIRE_FALSE(res.empty());
	REQUIRE(res[0] == res[1]);
}

TEST_CASE("twosat_graph: no clauses trivially satisfiable", "[twosat]") {
	twosat_graph sat(3);
	auto res = sat.solve();
	REQUIRE(res.size() == 3);
}

TEST_CASE("twosat_graph: stress vs brute", "[twosat][stress]") {
	std::mt19937 rng(1234);
	for (int t = 0; t < 60; ++ t) {
		int N = 2 + int(rng() % 5u);
		int M = 1 + int(rng() % 8u);
		std::vector<std::pair<int, int>> clauses;
		twosat_graph sat(N);
		for (int i = 0; i < M; ++ i) {
			int a = int(rng() % unsigned(N)); if (rng() & 1) a = ~ a;
			int b = int(rng() % unsigned(N)); if (rng() & 1) b = ~ b;
			clauses.push_back({a, b});
			sat.either(a, b);
		}
		auto res = sat.solve();
		bool brute_sat = false;
		std::vector<bool> assign(N, false);
		for (int mask = 0; mask < (1 << N); ++ mask) {
			for (int i = 0; i < N; ++ i) assign[i] = (mask >> i) & 1;
			bool ok = true;
			for (const auto& cl : clauses) {
				auto val = [&](const int& x) -> bool {
					if (x >= 0) return assign[x];
					return !assign[~ x];
				};
				if (!val(cl.first) && !val(cl.second)) { ok = false; break; }
			}
			if (ok) { brute_sat = true; break; }
		}
		REQUIRE(brute_sat == !res.empty());
		if (!res.empty()) {
			for (const auto& cl : clauses)
				REQUIRE((eval(res, cl.first) || eval(res, cl.second)));
		}
	}
}
