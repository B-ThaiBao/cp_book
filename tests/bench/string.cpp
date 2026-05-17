// Benchmarks for string/* (suffix_array, z_function, prefix_function, manacher, hash).
#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>

#include "ds/range_min_query.hpp"
#include "string/suffix_array.hpp"
#include "string/z_function.hpp"
#include "string/prefix_function.hpp"
#include "string/manacher.hpp"
#include "string/prefix_power_hash.hpp"
#include "mod/modnum.hpp"

namespace {
std::vector<int> make_intstr(int N, int alpha = 26, uint64_t seed = 42) {
	std::mt19937_64 rng(seed);
	std::vector<int> s(N);
	for (auto& c : s) c = 1 + int(rng() % unsigned(alpha));
	return s;
}
std::string make_str(int N, int alpha = 4, uint64_t seed = 43) {
	std::mt19937_64 rng(seed);
	std::string s(N, 'a');
	for (auto& c : s) c = char('a' + int(rng() % unsigned(alpha)));
	return s;
}
}

TEST_CASE("suffix_array: N=1e6 over alphabet 26", "[!benchmark][suffix_array]") {
	auto s = make_intstr(1'000'000, 26);
	BENCHMARK("suffix_array(N=1e6)") {
		suffix_array sa(s);
		return sa.size();
	};
}

TEST_CASE("z_function: N=5e6", "[!benchmark][z_function]") {
	auto s = make_intstr(5'000'000, 26);
	BENCHMARK("z_function(N=5e6)") {
		auto z = z_function(s);
		return z.size();
	};
}

TEST_CASE("prefix_function: N=5e6", "[!benchmark][prefix_function]") {
	auto s = make_intstr(5'000'000, 26);
	BENCHMARK("prefix_function(N=5e6)") {
		auto p = prefix_function(s);
		return p.size();
	};
}

TEST_CASE("manacher: N=2e6", "[!benchmark][manacher]") {
	auto s = make_str(2'000'000, 4);
	BENCHMARK("manacher.build(N=2e6)") {
		manacher m;
		m.build(s);
		return m.size();
	};
}
