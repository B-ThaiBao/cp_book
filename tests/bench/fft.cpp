// Benchmarks for fft/* (naive vs FFT multiplication).
#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>

#include "fft/fft.hpp"

namespace {
std::vector<long long> make_poly(int N, uint64_t seed) {
	std::mt19937_64 rng(seed);
	std::vector<long long> a(N);
	for (auto& x : a) x = (long long)(rng() % 1000ULL);
	return a;
}
}

TEST_CASE("fft: naive vs fft_complex (N=M=4096)", "[!benchmark][fft]") {
	auto a = make_poly(4096, 1);
	auto b = make_poly(4096, 2);

	BENCHMARK("naive_multiply 4096 x 4096") {
		auto r = fft::naive_multiply<long long>(a, b);
		return r.size();
	};

	BENCHMARK("fft_complex_multiply 4096 x 4096") {
		auto r = fft::fft_complex_multiply<long long>(a, b);
		return r.size();
	};
}

TEST_CASE("fft: fft_complex_multiply at N=M=2^18", "[!benchmark][fft]") {
	constexpr int N = 1 << 18;
	auto a = make_poly(N, 3);
	auto b = make_poly(N, 4);

	BENCHMARK("fft_complex_multiply 2^18 x 2^18") {
		auto r = fft::fft_complex_multiply<long long>(a, b);
		return r.size();
	};
}
