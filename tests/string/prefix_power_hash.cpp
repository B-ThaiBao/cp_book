#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "mod/modnum.hpp"
#include "string/prefix_power_hash.hpp"

namespace {
struct hash_init {
	hash_init() { base_power<hashnum>::base = hashnum(131); }
};
static hash_init _init;
} // namespace

TEST_CASE("hashnum: basic arithmetic", "[hashnum]") {
	hashnum a(123), b(456);
	REQUIRE(a + b == hashnum(579));
	REQUIRE(b - a == hashnum(333));
	REQUIRE(a * b == hashnum(123 * 456));
	REQUIRE(a == a);
	REQUIRE_FALSE(a == b);
}

TEST_CASE("hashnum: identities", "[hashnum]") {
	hashnum z(0), o(1), x(42);
	REQUIRE(x + z == x);
	REQUIRE(x * o == x);
	REQUIRE(x - x == z);
}

TEST_CASE("prefix_power_hash: equal substrings have equal hash", "[hash]") {
	prefix_power_hash<hashnum> H(std::string("abcabcabc"));
	REQUIRE(H.range_hash(0, 3) == H.range_hash(3, 6));
	REQUIRE(H.range_hash(0, 3) == H.range_hash(6, 9));
	REQUIRE(H.range_hash(0, 6) == H.range_hash(3, 9));
	REQUIRE_FALSE(H.range_hash(0, 4) == H.range_hash(1, 5));
}

TEST_CASE("prefix_power_hash: empty range", "[hash]") {
	prefix_power_hash<hashnum> H(std::string("abcdef"));
	REQUIRE(H.range_hash(2, 2) == H.range_hash(4, 4));
}

TEST_CASE("prefix_power_hash: find_lcp basic", "[hash]") {
	prefix_power_hash<hashnum> A(std::string("abcdefgh"));
	prefix_power_hash<hashnum> B(std::string("abcdxyzh"));
	REQUIRE(find_lcp(A, 0, B, 0) == 4);
	REQUIRE(find_lcp(A, 5, B, 5) == 0);
}

TEST_CASE("prefix_power_hash: find_lcp full match", "[hash]") {
	prefix_power_hash<hashnum> A(std::string("hello"));
	prefix_power_hash<hashnum> B(std::string("hello"));
	REQUIRE(find_lcp(A, 0, B, 0) == 5);
}

TEST_CASE("prefix_power_hash: find_lcp stress", "[hash][stress]") {
	std::mt19937 rng(101);
	for (int t = 0; t < 30; ++ t) {
		int N = 5 + int(rng() % 30u);
		std::string s;
		for (int i = 0; i < N; ++ i) s += char('a' + int(rng() % 3u));
		prefix_power_hash<hashnum> H(s);
		for (int a = 0; a < N; ++ a)
			for (int b = 0; b < N; ++ b) {
				int got = find_lcp(H, a, H, b);
				int naive = 0;
				while (a + naive < N && b + naive < N && s[a + naive] == s[b + naive]) ++ naive;
				REQUIRE(got == naive);
			}
	}
}

TEST_CASE("prefix_power_hash: compare equal-length", "[hash]") {
	prefix_power_hash<hashnum> A(std::string("apple"));
	prefix_power_hash<hashnum> B(std::string("apply"));
	REQUIRE(compare(A, 0, 5, B, 0, 5) == - 1);
	REQUIRE(compare(B, 0, 5, A, 0, 5) == 1);
	REQUIRE(compare(A, 0, 5, A, 0, 5) == 0);
}

TEST_CASE("prefix_power_hash: push_back appends", "[hash]") {
	prefix_power_hash<hashnum> A(std::string("abc"));
	A.push_back('d');
	A.push_back('e');
	prefix_power_hash<hashnum> B(std::string("abcde"));
	REQUIRE(A.range_hash(0, 5) == B.range_hash(0, 5));
}

TEST_CASE("prefix_power_hash: substring search stress", "[hash][stress]") {
	std::mt19937 rng(2024);
	for (int t = 0; t < 20; ++ t) {
		int N = 10 + int(rng() % 30u);
		int M = 1 + int(rng() % 4u);
		std::string s, p;
		for (int i = 0; i < N; ++ i) s += char('a' + int(rng() % 2u));
		for (int i = 0; i < M; ++ i) p += char('a' + int(rng() % 2u));
		prefix_power_hash<hashnum> SH(s), PH(p);
		auto target = PH.range_hash(0, M);
		std::vector<int> got, ref;
		for (int i = 0; i + M <= N; ++ i) {
			if (SH.range_hash(i, i + M) == target) got.push_back(i);
			if (s.compare(i, M, p) == 0) ref.push_back(i);
		}
		REQUIRE(got == ref);
	}
}
