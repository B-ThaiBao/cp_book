#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "string/prefix_tree.hpp"

namespace {
static std::vector<int> str_to_alpha(const std::string& s) {
	std::vector<int> r; r.reserve(s.size());
	for (char c : s) r.push_back(c - 'a');
	return r;
}
} // namespace

TEST_CASE("prefix_tree: insert / find", "[trie]") {
	prefix_tree t(26);
	t.make_root();
	t.insert(str_to_alpha("apple"));
	t.insert(str_to_alpha("app"));
	t.insert(str_to_alpha("apply"));
	REQUIRE(t.find(str_to_alpha("app")) >= 0);
	REQUIRE(t.find(str_to_alpha("apple")) >= 0);
	REQUIRE(t.find(str_to_alpha("banana")) < 0);
}

TEST_CASE("prefix_tree: erase decrements num_ends", "[trie]") {
	prefix_tree t(26);
	t.make_root();
	t.insert(str_to_alpha("apple"));
	int n = t.find(str_to_alpha("apple"));
	REQUIRE(t[n].num_ends == 1);
	t.erase(str_to_alpha("apple"));
	int n2 = t.find(str_to_alpha("apple"));
	REQUIRE(t[n2].num_ends == 0);
}

TEST_CASE("prefix_tree: insert duplicates", "[trie]") {
	prefix_tree t(26);
	t.make_root();
	t.insert(str_to_alpha("foo"));
	t.insert(str_to_alpha("foo"));
	t.insert(str_to_alpha("foo"));
	int n = t.find(str_to_alpha("foo"));
	REQUIRE(t[n].num_ends == 3);
}

TEST_CASE("prefix_tree: count_prefix", "[trie]") {
	prefix_tree t(26);
	t.make_root();
	t.insert(str_to_alpha("a"));
	t.insert(str_to_alpha("ab"));
	t.insert(str_to_alpha("abc"));
	t.insert(str_to_alpha("abcd"));
	REQUIRE(t.count_prefix(str_to_alpha("abcz")) == 3);
	REQUIRE(t.count_prefix(str_to_alpha("abcd")) == 4);
}

TEST_CASE("prefix_tree: count_prefix no match", "[trie]") {
	prefix_tree t(26);
	t.make_root();
	t.insert(str_to_alpha("abc"));
	REQUIRE(t.count_prefix(str_to_alpha("xyz")) == 0);
}

TEST_CASE("prefix_tree: num_starts tracks string count with prefix", "[trie]") {
	prefix_tree t(26);
	t.make_root();
	t.insert(str_to_alpha("apple"));
	t.insert(str_to_alpha("app"));
	t.insert(str_to_alpha("application"));
	int n_app = t.find(str_to_alpha("app"));
	REQUIRE(t[n_app].num_starts == 3);
	int n_appl = t.find(str_to_alpha("appl"));
	REQUIRE(t[n_appl].num_starts == 2);
}

TEST_CASE("prefix_tree: empty string is root", "[trie]") {
	prefix_tree t(26);
	t.make_root();
	t.insert(str_to_alpha(""));
	int n = t.find(str_to_alpha(""));
	REQUIRE(n >= 0);
	REQUIRE(t[n].num_ends == 1);
}

TEST_CASE("prefix_tree: stress insert/find", "[trie][stress]") {
	std::mt19937 rng(77);
	prefix_tree t(4);
	t.make_root();
	std::map<std::string, int> ref;
	std::vector<std::string> all;
	for (int i = 0; i < 200; ++ i) {
		int L = 1 + int(rng() % 6u);
		std::string s;
		for (int j = 0; j < L; ++ j) s += char('a' + int(rng() % 4u));
		t.insert(str_to_alpha(s));
		++ ref[s];
		all.push_back(s);
	}
	for (auto& [s, cnt] : ref) {
		int n = t.find(str_to_alpha(s));
		REQUIRE(n >= 0);
		REQUIRE(t[n].num_ends == cnt);
	}
}
