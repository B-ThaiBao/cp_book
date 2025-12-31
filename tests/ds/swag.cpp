#include <catch2/catch_test_macros.hpp>
#include <vector>
#include <deque>
#include <string>
#include <algorithm>
#include <numeric>
#include <random>

#include "misc/reverse_args.hpp"
#include "ds/swag.hpp"

template <class P, class F>
static P naive_fold(const std::deque<P>& a, const F& f) {
	REQUIRE(!a.empty());
	P res = a[0];
	for (size_t i = 1; i < a.size(); ++i) res = f(res, a[i]);
	return res;
}

template <class T, class P, class F, class Q>
static void assert_queue_state(const Q& q, const std::deque<T>& ref, const F& f) {
	REQUIRE(q.empty() == ref.empty());
	REQUIRE(q.size() == ref.size());
	if (ref.empty()) return;

	std::deque<P> tmp(ref.begin(), ref.end());
	P whole = naive_fold<P>(tmp, f);

	auto [front_val, front_agg] = q.front();
	auto [back_val, back_agg]   = q.back();

	REQUIRE(front_val == ref.front());
	REQUIRE(back_val  == ref.back());
	REQUIRE(front_agg == whole);
	REQUIRE(back_agg  == whole);
}

template <class T, class P, class F, class DQ>
static void assert_deque_state(const DQ& dq, const std::deque<T>& ref, const F& f) {
	REQUIRE(dq.empty() == ref.empty());
	REQUIRE(dq.size() == ref.size());
	if (ref.empty()) return;

	std::deque<P> tmp(ref.begin(), ref.end());
	P whole = naive_fold<P>(tmp, f);

	auto [front_val, front_agg] = dq.front();
	auto [back_val, back_agg]   = dq.back();

	REQUIRE(front_val == ref.front());
	REQUIRE(back_val  == ref.back());
	REQUIRE(front_agg == whole);
	REQUIRE(back_agg  == whole);
}

TEST_CASE("swag_stack: sum aggregate", "[swag_stack]") {
	auto f = [](int a, int b) { return a + b; };
	auto st = make_swag_stack<int, int>(f);

	std::vector<int> ref;
	for (int i = 1; i <= 10; ++i) {
		st.push_back(i);
		ref.push_back(i);

		int expected = std::accumulate(ref.begin(), ref.end(), 0);
		REQUIRE(st.back().first == i);
		REQUIRE(st.back().second == expected);
	}
}

TEST_CASE("swag_stack: min aggregate", "[swag_stack]") {
	auto f = [](int a, int b) { return std::min(a, b); };
	auto st = make_swag_stack<int, int>(f);

	std::vector<int> ref;
	for (int v : {5, 3, 7, 2, 9, 2, 4}) {
		st.push_back(v);
		ref.push_back(v);

		int expected = *std::min_element(ref.begin(), ref.end());
		REQUIRE(st.back().second == expected);
	}
}

TEST_CASE("swag_queue: sum push/pop", "[swag_queue]") {
	auto f = [](int a, int b) { return a + b; };
	auto q = make_swag_queue<int, int>(f);
	std::deque<int> ref;

	for (int v : {1, 2, 3, 4, 5}) {
		q.push_back(v);
		ref.push_back(v);
		assert_queue_state<int,int>(q, ref, f);
	}

	for (int i = 0; i < 3; ++i) {
		q.pop_front();
		ref.pop_front();
		assert_queue_state<int,int>(q, ref, f);
	}

	for (int v : {10, 20, 30}) {
		q.push_back(v);
		ref.push_back(v);
		assert_queue_state<int,int>(q, ref, f);
	}

	while (!ref.empty()) {
		q.pop_front();
		ref.pop_front();
		assert_queue_state<int,int>(q, ref, f);
	}
}

TEST_CASE("swag_queue: preserves order for string concat", "[swag_queue]") {
	auto f = [](const std::string& a, const std::string& b) { return a + b; };
	auto q = make_swag_queue<std::string, std::string>(f);
	std::deque<std::string> ref;

	for (auto s : {"A", "B", "C", "D"}) {
		q.push_back(std::string(s));
		ref.push_back(std::string(s));
		assert_queue_state<std::string, std::string>(q, ref, f);
	}

	q.pop_front(); ref.pop_front();
	assert_queue_state<std::string, std::string>(q, ref, f);

	q.push_back("E"); ref.push_back("E");
	assert_queue_state<std::string, std::string>(q, ref, f);
}

TEST_CASE("swag_deque: sum both ends", "[swag_deque]") {
	auto f = [](int a, int b) { return a + b; };
	auto dq = make_swag_deque<int, int>(f);
	std::deque<int> ref;

	dq.push_back(3);  ref.push_back(3);
	assert_deque_state<int,int>(dq, ref, f);

	dq.push_front(2); ref.push_front(2);
	assert_deque_state<int,int>(dq, ref, f);

	dq.push_front(1); ref.push_front(1);
	assert_deque_state<int,int>(dq, ref, f);

	dq.push_back(4);  ref.push_back(4);
	assert_deque_state<int,int>(dq, ref, f);

	dq.pop_front(); ref.pop_front();
	assert_deque_state<int,int>(dq, ref, f);

	dq.pop_back(); ref.pop_back();
	assert_deque_state<int,int>(dq, ref, f);
}

TEST_CASE("swag_deque: preserves order for string concat", "[swag_deque]") {
	auto f = [](const std::string& a, const std::string& b) { return a + b; };
	auto dq = make_swag_deque<std::string, std::string>(f);
	std::deque<std::string> ref;

	dq.push_back("B");  ref.push_back("B");
	assert_deque_state<std::string, std::string>(dq, ref, f);

	dq.push_front("A"); ref.push_front("A");
	assert_deque_state<std::string, std::string>(dq, ref, f);

	dq.push_back("C");  ref.push_back("C");
	assert_deque_state<std::string, std::string>(dq, ref, f);

	dq.push_front("0"); ref.push_front("0");
	assert_deque_state<std::string, std::string>(dq, ref, f);

	dq.pop_back(); ref.pop_back();
	assert_deque_state<std::string, std::string>(dq, ref, f);

	dq.pop_front(); ref.pop_front();
	assert_deque_state<std::string, std::string>(dq, ref, f);
}

TEST_CASE("swag_deque: random stress vs naive (min)", "[swag_deque][stress]") {
	auto f = [](int a, int b) { return std::min(a, b); };
	auto dq = make_swag_deque<int, int>(f);
	std::deque<int> ref;

	std::mt19937 rng(123456);
	std::uniform_int_distribution<int> valdist(-50, 50);
	std::uniform_int_distribution<int> opdist(0, 3); // 0 push_front, 1 push_back, 2 pop_front, 3 pop_back

	for (int step = 0; step < 2000; ++step) {
		int op = opdist(rng);
		if (op == 0) {
			int v = valdist(rng);
			dq.push_front(v);
			ref.push_front(v);
		} else if (op == 1) {
			int v = valdist(rng);
			dq.push_back(v);
			ref.push_back(v);
		} else if (op == 2) {
			if (!ref.empty()) { dq.pop_front(); ref.pop_front(); }
		} else {
			if (!ref.empty()) { dq.pop_back(); ref.pop_back(); }
		}
		assert_deque_state<int,int>(dq, ref, f);
	}
}
