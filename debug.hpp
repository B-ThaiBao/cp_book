#pragma once
#include <bits/stdc++.h>

void __print(short x) { std::cerr << x; }
void __print(int x) { std::cerr << x; }
void __print(long x) { std::cerr << x; }
void __print(long long x) { std::cerr << x; }
void __print(unsigned x) { std::cerr << x; }
void __print(unsigned long x) { std::cerr << x; }
void __print(unsigned long long x) { std::cerr << x; }
void __print(float x) { std::cerr << x; }
void __print(double x) { std::cerr << x; }
void __print(long double x) { std::cerr << x; }
void __print(char x) { std::cerr << '\'' << x << '\''; }
void __print(const char *x) { std::cerr << '\"' << x << '\"'; }
void __print(const std::string &x) { std::cerr << '\"' << x << '\"'; }
void __print(bool x) { std::cerr << (x ? "true" : "false"); }
template <typename A> void __print(const A &x);
template <typename A, typename B> void __print(const std::pair<A, B> &p);
template <typename... A> void __print(const std::tuple<A...> &t);
template <typename T> void __print(std::stack<T> s);
template <typename T> void __print(std::queue<T> q);
template <typename T, typename... U> void __print(std::priority_queue<T, U...> q);
template <typename A> void __print(const A &x) {
	bool first = true;
	std::cerr << '{';
	for (const auto &i : x) {
		std::cerr << (first ? "" : ","), __print(i);
		first = false;
	}
	std::cerr << '}';
}
template <typename A, typename B> void __print(const std::pair<A, B> &p) {
	std::cerr << '(';
	__print(p.first);
	std::cerr << ',';
	__print(p.second);
	std::cerr << ')';
}
template <typename T> void __print(std::stack<T> s) {
	std::vector<T> debugVector;
	while (!s.empty()) {
		T t = s.top();
		debugVector.push_back(t);
		s.pop();
	}
	std::reverse(debugVector.begin(), debugVector.end());
	__print(debugVector);
}
template <typename T> void __print(std::queue<T> q) {
	std::vector<T> debugVector;
	while (!q.empty()) {
		T t = q.front();
		debugVector.push_back(t);
		q.pop();
	}
	__print(debugVector);
}
template <typename T, typename... U> void __print(std::priority_queue<T, U...> q) {
	std::vector<T> debugVector;
	while (!q.empty()) {
		T t = q.top();
		debugVector.push_back(t);
		q.pop();
	}
	__print(debugVector);
}
void _print() { std::cerr << "]\n"; }
template <typename Head, typename... Tail>
void _print(const Head &H, const Tail &...T) {
	__print(H);
	if (sizeof...(T)) std::cerr << ", ";
	_print(T...);
}

#ifndef DEBUG
#define debug(...) std::cerr << "Line " << __LINE__ << ": [" << #__VA_ARGS__ << "] = ["; _print(__VA_ARGS__)
#else
#define debug(...) 2004
#endif

template <typename F> void benchmark(const std::string &func_name, const F &f) {
	auto start = std::chrono::high_resolution_clock::now();
	f();
	auto end = std::chrono::high_resolution_clock::now();
	double time_taken = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
	std::cerr << "Benchmark " << func_name << " : " << time_taken << " ns\n" << std::flush;
}
