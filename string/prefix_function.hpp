/**
 * NOTE: The suffix match with the prefix at position i is p[i - 1], p[p[i - 1] - 1], ...
 *                (this is on the decreasing order of the length)
**/
struct prefix_function : public std::vector<int> {
	using std::vector<int>::vector;

	template <typename String> inline void build(const String& S) {
		int N = int(S.size());
		(*this).resize(N);
		(*this)[0] = 0;
		for (int i = 1, j = 0; i < N; ++ i) {
			while (j > 0 && S[i] != S[j]){
				j = (*this)[j - 1];
			}
			if (S[i] == S[j]) ++ j;
			(*this)[i] = j;
		}
	}
	template <typename String> prefix_function(const String& S) { (*this).build(S); }
};

// Return 0-indexed positions of occurrences of s in w in O(n + m)
// NOTE: pi is prefix_function or prefix_automaton of s
template <typename String, typename Prefix>
std::vector<int> find_pos(const String& w, const String& s, const Prefix& pi) {
	std::vector<int> res; res.reserve(w.size() + 1 - s.size());
	int N = int(s.size());
	for (int i = 0, k = 0; i < int(w.size()); i ++) {
		while (k > 0 && (k == N || w[i] != s[k])){
			k = pi[k - 1];
		}
		if (w[i] == s[k]) ++ k;
		if (k == N){
			res.emplace_back(i - N + 1); // satisfied position
		}
	}	
	return res;
}

// Return num of occurrences of substring S[0 ... i] in string S: time O(n)
// NOTE: pi is prefix_function or prefix_automaton of s
template <typename Prefix> std::vector<int> count_prefix(const Prefix& pi) {
	int N = pi.size();
	std::vector<int> res(N, 1);
	for (int i = N - 1; i >= 0; -- i){
		if (pi[i] > 0){
			res[pi[i] - 1] += res[i];
		}
	}
	return res;
}

// Return compression of substring S[0 ... i] in string S: time O(1)
// NOTE: pi is prefix_function or prefix_automaton of s
template <typename Prefix>
int compress_prefix(const Prefix& pi, int N = - 1) {
	if (N == - 1) N = int(pi.size()) - 1;
	int K = N + 1 - pi[N];
	return ((N + 1) % K == 0) ? K : N + 1;
}

/**
 * PREFIX AUTOMATON:
 *   Make a huge contribution on DP problems.
 *   Solving DP problems, you want to add new char into string ==> you want to care about the
 *   suffix so KMP state is very powerful.
 * 
 *   Details, when you add character into a string and you know the KMP state of before ==>
 *   you can know the KMP state of this position.
 *   KMP state is a longest suffix of this string that also the prefix of other string.
 * 
 *   For more details, you can look at this solution:
 *     *  https://lqdoj.edu.vn/submission/2811916
 *     *  https://codeforces.com/contest/808/submission/221319609
**/
template <typename T, typename F = std::function<int(const T&)>>
struct prefix_automaton : public std::vector<int> {
	using std::vector<int>::vector;
	using std::vector<int>::operator [];

	int num_state = 0; // How many to change from this state to another state
	std::vector<int> next;
	F map; // Map each charater fix into [0 ... num_state)

	inline int& operator [] (const std::pair<int, T> &x) { return next[x.first * num_state + map(x.second)]; }
	inline const int& operator [] (const std::pair<int, T> &x) const { return next[x.first * num_state + map(x.second)]; }
	inline int& operator () (const int& x) { return (*this)[x]; }
	inline const int& operator () (const int& x) const { return (*this)[x]; }
	inline int& operator () (const int& x, const T& y) { return next[x * num_state + map(y)]; }
	inline const int& operator () (const int& x, const T& y) const { return next[x * num_state + map(y)]; }

	template <typename String>
	inline void build(const String& S, const int& state) {
		num_state = state;
		int N = int(S.size());
		(*this).resize(N);
		(*this)[0] = 0;
		for (int i = 1, j = 0; i < N; ++ i) {
			while (j > 0 && S[i] != S[j]) {
				j = (*this)[j - 1];
			}
			if (S[i] == S[j]) ++ j;
			(*this)[i] = j;
		}

		// Now, we try to build state for each character
		next.resize(num_state * (N + 1));
		for (int i = 0; i <= N; ++ i){
			for (int j = 0; j < num_state; ++ j) {
				auto mp = map(S[i]);
				if ((i > 0 && j != mp) || i == N) {
					next[i * num_state + j] = next[(*this)[i - 1] * num_state + j];
				} else {
					next[i * num_state + j] = i + (j == mp);
				}
			}
		}
	}
	template <typename String, typename L>
	prefix_automaton(const String& S, const int& state, const L& f) : map(f) {
		(*this).build(S, state);
	}
};

// NOTE: This is useful to use lambda (instead of std::function)
template <typename String, typename Lambda>
static inline auto build_prefix_automaton(const String& S, const int& state, const Lambda& f) {
	using char_t = typename String::value_type;
	prefix_automaton<char_t, typename std::decay<Lambda>::type> pa(S, state, f);
	return pa;
}