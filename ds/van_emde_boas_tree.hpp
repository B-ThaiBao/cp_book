/**
 * Main idea: https://www.youtube.com/watch?v=hmReJCupbNU
 * 
 * Code from: https://judge.yosupo.jp/submission/76516
 * 
 * Details:
 *   * van_emde_boas_tree<B> try to maintain a set of int in
 *   * range [0, 2 ^ B) by using array indicies (ensure that
 *   * 2 ^ B fixed into array-index)
 * 
 * Supports:
 *   * find_next(i): returns minimum j >= i in set, or 2 ^ B if not exist
 *   * find_prev(i): returns maximum j <= i in set, or - 1 if not exist
 *   * insert(i), erase(i): insert / erase i into the set
 *   * empty(): returns TRUE if the set is empty
 *   * clear(): empties the set
 * 
 * NOTE: All operations except empty and clear are O(log B) = O(log log 2^B)
**/
template <const int B, typename ENABLE = void> struct van_emde_boas_tree {
	static constexpr int K = B >> 1, R = (B + 1) >> 1, M = (1 << B);
	static constexpr int S = 1 << K, MASK = (1 << R) - 1;

	std::array<van_emde_boas_tree<R>, S> c;
	van_emde_boas_tree<K> act;
	int mi = M, ma = - 1;

	inline bool empty() const { return ma < mi; }
	inline int max() const { return ma; }
	inline int min() const { return mi; }

	inline void insert(int i) {
		if (i <= mi) {
			if (i == mi) return;
			std::swap(mi, i);
			if (i == M) ma = mi; // we were empty
			if (i >= ma) return; // we had mi == ma
		} else if (i >= ma) {
			if (i == ma) return;
			std::swap(ma, i);
			if (i <= mi) return; // we had mi == ma
		}
		int j = i >> R;
		if (c[j].empty()) act.insert(j);
		c[j].insert(i & MASK);
	}
	inline int find_next(const int& i) const {
		if (i <= mi) return mi;
		if (i > ma) return M;
		int j = i >> R, x = i & MASK;
		int res = c[j].find_next(x);
		if (res <= MASK) return (j << R) + res;
		j = act.find_next(j + 1);
		return (j >= S) ? ma : ((j << R) + c[j].find_next(0));
	}
	inline int find_prev(const int& i) const {
		if (i >= ma) return ma;
		if (i < mi) return - 1;
		int j = i >> R, x = i & MASK;
		int res = c[j].find_prev(x);
		if (res >= 0) return (j << R) + res;
		j = act.find_prev(j - 1);
		return (j < 0) ? mi : ((j << R) + c[j].find_prev(MASK));
	}
	inline void erase(int i) {
		if (i <= mi) {
			if (i < mi) return;
			i = mi = find_next(mi + 1);
			if (i >= ma) {
				if (i > ma) ma = -1;  // we had mi == ma
				return;               // after erase we have mi == ma
			}
		} else if (i >= ma) {
			if (i > ma) return;
			i = ma = find_prev(ma - 1);
			if (i <= mi) return;  // after erase we have mi == ma
		}
		int j = i >> R;
		c[j].erase(i & MASK);
		if (c[j].empty()) act.erase(j);
	}
	inline bool find(const int& i) const {
		if (i < 0 || i >= M) return false;
		if (i == mi || i == ma) return true;
		return c[i >> R].find(i & MASK);
	}

	inline void clear() {
		mi = M, ma = -1;
		act.clear();
		for (int i = 0; i < S; ++ i) c[i].clear();
	}
};

template <const int B> struct van_emde_boas_tree<B, std::enable_if_t<(B <= 6)>> {
	static constexpr int M = 1 << B;
	uint64_t act = 0;

	inline bool empty() const { return act == 0; }
	inline void clear() { act = 0; }
	inline int min() const { return __builtin_ctzll(act); }
	inline int max() const { return __builtin_clzll(act); }
	inline int find_next(const int& i) const {
		if (i == M) return M;
		uint64_t tmp = act >> i;
		return (tmp ? i + __builtin_ctzll(tmp) : M);
	}
	inline int find_prev(const int& i) const {
		if (i == - 1) return - 1;
		uint64_t tmp = act << (63 - i);
		return (tmp ? i - __builtin_clzll(tmp) : -1);
	}
	inline void insert(const int& i) { act |= 1ULL << i; }
	inline void erase(const int& i) { act &= ~(1ULL << i); }
	inline bool find(const int& i) const {
		if (i < 0 || i >= M) return false;
		return (act >> i) & 1;
	}
};