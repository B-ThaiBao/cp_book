/**
 * SPARSE_TABLE !!!
 * 
 * Usage:
 *   - Create the table size N: sparse_table<T> st(N)
 *   - Init the value for index i: st(0, i) = ???
 *   - Build the table: st.build(...);
 *   - Get range [L, R): st.range(L, R);
 *   - Iterate range from [L, R): st.for_range(L, R, ...);
 *   - Iterate range [L, R) in the reverse order: st.for_reverse_range(L, R)
 * 
 * Example: https://codeforces.com/contest/380/submission/247837043
**/
template <typename T> struct sparse_table {
	struct point_t {
		int dist = 0;
		int idx = 0;

		point_t() {}
		point_t(const int& x, const int& y) : dist(x), idx(y) {}

		constexpr point_t c(const bool& z) const {
			return point_t{dist - 1, idx + (z ? (1 << (dist - 1)) : 0)};
		}
		constexpr point_t operator [] (bool z) const { return c(z); }

		constexpr int range(const bool& z) const { return z ? idx : (idx + (1 << dist)); }
		constexpr std::array<int, 2> range() const {
			return std::array<int, 2>{idx, idx + (1 << dist)};
		}

		template <typename U>
		friend U& operator << (U& o, const point_t& p) {
			return o << '[' << p.idx << " .. " << (p.idx + (1 << p.dist)) << ')';
		}
		friend void __print(const point_t& p) { std::cerr << p; }
	};

	int N, log;
	std::vector<std::vector<T>> data;

	static constexpr int log_2(const int& F) {
		return 31 - __builtin_clz(F);
	}

	constexpr size_t size() const { return N; }
	sparse_table() {}
	sparse_table(const int& N_) : N(N_) {
		log = log_2(N);
		data.resize(log + 1);
		data[0].resize(N);
	}

	template <typename F> inline void build(const F& f) {
		for (int i = 1; i <= log; ++ i) {
			data[i].resize(N - (1 << i) + 1);
			for (int u = 0; u + (1 << i) <= N; ++ u) {
				f(point_t(i, u));
			}
		}
	}

	// NOTE: [L, R)
	constexpr std::array<point_t, 2> range(const int& L, const int& R) const {
		int K = log_2(R - L);
		return std::array<point_t, 2>{point_t(K, L), point_t(K, R - (1 << K))};
	}
	template <typename F>
	constexpr void for_range(const int& L, const int& R, const F& f) const {
		int K = R - L;
		for (int j = 0, l = L; (1 << j) <= K; ++ j) {
			if ((K >> j) & 1) {
				f(point_t(j, l));
				l += (1 << j);
			}
		}
	}
	template <typename F>
	constexpr void for_reverse_range(const int& L, const int& R, const F& f) const {
		int K = R - L;
		for (int j = log_2(K), r = R; j >= 0; -- j) {
			if ((K >> j) & 1) {
				f(point_t(j, r - (1 << j)));
				r -= (1 << j);
			}
		}
	}

	inline T& operator [] (const point_t& p) { return data[p.dist][p.idx]; }
	inline const T& operator [] (const point_t& p) const { return data[p.dist][p.idx]; }
	inline T& operator () (const point_t& p) { return data[p.dist][p.idx]; }
	inline const T& operator () (const point_t& p) const { return data[p.dist][p.idx]; }
		
	inline T& operator [] (const std::array<int, 2>& p) { return data[p[0]][p[1]]; }
	inline const T& operator [] (const std::array<int, 2>& p) const { return data[p[0]][p[1]]; }
	inline T& operator () (const int& x, const int& y) { return data[x][y]; }
	inline const T& operator () (const int& x, const int& y) const { return data[x][y]; }
};

/**
 * DISJOINT_SPARSE_TABLE !!!
 * 
 * Disjoint sparse table with type T and N elements provides seg_tree
 * with log_2(N) height, at each height H, exist the mid point M, that
 * can be considered as a prefix and suffix array at point M. In other
 * words, each query can be non-associative based on suffix and prefix
 * properties at the middle of range [L, R]
 * 
 * NOTE: This build the table based on the front from [L, M] and the
 * back from [M + 1, R]
 * 
 * Usage:
 *   Create the table and data structure in O(N * log_2(N)):
 *     disjoint_sparse_table<T> dst(A, front, back): front is add_front and
 *             back is add_back
 *     disjoint_sparse_table<T> dst(A, f): f is combine(prefix, new_node).
 *     disjoint_sparse_table<T> dst; dst.build(A, f): f is no change.
 * 
 *   Iterate thourgh the range [L, R] only in O(1):
 *     If L != R:
 *       int d = dst.depth(L, R);
 *       query = merge(dst(d, L), dst(d, R));
 *     Else:
 *       query = dst(0, L)
 * 
 * Submission: https://codeforces.com/contest/380/submission/247750050
**/
template <typename T> struct disjoint_sparse_table {
	int N, log;
	std::vector<T> data;

	disjoint_sparse_table() {}
	disjoint_sparse_table(const int& n) : N(n) {
		log = (n > 1) ? 33 - __builtin_clz(n - 1) : 1;
		data.resize(log * (1 << (log - 1)));
	}
	template <typename Q>
	disjoint_sparse_table(const int& n, const Q& val) : N(n) {
		log = (n > 1) ? 33 - __builtin_clz(n - 1) : 1;
		data.assign(log * (1 << (log - 1)), val);
	}

	template <typename F, typename B>
	inline void build(const F& build_front, const B& build_back) {
		for (int h = 1, range, half; h < log; ++ h) {
			range = (1 << h), half = range >> 1; // middle point
			for (int i = half; i < N; i += range){
				(*this)(h, i - 1) = (*this)(0, i - 1);
				for (int j = i - 2; j >= i - half; -- j){
					// From j + 1 --> build j
					build_front(h, j);
				}
				(*this)(h, i) = (*this)(0, i);
				for (int j = i + 1; j < i + half; ++ j){
					// From j - 1 --> build j
					build_back(h, j);
				}
			}
		}
	}
	template <typename F>
	inline void build(const F& f) {
		for (int h = 1, range, half; h < log; ++ h) {
			range = (1 << h), half = range >> 1; // middle point
			for (int i = half; i < N; i += range){
				(*this)(h, i - 1) = (*this)(0, i - 1);
				for (int j = i - 2; j >= i - half; -- j){
					// From j + 1 --> build j
					f(h, j);
				}
				(*this)(h, i) = (*this)(0, i);
				for (int j = i + 1; j < i + half; ++ j){
					// From j - 1 --> build j
					f(h, j);
				}
			}
		}
	}

	inline T& operator () (const int& d, const int& x) { return data[d * N + x]; }
	inline const T& operator () (const int& d, const int& x) const { return data[d * N + x]; }
	inline T& operator [] (const std::array<int, 2>& r) { return data[r[0] * N + r[1]]; }
	inline const T& operator [] (const std::array<int, 2>& r) const { return data[r[0] * N + r[1]]; }

	// NOTE: lo != hi is very important and the range is [lo, hi]
	static constexpr int depth(const int& lo, const int& hi) {
		return 32 - __builtin_clz(lo ^ hi);
	}
};
