/**
 * SEG_TREE
 *
 * Description: The fundamental ideas are that we decompose any range into O(log n)
 *  nodes of the tree, and you can use pass a callback to many iterator methods to
 *  iterate over that subset or their parents.
 * Copied code from: https://codeforces.com/contest/1916/submission/239673725
 *
 * Usage:
 *  For point_t p (node number on tree):
 *   * Get children: p.c(0), p.c(1)
 *   * For all nodes which range's node contain range's p and iterate up to root: for_ancestor_up
 *   * Same but iterate from root to leaf: for_ancestor_down
 *
 *  For range_t r (multiple nodes on tree that represent the range in array):
 *   * Iterate decomposed range in each nodes and from left to right: for_range
 *   * Same but end function when condition is true: for_range_with_condition
 *   * Same but from right to left: for_reverse_range
 *   * And also end with true condition: for_reverse_range_with_condition
 *   * Iterate over the range from outside - in: for_each or for_each_with_side
 *
 *   * Iterate up all ancestors of all decomposed nodes: for_ancestor_up
 *   * Same but iterate down: for_ancestor_down
 *
 *  Based on node distribution, we have two type of tree:
 *   * in_order_tree: normal tree but compressed to fit into 2 * N nodes
 *   * circular_tree: weird tree but really fast (based on Codeforce blog)
 *
 * More deatail usages:
 *  Make a function update for internal nodes in tree: auto update = ...
 *  Make a function downdate for internal nodes in tree (lazy propagation): auto downdate = ...
 *
 * With point_update/range_query:
 *  * Update at index i:
 *     seg_tree::point_t pt = seg.point(i);
 *     ** Change at leaf node pt (i.e segs[pt] = val) **
 *     pt.for_ancestor_up(update)
 *
 *  * Get query in range [L, R]:
 *     seg_tree::range_t rng = seg.range(L, R + 1);
 *     ** Iterate nodes in range (i.e rng.for_range(...)) **
 *
 * With range_update/point_query:
 *  * Update range [L, R]:
 *     seg_tree::range_t r = seg.range(L, R + 1);
 *     ** Update at some nodes lazily (i.e rng.for_range(...)) **
 *
 *  * Get query at index i:
 *     seg_tree::point_t pt = seg.point(i);
 *     ** Go up to root and get info from internal nodes **
 *
 * With range_update/range_query (lazy propagation):
 *  * Update range [L, R]:
 *     seg_tree::range_t rng = seg.range(L, R + 1);
 *     rng.for_ancestor_down(downdate);
 *     ** Iterate all internal nodes inside range and update it **
 *     rng.for_ancestor_up(update);
 *
 *  * Get query in range [L, R]:
 *     seg_tree::range_t rng = seg.range(L, R + 1);
 *     rng.for_ancestor_down(downdate);
 *     ** Iterate all internal nodes inside range and get infor it **
 *
 * For binary_search on the seg_tree in O(log_2(N)):
 *  seg_tree::range_t rng = seg.range(L, R + 1);
 *  seg_tree::point_t pt = rng.for_range_with_condition(...);
 *  downdate(pt); // very important, if not this node can node downdate anything
 *  int idx = seg.for_descendant(pt, [&](seg_tree::point_t p) {
 *      downdate(p);
 *      check_point(p);
 *  });
**/
namespace seg_tree {

inline int floor_log_2(const int &a) {
	return a ? 31 - __builtin_clz(a) : - 1;
}
inline int ceil_log_2(const int &a) {
	return a ? floor_log_2((a << 1) - 1) : - 1;
}
inline int next_pow_2(const int &a) {
	return 1 << ceil_log_2(a);
}

struct point_t {
	int a;
	point_t() : a(0) {}
	point_t(const int& a_) : a(a_) {}

	explicit operator bool () { return bool(a); }
	operator int() const { return a; }

	inline point_t c(bool z) const { return point_t((a << 1) | z); }
	inline point_t operator [] (bool z) const { return c(z); }
	inline point_t p() const { return point_t(a >> 1); }

	template <typename U>
	friend U& operator >> (U& o, point_t& p) { return o >> int(p); }
	template <typename U>
	friend U& operator << (U& o, const point_t& p) { return o << int(p); }
	friend void __print(const point_t& p) { std::cerr << '(' << int(p) << ')'; }

	template <typename F> void for_each(const F& f) const {
		for (int v = a; v > 0; v >>= 1) {
			f(point_t(v));
		}
	}

	template <typename F> void for_ancestor_down(const F& f) const {
		for (int L = floor_log_2(a); L > 0; L--) {
			f(point_t(a >> L));
		}
	}

	template <typename F> void for_ancestor_up(const F& f) const {
		for (int v = a >> 1; v > 0; v >>= 1) {
			f(point_t(v));
		}
	}

	point_t& operator ++ () { ++ a; return *this; }
	point_t operator ++ (int) { return point_t(a ++); }
	point_t& operator -- () { -- a; return *this; }
	point_t operator -- (int) { return point_t(a --); }
};

struct range_t {
	int a, b;
	range_t() : a(1), b(1) {}
	range_t(const int& a_, const int& b_) : a(a_), b(b_) {}
	range_t(const std::array<int, 2>& r) : range_t(r[0], r[1]) {}

	explicit operator std::array<int, 2>() const { return {a, b}; }
	inline const int& operator[] (bool z) const { return z ? b : a; }

	template <typename U>
	friend U& operator << (U& o, range_t& r) { return o >> r.a >> r.b; }
	template <typename U>
	friend U& operator << (U& o, const range_t& r) { return o << "[" << r.a << ".." << r.b << ")"; }
	friend void __print(const range_t& r) { std::cerr << "[" << r.a << " ... " << r.b << ")"; }

	// Iterate over the range from outside-in.
	template <typename F> void for_each(const F& f) const {
		for (int x = a, y = b; x < y; x >>= 1, y >>= 1) {
			if (x & 1) f(point_t(x ++));
			if (y & 1) f(point_t(-- y));
		}
	}
	// NOTE: Calls f(point_t a, bool is_right)
	template <typename F> void for_each_with_side(const F& f) const {
		for (int x = a, y = b; x < y; x >>= 1, y >>= 1) {
			if (x & 1) f(point_t(x ++), false);
			if (y & 1) f(point_t(-- y), true);
		}
	}

	// Iterate over the range from left to right.
	//    Calls f(point_t)
	template <typename F> void for_range(const F& f) const {
		int anc_depth = floor_log_2((a - 1) ^ b);
		int anc_msk = (1 << anc_depth) - 1;
		for (int v = (- a) & anc_msk; v; v &= v - 1) {
			int i = __builtin_ctz(v);
			f(point_t(((a - 1) >> i) + 1));
		}
		for (int v = b & anc_msk; v; ) {
			int i = floor_log_2(v);
			f(point_t((b >> i) - 1));
			v ^= (1 << i);
		}
	}

	// Iterate over the range from right to left.
	//    Calls f(point_t)
	template <typename F> point_t for_reverse_range_with_condition(const F& f) const {
		int anc_depth = floor_log_2((a - 1) ^ b);
		int anc_msk = (1 << anc_depth) - 1;
		for (int v = b & anc_msk; v; v &= v - 1) {
			int i = __builtin_ctz(v);
			if (f(point_t((b >> i) - 1))) return point_t((b >> i) - 1);
		}
		for (int v = (- a) & anc_msk; v; ) {
			int i = floor_log_2(v);
			if (f(point_t(((a - 1) >> i) + 1))) return point_t(((a - 1) >> i) + 1);
			v ^= (1 << i);
		}
		return point_t(- 1);
	}

	// Iterate over the range from left to right.
	//    Calls f(point_t)
	template <typename F> point_t for_range_with_condition(const F& f) const {
		int anc_depth = floor_log_2((a - 1) ^ b);
		int anc_msk = (1 << anc_depth) - 1;
		for (int v = (- a) & anc_msk; v; v &= v - 1) {
			int i = __builtin_ctz(v);
			if (f(point_t(((a - 1) >> i) + 1))) return point_t(((a - 1) >> i) + 1);
		}
		for (int v = b & anc_msk; v; ) {
			int i = floor_log_2(v);
			if (f(point_t((b >> i) - 1))) return point_t((b >> i) - 1);
			v ^= (1 << i);
		}
		return point_t(- 1);
	}

	// Iterate over the range from right to left.
	//    Calls f(point_t)
	template <typename F> void for_reverse_range(const F& f) const {
		int anc_depth = floor_log_2((a - 1) ^ b);
		int anc_msk = (1 << anc_depth) - 1;
		for (int v = b & anc_msk; v; v &= v - 1) {
			int i = __builtin_ctz(v);
			f(point_t((b >> i) - 1));
		}
		for (int v = (- a) & anc_msk; v; ) {
			int i = floor_log_2(v);
			f(point_t(((a - 1) >> i) + 1));
			v ^= (1 << i);
		}
	}

	template <typename F> void for_ancestor_down(const F& f) const {
		int x = a, y = b;
		if ((x ^ y) > x) { x <<= 1, std::swap(x, y); }
		int dx = __builtin_ctz(x);
		int dy = __builtin_ctz(y);
		int anc_depth = floor_log_2((x - 1) ^ y);
		for (int i = floor_log_2(x); i > dx; -- i) {
			f(point_t(x >> i));
		}
		for (int i = anc_depth; i > dy; -- i) {
			f(point_t(y >> i));
		}
	}

	template <typename F> void for_ancestor_up(const F& f) const {
		int x = a, y = b;
		if ((x ^ y) > x) { x <<= 1, std::swap(x, y); }
		int dx = __builtin_ctz(x);
		int dy = __builtin_ctz(y);
		int anc_depth = floor_log_2((x - 1) ^ y);
		for (int i = dx + 1; i <= anc_depth; ++ i) {
			f(point_t(x >> i));
		}
		for (int v = y >> (dy + 1); v; v >>= 1) {
			f(point_t(v));
		}
	}
};

struct in_order_tree {
	int N, S;
	in_order_tree() : N(0), S(0) {}
	in_order_tree(const int& N_) : N(N_), S(N ? next_pow_2(N) : 0) {}

	inline size_t size() const { return N << 1; }
	inline bool is_leaf(const point_t& pt) const { return int(pt) >= N; }
	inline point_t point(int a) const {
		// With a is a index in array, attempt to find node in tree
		// In seg_tree, leaf nodes sometimes in 1 or 2 layer
		// If in lowest layer = a + S and in higher layer = a + S - N (larger index)
		a += S;
		return point_t(a >= (N << 1) ? a - N : a);
	}
	inline range_t range(int a, int b) const {
		// With range [a, b), we try to find a virtual (or non - virtual) range
		// (lowest nodes) that wrap [a, b)
		if (N == 0) return range_t();
		a += S, b += S;
		if (a >= (N << 1)) a = (a - N) << 1;
		if (b >= (N << 1)) b = (b - N) << 1;
		return range_t(a, b);

	}
	inline range_t range(const std::array<int, 2>& p) const { return range(p[0], p[1]); }
	inline int leaf_index(const point_t& pt) const {
		// pt is the node of the leaf in tree
		// We try to get the index of that leaf in the real array
		int a = int(pt);
		return (a < S ? a + N : a) - S;
	}
	inline std::array<int, 2> node_bound(const point_t& pt) const {
		// With node in seg_tree, we find the [L, R) that is the range in real array
		int a = int(pt);
		int l = __builtin_clz(a) - __builtin_clz((N << 1) - 1);
		int x = a << l, y = (a + 1) << l;
		return {(x >= (N << 1) ? (x >> 1) + N : x) - S, (y >= (N << 1) ? (y >> 1) + N : y) - S};
	}
	inline int node_split(const point_t& pt) const {
		// With node pt, 2 child is c[0] and c[1] with [L, M) and [M, R)
		// We try to find the value M (middle point)
		int a = int(pt);
		int l = __builtin_clz((a << 1) + 1) - __builtin_clz((N << 1) - 1);
		int x = ((a << 1) + 1) << l;
		return (x >= (N << 1) ? (x >> 1) + N : x) - S;
	}
	inline int node_size(const point_t& pt) const {
		auto bounds = node_bound(pt);
		return bounds[1] - bounds[0];
	}

	// NOTE: point_t pt wasn't downdate, must downdate before use this function
	// For each point_t, we must downdate before or after check_point
	// We don't need to update anything after find the leaf index
	template <typename F>
	inline int for_descendant(const point_t& pt, const F& check_point) const {
		int a = int(pt);
		while (a < N) {
			if (check_point(point_t((a << 1)))) a <<= 1;
			else a = (a << 1) + 1;
		}
		return leaf_index(point_t(a));
	}
	template <typename F>
	inline int for_reverse_descendant(const point_t& pt, const F& check_point) const {
		int a = int(pt);
		while (a < N) {
			if (check_point(point_t((a << 1) + 1))) a = (a << 1) + 1;
			else a <<= 1;
		}
		return leaf_index(point_t(a));
	}
};

// Main idea: https://codeforces.com/blog/entry/18051
struct circular_tree {
	int N;
	circular_tree() : N(0) {}
	circular_tree(const int& N_) : N(N_) {}

	inline size_t size() const { return N << 1; }
	inline bool is_leaf(const point_t& pt) const { return int(pt) >= N; }
	inline point_t point(int a) const {
		// With a is a index in array, attempt to find node in tree
		return point_t(N + a);
	}
	inline range_t range(int a, int b) const {
		if (N == 0) return range_t();
		return range_t(N + a, N + b);

	}
	inline range_t range(const std::array<int, 2>& p) const { return range(p[0], p[1]); }
	inline int leaf_index(const point_t& pt) const {
		int a = int(pt);
		return a - N;
	}
	// NOTE: Returns {x,y} so that 0 <= x < N and 1 <= y <= N
	// If the point is non-wrapping (all leaf are fix into range), then 0 <= x < y <= N
	inline std::array<int, 2> node_bound(const point_t& pt) const {
		// With node in seg_tree, we find the [L, R) that is the range in real array
		int a = int(pt);
		int l = __builtin_clz(a) - __builtin_clz((N << 1) - 1);
		int x = a << l, y = (a + 1) << l;
		return {(x >= (N << 1) ? (x >> 1) : x) - N, (y > (N << 1) ? (y >> 1) : y) - N};
	}
	// NOTE: Returns the split point of the node, such that 1 <= s <= N.
	inline int node_split(const point_t& pt) const {
		// With node pt, 2 child is c[0] and c[1] with [L, M) and [M, R)
		// We try to find the value M (middle point)
		return node_bound(pt.c(0))[1];
	}
	inline int node_size(const point_t& pt) const {
		auto bounds = node_bound(pt);
		int r = bounds[1] - bounds[0];
		return r > 0 ? r : r + N;
	}
};

} // namespace seg_tree
