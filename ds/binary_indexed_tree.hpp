/**
 * BINARY_INDEXED_TREE (BIT) !!!
 *
 * Usage:
 * A binary indexed tree with N nodes of type T purely supports for ranges:
 *   - For prefix_tree of index i: [x, i) with x < i
 *   - For suffix_tree of index i: [i, x) with x > i
 *
 * For that reason, it provides 2 following functions for virtual array:
 *   - prefix(i):
 *       + For prefix tree: list of prefix ranges which can aggregate for query [0, i)
 *       + For suffix tree: list of prefix ranges which contain index i
 *   - suffix(i):
 *       + For suffix tree: list of suffix ranges which can aggregate for query [i, N)
 *       + For prefix tree: list of suffix ranges which contain index i
 * such that each range is not intersect with each other and the list of ranges at most log_2(N)
 *
 * Furthermore, this template provides some optimized function like:
 *   - build_prefix(): build prefix tree only in O(N)
 *   - build_suffix(): same but for suffix tree
 *   - lower_bound(): return the first index that sum[0 ... pos) >= x only in O(logN)
 *   - upper_bound(): return the first index that sum[0 ... pos) > x only in O(logN)
 *
 * With 0-indexed data, we can use this tree for following operations:
 *   - For prefix_query + point_update:
 *       + prefix_query [0, i): for_prefix(i)
 *       + point_update [i]: for_suffix(i)
 *   - For suffix_query + point_update:
 *       + suffix_query [i, N): for_suffix(i)
 *       + point_update [i]: for_prefix(i + 1) (1-indexed)
 *   - For point_query + range_update: Apply difference array and prefix_query + point_update (no change)
 *   - For range_query + range_update: Using range_update_sum_query below (using two tree method)
**/
template <typename T> struct binary_indexed_tree : public std::vector<T> {
	using std::vector<T>::vector;

	template <typename Vec, typename F>
	void build_prefix(const Vec& A, F&& set_num) {
		if (this->size() != A.size()) {;
			this->resize(A.size());
		}
		int N = int(A.size());
		for (int i = 0; i < N; ++i) {
			set_num(this->at(i), A[i]);
			int par = i | (i + 1);
			if (par < N) set_num(this->at(par), this->at(i));
		}
	}
	template <typename Vec, typename F>
	void build_suffix(const Vec& A, F&& set_num) {
		if (this->size() != A.size()) {;
			this->resize(A.size());
		}
		int N = int(A.size());
		for (int i = N - 1; i >= 0; --i) {
			set_num(this->at(i), A[i]);
			int par = (i & (i + 1)) - 1;
			if (par >= 0) set_num(this->at(par), this->at(i));
		}
	}

	template <typename B, typename E = B> struct iterator_range {
		B beg;
		E en;

		iterator_range() : beg(), en() {}
		iterator_range(const B &b, const E &e) : beg(b), en(e) {}
		iterator_range(B &&b, E &&e) : beg(b), en(e) {}
		B begin() const { return beg; }
		E end() const { return en; }
	};

	struct prefix_iterator {
		T* data;
		int idx;

		prefix_iterator() : data(nullptr), idx(0) {}
		prefix_iterator(T* d, const int& i) : data(d), idx(i) {}

		friend inline bool operator != (const prefix_iterator& i, const prefix_iterator& j) {
			return i.idx > j.idx;
		}
		prefix_iterator& operator ++ () {
			idx &= idx - 1;
			return *this;
		}
		T& operator * () const { return data[idx - 1]; }
	};
	using prefix_range = iterator_range<prefix_iterator>;

	prefix_range prefix(const int& x) {
#ifdef _GLIBCXX_DEBUG
		assert(0 <= x && x <= int(this->size()));
#endif
		return prefix_range{prefix_iterator{this->data(), x}, prefix_iterator{nullptr, 0}};
	}

	struct const_prefix_iterator{
		const T* data;
		int idx;

		const_prefix_iterator() : data(nullptr), idx(0) {}
		const_prefix_iterator(const T* d, const int& i) : data(d), idx(i) {}

		friend inline bool operator != (const const_prefix_iterator& i, const const_prefix_iterator& j) {
			return i.idx > j.idx;
		}
		const_prefix_iterator& operator ++ () {
			idx &= idx - 1;
			return *this;
		}
		const T& operator * () const{ return data[idx - 1]; }
	};
	using const_prefix_range = iterator_range<const_prefix_iterator>;

	const_prefix_range prefix(const int& x) const {
#ifdef _GLIBCXX_DEBUG
		assert(0 <= x && x <= int(this->size()));
#endif
		return const_prefix_range{const_prefix_iterator{this->data(), x}, const_prefix_iterator{nullptr, 0}};
	}

	struct suffix_iterator {
		T* data;
		int idx;

		suffix_iterator() : data(nullptr), idx(0) {}
		suffix_iterator(T* d, const int& i) : data(d), idx(i) {}

		friend inline bool operator != (const suffix_iterator& i, const suffix_iterator& j) {
			return i.idx < j.idx;
		}
		suffix_iterator& operator ++ () {
			idx |= idx + 1;
			return *this;
		}
		T& operator * () const{ return data[idx]; }
	};
	using suffix_range = iterator_range<suffix_iterator>;

	suffix_range suffix(const int& x) {
#ifdef _GLIBCXX_DEBUG
		assert(0 <= x && x <= int(this->size()));
#endif
		return suffix_range{suffix_iterator{this->data(), x}, suffix_iterator{nullptr, int(this->size())}};
	}

	struct const_suffix_iterator {
		const T* data;
		int idx;

		const_suffix_iterator() : data(nullptr), idx(0) {}
		const_suffix_iterator(const T* d, const int& i) : data(d), idx(i) {}

		friend inline bool operator != (const const_suffix_iterator& i, const const_suffix_iterator& j) {
			return i.idx < j.idx;
		}
		const_suffix_iterator& operator ++ () {
			idx |= idx + 1;
			return *this;
		}
		const T& operator * () const { return data[idx]; }
	};
	using const_suffix_range = iterator_range<const_suffix_iterator>;

	const_suffix_range suffix(const int& x) const {
#ifdef _GLIBCXX_DEBUG
		assert(0 <= x && x <= int(this->size()));
#endif
		return const_suffix_range{const_suffix_iterator{this->data(), x}, const_suffix_iterator{nullptr, int(this->size())}};
	}

	static constexpr int log_2(const int& P) { return 31 - __builtin_clz(P); }

	template <typename V, typename F = std::plus<T>>
	inline int lower_bound(const V& val, const F& op = F()) const {
		T sum{};
		int pos = 0;
		int N = int(this->size());
		for (int i = log_2(N); i >= 0; --i) {
			if (pos + (1 << i) <= N) {
				auto cur = op(sum, this->at(pos + (1 << i) - 1));
				if (cur < T(val)) {
					sum = cur;
					pos += (1 << i);
				}
			}
		}
		return pos + 1;
	}
	template <typename V, typename F = std::plus<T>>
	inline int upper_bound(const V& val, const F& op = F()) const {
		T sum{};
		int pos = 0;
		int N = int(this->size());
		for (int i = log_2(N); i >= 0; --i) {
			if (pos + (1 << i) <= N) {
				auto cur = op(sum, this->at(pos + (1 << i) - 1));
				if (cur <= T(val)) {
					sum = cur;
					pos += (1 << i);
				}
			}
		}
		return pos + 1;
	}
};

template <typename T> struct range_update_sum_query : binary_indexed_tree<std::array<T, 2>>{
	range_update_sum_query() {}
	range_update_sum_query(const int& N) {
		this->assign(N, {0, 0});
	}
	template <typename Vector> range_update_sum_query(const Vector& A) {
		int N = int(A.size());
		this->assign(N, {0, 0});
		std::vector<std::array<T, 2>> D(N);
		for (int i = 0; i < N; i++) {
			D[i][1] = (i == 0) ? T(A[i]) : T(A[i] - A[i - 1]);
			D[i][0] = T(N - i - 1) * D[i][1];
		}
		this->build_prefix(D, [](auto &p, const auto &c) -> void {
			p[0] += c[0];
			p[1] += c[1];
		});
	}
	range_update_sum_query(const int& N, const T& t) : range_update_sum_query(std::vector<T>(N, t)) {}

	template <typename B> void range_update(const int& l, const int& r, const B& add) {
#ifdef _GLIBCXX_DEBUG
		assert(0 <= l && l <= r && r <= int(this->size()));
#endif
		int N = int(this->size());
		for (auto& x : this->suffix(l)) {
			x[0] += T(add) * (N - l - 1);
			x[1] += T(add);
		}
		for (auto& x : this->suffix(r)) {
			x[0] -= T(add) * (N - r - 1);
			x[1] -= T(add);
		}
	}
	template <typename B> void point_update(const int& i, const B& add) {
#ifdef _GLIBCXX_DEBUG
		assert(0 <= i && i < int(this->size()));
#endif
		this->range_update(i, i, add);
	}

	T prefix_sum(const int& r) const {
#ifdef _GLIBCXX_DEBUG
		assert(0 <= r && r <= int(this->size()));
#endif
		if (r == 0) return T(0);
		T sum{};
		for (const auto& x : this->prefix(r)) {
			sum += x[0];
			sum -= T(int(this->size()) - r - 1) * x[1];
		}
		return sum;
	}
	T range_query(const int& l, const int& r) const {
#ifdef _GLIBCXX_DEBUG
		assert(0 <= l && l <= r && r <= int(this->size()));
#endif
		return this->prefix_sum(r) - this->prefix_sum(l);
	}
};
