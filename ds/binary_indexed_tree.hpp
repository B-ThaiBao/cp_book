/**
 * BINARY_INDEXED_TREE (BIT) !!!
 * 
 * Usage:
 * This tree provides two main functions:
 * 		- prefix(i) : all range from i -> 0 (using by for (auto x : p.prefix(i)))
 * 		- suffix(i) : all range from i -> N (using by for (auto x : p.suffix(i)))
 * 
 * Binary_indexed_tree store information is quite complex, but 3 things important
 * to know is:
 * 
 * 		- All range from query can not intesect with each other when using for
 * 		==> using for commutative operator with easy
 * 
 * 		- At node i, the tree can store range at least (prefix) or at most (suffix)
 * 
 * 		- Can be considered as a virtual array
 * 	
 * Furthermore, this template provides some optimized function like:
 * 
 * 		- build_prefix(): build prefix tree only in O(N)
 * 		- build_suffix(): same but for suffix tree
 * 		- lower_bound(): return the first index that sum[0 ... pos] >= x only in O(logN)
 * 		- upper_bound(): return the first index that sum[0 ... pos] > x only in O(logN)
 * 
 * How to use for update and get query:
 * 		* For prefix_query + point_update: for_prefix() for query and for_suffix() for update
 * 		* For suffix_query + point_update: for_suffix() for query and for_prefix() for update
 * 		* For point_query + range_update: Apply difference array and prefix_query + point_update (no change)
 * 		* For range_query + range_update: Using range_update_sum_query below (using two tree method)
**/
template <typename T> struct binary_indexed_tree : public std::vector<T> {
	using std::vector<T>::vector;

	template <typename Vec, typename F>
	void build_prefix(const Vec& A, F&& set_num) {
		if ((*this).size() != A.size()) {;
			(*this).resize(A.size());
		}
		int N = (int) A.size();
		for (int i = 0; i < N; ++ i) {
			set_num((*this)[i], A[i]);
			int par = i | (i + 1);
			if (par < N) set_num((*this)[par], (*this)[i]);
		}
	}
	template <typename Vec, typename F>
	void build_suffix(const Vec& A, F&& set_num) {
		if ((*this).size() != A.size()) {;
			(*this).resize(A.size());
		}
		int N = (int) A.size();
		for (int i = N - 1; i >= 0; -- i) {
			set_num((*this)[i], A[i]);
			int par = (i & (i + 1)) - 1;
			if (par >= 0) set_num((*this)[par], (*this)[i]);
		}
	}

	template <typename B, typename E = B> struct iterator_range {
		B beg;
		E en;

		iterator_range() : beg(), en(){}
		iterator_range(const B& b, const E& e) : beg(b), en(e) {}
		iterator_range(B&& b, E&& e) : beg(b), en(e) {}
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
			idx = (idx & (idx + 1)) - 1; return *this;
		}
		T& operator * () const { return data[idx]; }
	};
	using prefix_range = iterator_range<prefix_iterator>;

	prefix_range prefix(int x) {
		int N = int(this -> size());
		if (x >= N) x = N - 1; 
		return prefix_range{prefix_iterator{(*this).data(), x}, prefix_iterator{nullptr, - 1}};
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
			idx = (idx & (idx + 1)) - 1; return *this;
		}
		const T& operator * () const{ return data[idx]; }
	};
	using const_prefix_range = iterator_range<const_prefix_iterator>;

	const_prefix_range prefix(int x) const {
		int N = int(this -> size());
		if (x >= N) x = N - 1; 
		return const_prefix_range{const_prefix_iterator{(*this).data(), x}, const_prefix_iterator{nullptr, - 1}};
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
			idx |= (idx + 1); return *this;
		}
		T& operator * () const{ return data[idx]; }
	};
	using suffix_range = iterator_range<suffix_iterator>;

	suffix_range suffix(int x) {
		if (x < 0) x = 0;
		return suffix_range{suffix_iterator{(*this).data(), x}, suffix_iterator{nullptr, int((*this).size())}};
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
			idx |= (idx + 1); return *this;
		}
		const T& operator * () const { return data[idx]; }
	};
	using const_suffix_range = iterator_range<const_suffix_iterator>;

	const_suffix_range suffix(int x) const {
		if (x < 0) x = 0;
		return const_suffix_range{const_suffix_iterator{(*this).data(), x}, const_suffix_iterator{nullptr, int((*this).size())}};
	}

	static constexpr int log_2(const int& P) { return 31 - __builtin_clz(P); }

	int lower_bound(const T& V) const {
		T sum{};
		int pos = 0;
		int N = this -> size();
		for (int i = log_2(N); i >= 0; -- i) {
			if (pos + (1 << i) <= N && sum + (*this)[pos + (1 << i) - 1] < V) {
				sum = sum + (*this)[pos + (1 << i) - 1];
				pos += (1 << i);
			}
		}
		return pos;
	}
	int upper_bound(const T& V) const {
		T sum{};
		int pos = 0;
		int N = this -> size();
		for (int i = log_2(N); i >= 0; -- i) {
			if (pos + (1 << i) <= N && sum + (*this)[pos + (1 << i) - 1] <= V) {
				sum = sum + (*this)[pos + (1 << i) - 1];
				pos += (1 << i);
			}
		}
		return pos;
	}
};

template <typename T> struct range_update_sum_query{
	int N;
	// Merge them in only one tree (save memory and much more faster due to using arrays)
	binary_indexed_tree<std :: array<T, 2>> data;

	range_update_sum_query() : N(N), data(){}
	range_update_sum_query(const int& N) : N(N), data(N, {0, 0}){}
	range_update_sum_query(const int& N, const T& t) : N(N), data(N, {t, t}){}
	template <typename Vec>
	range_update_sum_query(const Vec& A){
		N = (int) A.size(); data.data.assign(N, {0, 0});
		std :: vector<std :: array<T, 2>> D(N);
		for (int i = 0; i < N; i ++){
			D[i][1] = (i == 0) ? A[i] : (A[i] - A[i - 1]);
			D[i][0] = (N - i - 1) * D[i][1];
		}
		data.build_prefix(D, [](std :: array<T, 2>& p, const std :: array<T, 2>& c) -> void{
			p[0] += c[0];
			p[1] += c[1];
		});
	}

	template <typename B>
	void range_update(const int& L, const int& R, const B& add){
		for (auto& x : data.suffix(L)){
			x[0] += (T) add * (N - L - 1);
			x[1] += (T) add;
		}
		for (auto& x : data.suffix(R + 1)){
			x[0] -= (T) add * (N - R - 2);
			x[1] -= (T) add;
		}
	}
	template <typename B>
	void point_update(const int& pt, const B& add){
		this -> range_update(pt, pt, add);
	}

	T point_query(const int& pt) const{
		T sum{};
		for (const auto& x : data.prefix(pt)){
			sum += x[0];
			sum -= (T) (N - pt - 2) * x[1];
		}
		return sum;
	}
	T range_query(const int& L, const int& R) const{
		return point_query(R) - ((L == 0) ? T() : point_query(L - 1));
	}

	// FOR DEBUG
	friend void _print(range_update_sum_query<T>& p){
		#ifdef LOCAL
		std :: vector<T> A(p.N);
		for (int i = 0; i < p.N; i ++){
			A[i] = p.range_query(i, i);
		}
		_print(A);
		#endif
	}
};