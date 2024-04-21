// NOTE: This can store more infor in each set (i.e max, min, sum, ...)
// Furthermore, trick from: https://codeforces.com/blog/entry/75369 (dealing with queries)
// You can store the query-index (query-form: a and b) for each node and do operations:
//   - Merge: If two sets have the same query-index (called i) so i-th query can be answer.
//            The remaining query-index can be merge like normal and meet in the next merging
//   - Find_par: unchange
struct disjoint_set {
	int N;
	std::vector<int> data;

	disjoint_set() {}
	disjoint_set(const int& N) : N(N), data(N, - 1) {}

	int size(const int& x) { return - data[find_par(x)]; }
	const int& size() const { return N; }

	void init(const int& N) {
		this -> N = N;
		data.assign(N, - 1);
	}

	int find_par(const int& x) {
		return data[x] < 0 ? x : data[x] = find_par(data[x]);
	}
	template <typename F> int find_par(const int& x, const F& f) {
		if (data[x] < 0) return x;
		int root = find_par(data[x], f);
		// TODO: We already go par, now we try to get from par (like dist, ...)
		f(data[x], x);
		return data[x] = root;
	}

	bool merge_par(const int& a, const int& b) {
		data[a] += data[b];
		data[b] = a; -- N;
		return true;
	}
	template <typename F>
	bool merge_par(const int& a, const int& b, const F& f) {
		f(a, b);
		data[a] += data[b];
		data[b] = a; -- N;
		return true;
	}

	bool merge(int a, int b) {
		a = find_par(a); 
		b = find_par(b);
		if (a == b) return false;
		if (data[a] > data[b]) std::swap(a, b);
		merge_par(a, b);
		return true;
	}
	template <typename F> bool merge(int a, int b, const F& f) {
		a = find_par(a); 
		b = find_par(b);
		if (a == b) return false;
		if (data[a] > data[b]) std::swap(a, b);
		merge_par(a, b, f);
		return true;
	}
	template <typename P, typename F>
	bool merge(int a, int b, const P& get, const F& f) {
		a = find_par(a ,get);
		b = find_par(b, get);
		if (a == b) return false;
		if (data[a] > data[b]) std::swap(a, b);
		merge_par(a, b, f);
		return true;
	}
};

struct bipartite_graph : public disjoint_set {
	/**
	 * TODO: Try to keep the bipartite of components
	 * by using disjoint set and based on the property:
	 *  The bipartite graph doesn't allow odd cycle !!!
	**/
	using disjoint_set::data;
	std::vector<bool> is;
	std::vector<bool> p;

	bipartite_graph() {}
	bipartite_graph(const int& N) : disjoint_set(N) {
		is.assign(N, true);
		p.assign(N, true);
	}

	bool is_same(const int& x, const int& y) {
		auto from_par = [&](const int& par, const int& c) {
			p[c] = p[c] ^ p[par];
		};
		return (*this).find_par(x, from_par) == (*this).find_par(y, from_par);
	}

	void add_edge(int a, int b) {
		auto from_par = [&](const int& par, const int& c) {
			p[c] = p[c] ^ p[par];
		};
		a = (*this).find_par(a, from_par);
		b = (*this).find_par(b, from_par);
		if (a == b) {
			if (p[a] == p[b]) is[a] = false; // Same component but odd cycle
		}
		else {
			if (data[a] > data[b]) std::swap(a, b);
			p[b] = p[a] ^ p[b] ^ 1;
			is[a] = is[a] && is[b];
			data[a] += data[b];
			data[b] = a;
		}
	}

	bool is_bipartite(const int& v) {
		auto from_par = [&](const int& par, const int& c) {
			p[c] = p[c] ^ p[par];
		};
		return is[(*this).find_par(v, from_par)];
	}
};

// NOTE: This supports the pointer point the result at the moment
// When the leader cannot be the result anymore, we merge them with the others
struct disjoint_set_ancestor {
	int N;
	std::vector<int> data;

	disjoint_set_ancestor() {}
	disjoint_set_ancestor(const int& N) : N(N), data(N) {
		std::iota(data.begin(), data.end(), 0);
	}
	const int& size() const { return N; }

	void init(const int& N) {
		this -> N = N;
		data.resize(N);
		std::iota(data.begin(), data.end(), 0);
	}

	int find_par(const int& x) {
		return data[x] == x ? x : data[x] = find_par(data[x]);
	}
	template <typename F> int find_par(const int& x, const F& f) {
		if (data[x] == x) return x;
		int root = find_par(data[x], f);
		// TODO: We already go par, now we try to get from par (like dist, ...)
		f(data[x], x);
		return data[x] = root;
	}

	bool merge_par(const int& a, const int& b) {
		data[b] = a; -- N;
		return true;
	}
	template <typename F>
	bool merge_par(const int& a, const int& b, const F& f) {
		f(a, b);
		data[b] = a; -- N;
		return true;
	}

	bool merge(int a, int b) {
		a = find_par(a); 
		b = find_par(b);
		if (a == b) return false;
		merge_par(a, b);
		return true;
	}
	template <typename F> bool merge(int a, int b, const F& f) {
		a = find_par(a); 
		b = find_par(b);
		if (a == b) return false;
		merge_par(a, b, f);
		return true;
	}
	template <typename P, typename F>
	bool merge(int a, int b, const P& get, const F& f) {
		a = find_par(a ,get);
		b = find_par(b, get);
		if (a == b) return false;
		merge_par(a, b, f);
		return true;
	}
};