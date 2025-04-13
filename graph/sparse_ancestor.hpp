/**
 * SPARSE ANCESTOR !!!
 *
 * Usage:
 *   * combine(i, u): merge node {i - 1, u} and {i - 1, sa[i - 1][u]} by nested loops (lg and N)
 *   * Find kth ancestor and path aggregate: sa.find_ancestor(u, k)
 *   * Find lca and path aggregate: sa.find_lca(u, v)
**/
struct sparse_ancestor : public std::vector<std::vector<int>> {
	using std::vector<std::vector<int>>::size;

	static constexpr int log_2(const int& N) { return 31 - __builtin_clz(N); }
	inline void build_sparse() {
		for (int j = 1; j < int(this->size()); ++j) {
			for (int i = 0; i < int(this->at(j).size()); ++i) {
				if (this->at(j - 1).at(i) < 0) {
					this->at(j).at(i) = -1;
				} else {
					this->at(j).at(i) = this->at(j - 1).at(this->at(j - 1).at(i));
				}
			}
		}
	}
	inline size_t size(const bool& z) const {
		return z ? this->at(0).size() : this->size();
	}

	sparse_ancestor() {}
	explicit sparse_ancestor(std::vector<int>&& par_) {
		int N = int(par_.size());
		int lg = log_2(N);
		this->resize(lg + 1, std::vector<int>(N, -1));
		this->at(0) = std::forward<std::vector<int>>(par_);
		this->build_sparse();
	}
	explicit sparse_ancestor(const std::vector<int>& par_) {
		int N = int(par_.size());
		int lg = log_2(N);
		this->resize(lg + 1, std::vector<int>(N, -1));
		this->at(0) = par_;
		this->build_sparse();
	}

	template <typename F> inline int find_ancestor(int u, int depth, const F& f) const {
		int lg = log_2(depth);
		for (int i = 0; i <= lg; ++i) {
			if ((depth >> i) & 1) {
				f(i, u);
				u = this->at(i).at(u);
				if (u == -1) return -1;
			}
		}
		return u;
	}
	inline int find_ancestor(int u, int depth) const {
		int lg = log_2(depth);
		for (int i = 0; i <= lg; ++i) {
			if ((depth >> i) & 1) {
				u = this->at(i).at(u);
				if (u == -1) return -1;
			}
		}
		return u;
	}

	template <typename Vector, typename F>
	inline int find_lca(int u, int v, const Vector& depth, const F& f) const {
		if (depth[u] < depth[v]) std::swap(u, v);
		int d = depth[u] - depth[v];
		int lg = log_2(d);
		for (int i = 0; i <= lg; ++i) {
			if ((d >> i) & 1) {
				f(i, u);
				u = this->at(i).at(u);
			}
		}
		if (u == v) return u;
		for (int i = log_2(depth[v]); i >= 0; --i) {
			if (this->at(i).at(u) != this->at(i).at(v)) {
				f(i, u);
				f(i, v);
				u = this->at(i).at(u);
				v = this->at(i).at(v);
			}
		}
		f(0, u);
		f(0, v);
		return this->at(0).at(u);
	}
	template <typename Vector>
	inline int find_lca(int u, int v, const Vector& depth) const {
		if (depth[u] < depth[v]) std::swap(u, v);
		int d = depth[u] - depth[v];
		int lg = log_2(d);
		for (int i = 0; i <= lg; ++i) {
			if ((d >> i) & 1) {
				u = this->at(i).at(u);
			}
		}
		if (u == v) return u;
		for (int i = log_2(depth[v]); i >= 0; --i) {
			if (this->at(i).at(u) != this->at(i).at(v)) {
				u = this->at(i).at(u);
				v = this->at(i).at(v);
			}
		}
		return this->at(0).at(u);
	}
};
