struct indexed_set {
	static constexpr uint32_t BUCKET = 64;
	static constexpr uint32_t LOG_BUCKET = 6;

	int32_t N, lg;
	std::vector<std::vector<uint64_t>> data;

	static inline int div_bucket(const int& i) { return i >> LOG_BUCKET; }
	static inline int mod_bucket(const int& i) { return (i & (BUCKET - 1)); }
	static inline int low_bit(const uint64_t& x) { return x == 0 ? - 1 : __builtin_ctzll(x); }
	static inline int top_bit(const uint64_t& x) { return x == 0 ? - 1 : 63 - __builtin_clzll(x); }

	inline void build(int M) {
		N = M;
		do {
			int32_t sz = (M + BUCKET - 1) >> LOG_BUCKET;
			data.emplace_back(std::vector<uint64_t>(sz));
			M = sz;
		} while (M > 1);
		lg = int(data.size());
	}
	template <typename Vector>
	inline void build(const int& M, const Vector& vals) {
		build(M);
		for (int i = 0; i < N; ++ i) {
			data[0][ i >> LOG_BUCKET] |= uint64_t(bool(vals(i))) << (i & (BUCKET - 1));
		}
		for (int h = 0; h < lg - 1; ++ h) {
			for (int i = 0; i < int(data[h].size()); ++ i) {
				data[h + 1][i >> LOG_BUCKET] |= uint64_t(bool(data[h][i])) << (i & (BUCKET - 1));
			}
		}
	}

	indexed_set() {}
	indexed_set(const int& M) { build(M); }
	template <typename Vector>
	indexed_set(const int& M, const Vector& vals) { build(M, vals); }

	inline void insert(int i) {
		for (int h = 0; h < lg; ++ h) {
			data[h][i >> LOG_BUCKET] |= uint64_t(1) << (i & (BUCKET - 1)), i >>= LOG_BUCKET;
		}
	}
	inline void erase(int i) {
		uint64_t x = 0;
		for (int h = 0; h < lg; ++ h) {
			data[h][i >> LOG_BUCKET] &= ~(uint64_t(1) << (i & (BUCKET - 1)));
			data[h][i >> LOG_BUCKET] |= x << (i & (BUCKET - 1));
			x = bool(data[h][i >> LOG_BUCKET]);
			i >>= LOG_BUCKET;
		}
	}
	inline bool find(const int& i) const {
		return (data[0][i >> LOG_BUCKET] >> (i & (BUCKET - 1))) & 1;
	}
	inline int find_next(int i) const {
		if (i < 0) i = 0;
		for (int h = 0; h < lg; ++ h) {
			if ((i >> LOG_BUCKET) == int(data[h].size())) break;
			uint64_t d = data[h][i >> LOG_BUCKET] >> (i & (BUCKET - 1));
			if (!d) {
				i = (i >> LOG_BUCKET) + 1;
				continue;
			}
			i += low_bit(d);
			for (int g = h - 1; g >= 0; -- g) {
				i *= BUCKET;
				i += low_bit(data[g][i >> LOG_BUCKET]);
			}
			return i;
		}
		return N;
	}
	inline int find_prev(int i) const {
		if (i >= N) i = N - 1;
		for (int h = 0; h < lg; ++ h) {
			if (i == - 1) break;
			uint64_t d = data[h][i >> LOG_BUCKET] << (63 - (i & (BUCKET - 1)));
			if (!d) {
				i = (i >> LOG_BUCKET) - 1;
				continue;
			}
			i -= __builtin_clzll(d);
			for (int g = h - 1; g >= 0; -- g) {
				i *= BUCKET;
				i += top_bit(data[g][i >> LOG_BUCKET]);
			}
			return i;
		}
		return - 1;
	}
	template <typename F>
	inline void for_each(const int& L, const int& R, const F& fun) {
		for (int x = find_next(L); x < R; x = find_next(x + 1)) fun(x);
	}
};
