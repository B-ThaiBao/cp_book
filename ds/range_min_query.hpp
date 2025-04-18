/**
 * Mostly inspired here: https://codeforces.com/blog/entry/78931
 *
 * Usage: range_min_query<T> rmq_min or range_max_query<T> rmq_max
 *   * Build from array A: build(A) or build(std::move(A))
 *   * Get query: range_query(L, R) --> return pair {idx, value (T)}
 *
 * Return min/max query in O(1) with precomputed in O(N) by principles:
 *   * Leftmost min query : std::less<>
 *   * Rightmost min query: std::less_equal<>
 *   * Leftmost max query : std::greater<>
 *   * Rightmost max query: std::greater_equal<>
**/
template <typename T, typename Comp = std::less<T>> struct range_min_query {
	static constexpr int32_t BUCKET_SIZE = 32;
	static constexpr int32_t BUCKET_SIZE_LOG = 5;

	int32_t N;
	std::vector<T> data;
	std::vector<uint32_t> mask;
	std::vector<int32_t> sparse;
	Comp comp;

	static constexpr int32_t floor_log_2(const uint32_t& x) { return 31 - __builtin_clz(x); }
	static constexpr int32_t least_bit(const uint32_t& x) { return __builtin_ctz(x); }
	inline int32_t index(const int32_t& i, const int32_t& j) const {
		return comp(data[j], data[i]) ? j : i;
	}
	int32_t small_range_query(const int32_t& R, const int32_t& num = BUCKET_SIZE) const {
		return R - floor_log_2(mask[R] << (BUCKET_SIZE - num) >> (BUCKET_SIZE - num));
	}

	void build_mask_and_sparse() {
		mask.assign(N, uint32_t());
		sparse.assign(N, int32_t());
		for (uint32_t i = 0, num = 0; i < uint32_t(N); mask[i++] = (num |= 1)) {
			num <<= 1;
			while (num > 0 && comp(data[i], data[i - least_bit(num)])) {
				num ^= (num & (- num));
			}
		}
		// TODO: Try to build a block sparse table
		int32_t num_blocks = N >> BUCKET_SIZE_LOG;
		for (int32_t i = 0; i < num_blocks; ++i) {
			sparse[i] = small_range_query(BUCKET_SIZE * i + BUCKET_SIZE - 1);
		}
		for (int32_t j = 1; (1 << j) <= num_blocks; ++j) {
			for (int32_t i = 0; i + (1 << j) <= num_blocks; i ++) {
				sparse[num_blocks * j + i] = index(sparse[num_blocks * (j - 1) + i], sparse[num_blocks * (j - 1) + i + (1 << (j - 1))]);
			}
		}
	}

	range_min_query() {}
	explicit range_min_query(std::vector<T>&& v, const Comp& com = Comp()) : comp(com) {
		N = int32_t(v.size());
		data = std::forward<std::vector<T>>(v);
		this->build_mask_and_sparse();
	}
	template <typename Vec> explicit range_min_query(const Vec& v, const Comp& com = Comp()) : comp(com) {
		N = int32_t(v.size());
		data = v;
		this->build_mask_and_sparse();
	}
	void build(std::vector<T>&& v, const Comp& com = Comp()) {
		N = int32_t(v.size());
		data = std::forward<std::vector<T>>(v);
		this->comp = com;
		this->build_mask_and_sparse();
	}
	template <typename Vec> void build(const Vec& v, const Comp& com = Comp()) {
		N = int32_t(v.size());
		data = v;
		this->comp = com;
		this->build_mask_and_sparse();
	}

	// Return the index and the value of the minimum element in [L, R)
	std::pair<int32_t, T> range_query(const int32_t& L, const int32_t& R) const {
#ifdef _GLIBCXX_DEBUG
		assert(0 <= L && L < R && R <= N);
#endif
		if (R - L <= BUCKET_SIZE) {
			int32_t idx = small_range_query(R - 1, R - L);
			return std::make_pair(idx, data[idx]);
		}

		int32_t res = small_range_query(L + BUCKET_SIZE - 1);
		int32_t bucket_l = (L >> BUCKET_SIZE_LOG) + 1, bucket_r = ((R - 1) >> BUCKET_SIZE_LOG) - 1;
		if (bucket_l <= bucket_r) {
			const int32_t j = floor_log_2(bucket_r - bucket_l + 1);
			const int32_t num_blocks = N >> BUCKET_SIZE_LOG;
			res = index(res, sparse[num_blocks * j + bucket_l]);
			res = index(res, sparse[num_blocks * j + bucket_r - (1 << j) + 1]);
		}
		res = index(res, small_range_query(R - 1));
		return std::make_pair(res, data[res]);
	}
};
