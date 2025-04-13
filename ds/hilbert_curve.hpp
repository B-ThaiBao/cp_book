namespace hilbert_curve {

static inline int floor_log_2(const uint64_t& N) { return 63 - __builtin_clzll(N); }

static inline uint64_t hilbert_order(uint64_t x, uint64_t y) {
	const uint64_t logn = floor_log_2(std::max(x, y) << 1 | 1) | 1;
	const uint64_t maxn = (1LLU << logn) - 1;
	uint64_t res = 0;
	for (uint64_t s = 1LLU << (logn - 1); s != 0; s >>= 1) {
		bool rx = x & s, ry = y & s;
		res = (res << 2) | (rx ? ry ? 2 : 1 : ry ? 3 : 0);
		if (!rx) {
			if (ry) x ^= maxn, y ^= maxn;
			std::swap(x, y);
		}
	}
	return res;
}

template <typename Begin, typename End, typename P, typename Q, typename R, typename S, typename F>
static inline void for_each(Begin begin, End end, const P& add_left, const Q& add_right, const R& delete_left, const S& delete_right, F&& f) {
	if (begin == end) return;
	// At every step, we maintain [lo, hi)
	int lo = begin->L, hi = (begin->L);
	for (auto it = begin; it != end; ++it) {
		while (lo > it->L) --lo, add_left(lo);
		while (hi < it->R) add_right(hi), ++hi;
		while (lo < it->L) delete_left(lo), ++lo;
		while (hi > it->R) --hi, delete_right(hi);
		f(*it);
	}
}

} // namespace hilbert_curve
