template <typename T, typename U> static inline T legendre_power(T N, const U& P) {
#ifdef _GLIXX_DEBUG
	assert(P > 0 && N >= 0);
#endif
	T res = 0;
	while (N > 0) {
		N /= P;
		res += N;
	}
	return res;
}
