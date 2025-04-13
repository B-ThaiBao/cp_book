// Babylonian sqrt method: https://www.codeabbey.com/index/wiki/square-root-approximation
template <typename T> static constexpr T floor_sqrt(const T& N) {
#ifdef _GLIBCXX_DEBUG
	assert(N >= 0);
#endif
	if (N == 0) return 0;
	T a = N;
	while (true) {
		T b = N / a;
		if (a - b <= 1) return b;
		a = (a + b + 1) >> 1;
	}
}

template <typename T> static constexpr T ceil_sqrt(const T& N) {
#ifdef _GLIBCXX_DEBUG
	assert(N >= 0);
#endif
	if (N == 0) return 0;
	T a = N;
	while (true) {
		T b = N / a;
		if (a - b <= 1) return b * b < N ? b + 1 : b;
		a = (a + b + 1) >> 1;
	}
}
