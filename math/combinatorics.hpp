// NOTE: mnum must be inversible (to apply this one).
// If not, try another ways (i.e: Pascal triangle, Lucas theorem, ...)
template <typename mnum> struct combinatorics {
	static std::vector<mnum> factor;
	static std::vector<mnum> ifactor;

	static inline void ensure_fact(const int& N) {
		int M = int(factor.size()) - 1;
		if (N <= M) return;

		factor.resize(N + 1); ifactor.resize(N + 1);
		for (int i = M + 1; i <= N; ++ i) {
			factor[i] = factor[i - 1] * i;
		}
		ifactor[N] = mnum(1) / factor[N];
		for (int i = N - 1; i > M; -- i) {
			ifactor[i] = ifactor[i + 1] * (i + 1);
		}
	}

	combinatorics() {}
	combinatorics(const int& N) { ensure_fact(N); }

	static inline mnum fact(const int& N) {
		if (N < 0) return mnum(0);
		ensure_fact(N); return factor[N];
	}
	static inline mnum inv_fact(const int& N) {
		if (N < 0) return mnum(0);
		ensure_fact(N); return ifactor[N];
	}

	template <typename T, typename P>
	static inline mnum choose(const T& K, const P& N) {
		if (K < 0 || K > N) return mnum(0);
		ensure_fact(N);
		return factor[N] * ifactor[K] * ifactor[N - K];
	}
	template <typename T, typename P>
	static inline mnum inv_choose(const T& K, const P& N) {
		if (K < 0 || K > N) return mnum(0);
		ensure_fact(N);
		return ifactor[N] * factor[K] * factor[N - K];
	}

	template <typename T, typename P>
	static inline mnum permute(const T& K, const P& N) {
		if (K < 0 || K > N) return mnum(0);
		ensure_fact(N);
		return factor[N] * ifactor[N - K];
	}
	template <typename T, typename P>
	static inline mnum inv_permute(const T& K, const P& N) {
		if (K < 0 || K > N) return mnum(0);
		ensure_fact(N);
		return ifactor[N] * factor[N - K];
	}

	static inline mnum catalan(const int& N) {
		if (N < 0) return mnum(0);
		if (N == 0) return mnum(1);
		ensure_fact(N << 1);
		return factor[N << 1] * ifactor[N + 1] * ifactor[N];
	}
	static inline mnum inv_catalan(const int& N) {
		if (N < 0) return mnum(0);
		if (N == 0) return mnum(1);
		ensure_fact(N << 1);
		return ifactor[N << 1] * factor[N + 1] * factor[N];
	}
};
template <typename mnum> std::vector<mnum> combinatorics<mnum>::factor
						= std::vector<mnum>(1, mnum(1));
template <typename mnum> std::vector<mnum> combinatorics<mnum>::ifactor
						= std::vector<mnum>(1, mnum(1));