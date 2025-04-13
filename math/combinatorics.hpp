/**
 * COMBINATORICS:
 *  * NOTE: mnum must be inversible (to apply this one).
 *  * If not, try another ways (i.e: Pascal triangle, Lucas theorem, ...)
 *
 * CATALAN NUMBERS: Counting problem has a solution expressed by combining two
 * identical subproblems with total size = N - 1, it means this problem has the
 * same solution as all other catalan problems.
 *  * Formula based on dp:
 *     * Base case : dp(0) = dp(1) = 1
 *     * Transition: dp(N) = sum(dp(K) * dp(N - 1 - K)) with 0 <= K < N
 *                or dp(N) = sum(dp(a) * dp(b)) with all a + b == N - 1
 *
 *  * To solve the probs with prefix known beforehand, just image the grid:
 *     * After prefix, assume that the point is (x, y) (without prefix it is (0, 0))
 *     * Total path is: (x, y) -> (N, M); invalid path is: (x, y) -> (N + 1, M - 1)
 *     ---> Valid path is: choose(N + M - x - y, N - x) - choose(N + M - x - y, N + 1 - x)
 *
 * CATALAN CONVOLUTIONS: sum of products of M catalan numbers where the
 * index of product sums up to S.
 *  * Formula based on Catalan numbers:
 *     * dp(M, S) = sum(C(A[1]) * C(A[2]) * ... * C(A[M])) with A[1] + A[2] + ... + A[M] == S
 *
 * Usage:
 *   * catalan(N, M, K): num ways to place N "+1" s and M "-1" s starting with K "+1" at prefix
 *                       in such that all vals of prefix sum >= 0 (+1 always >= -1)
 *   * catalan(S, S + M - 1, M - 1): catalan convolutions dp(M, S)
 *   * Besides catalan ways, you can change the problems into grid problems like above.
**/
template <typename num> struct combinatorics : std::vector<std::array<num, 2>> {
	combinatorics() {}
	combinatorics(const int& N) : std::vector<std::array<num, 2>>(N) {
		this->at(0)[0] = num(1);
		for (int i = 1; i < N; ++i) this->at(i)[0] = this->at(i - 1)[0] * num(i);
		this->at(N - 1)[1] = num(1) / this->at(N - 1)[0];
		for (int i = N - 1; i >= 1; --i) this->at(i - 1)[1] = this->at(i)[1] * num(i);
	}

	template <typename T> inline num fact(const T& N) const {
		if (N < 0) return num(0);
#ifdef _GLIBCXX_DEBUG
		assert(0 <= int(N) && int(N) < int(this->size()));
#endif
		return this->at(N)[0];
	}
	template <typename T> inline num inv_fact(const T& N) const {
		if (N < 0) return num(0);
#ifdef _GLIBCXX_DEBUG
		assert(0 <= int(N) && int(N) < int(this->size()));
#endif
		return this->at(N)[1];
	}

	template <typename T, typename U> inline num choose(const T& N, const U& R) const {
		if (R < 0 || R > N) return num(0);
#ifdef _GLIBCXX_DEBUG
		assert(0 <= int(R) && int(R) <= int(N) && int(N) < int(this->size()));
#endif
		return this->at(N)[0] * this->at(R)[1] * this->at(N - R)[1];
	}
	template <typename T, typename U> inline num inv_choose(const T& N, const U& R) const {
		if (R < 0 || R > N) return num(0);
#ifdef _GLIBCXX_DEBUG
		assert(0 <= int(R) && int(R) <= int(N) && int(N) < int(this->size()));
#endif
		return this->at(N)[1] * this->at(R)[0] * this->at(N - R)[0];
	}

	template <typename T, typename U> inline num permute(const T& N, const U& R) const {
		if (R < 0 || R > N) return num(0);
#ifdef _GLIBCXX_DEBUG
		assert(0 <= int(R) && int(R) <= int(N) && int(N) < int(this->size()));
#endif
		return this->at(N)[0] * this->at(N - R)[1];
	}
	template <typename T, typename U> inline num inv_permute(const T& N, const U& R) const {
		if (R < 0 || R > N) return num(0);
#ifdef _GLIBCXX_DEBUG
		assert(0 <= int(R) && int(R) <= int(N) && int(N) < int(this->size()));
#endif
		return this->at(N)[1] * this->at(N - R)[0];
	}

	template <typename T> inline num catalan(const T& N) const {
		if (N < 0) return num(0);
#ifdef _GLIBCXX_DEBUG
		assert(0 <= int(N) && int(N << 1) < int(this->size()) && int(N + 1) < int(this->size()));
#endif
		return this->at(N << 1)[0] * this->at(N + 1)[1] * this->at(N)[1];
	}
	template <typename T> inline num inv_catalan(const T& N) const {
		if (N < 0) return num(0);
#ifdef _GLIBCXX_DEBUG
		assert(0 <= int(N) && int(N << 1) < int(this->size()) && int(N + 1) < int(this->size()));
#endif
		return this->at(N << 1)[1] * this->at(N + 1)[0] * this->at(N)[0];
	}

	template <typename T, typename U> inline num catalan(const T& N, const U& M) const {
		if (M > N) return num(0);
		return this->choose(N + M, M) - this->choose(N + M, M - 1);
	}

	template <typename T, typename U, typename V>
	inline num catalan(const T& N, const U& M, const V& K) const {
		if (M > N + K) return num(0);
		return this->choose(N + M, M) - this->choose(N + M, M - K - 1);
	}
};
