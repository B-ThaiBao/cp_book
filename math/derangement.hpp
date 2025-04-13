// Count derangements [0, N) by recurive formula: dp[i] = (i - 1) * (dp[i - 1] + dp[i - 2])
// Beware of dp[0], in this case assuming dp[0] = 1
template <typename T> static inline std::vector<T> count_derangement(const int& N) {
#ifdef _GLIBCXX_DEBUG
	assert(N >= 0);
#endif
	if (N == 0) return {};
	if (N == 1) return {T(0)};
	std::vector<T> dp(N);
	dp[0] = T(1); dp[1] = T(0);
	for (int i = 2; i < N; ++i) dp[i] = T(i - 1) * (dp[i - 1] + dp[i - 2]);
	return dp;
}
