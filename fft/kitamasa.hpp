/**
 * KITASAMA !!
 *
 * * Computes the k-th term of a linear recurrence in O(N^2 log k).
 *
 * * Input:
 *   * A: first N terms, A[0..N-1]
 *   * C: recurrence coefficients of length N
 *        A[t] = C[0]*A[t-1] + C[1]*A[t-2] + ... + C[N-1]*A[t-N]
 *   * k: target index (0-based)
 * * Returns: A[k]
 *
 * * Usage: auto kth = kitamasa(init, bm, k);
 *
 * * Notes:
 *   * This implementation builds a polynomial `pol` such that:
 *         A[k] = sum_{i=0..N-1} pol[i] * A[i]
 *     using repeated squaring with reduction by the recurrence.
 *
 * * Requirements:
 *   * num must support +, -, *, / (field arithmetic, typically mod prime)
 *   * A.size() >= C.size()
 *
 * * Complexity:
 *   * Time :  O(N^2 log k)
 *   * Space: O(N)
 *
 * * TODO: If you need faster multiplication (e.g. O(N log N log k) or similar), use Bostan–Mori (with NTT/FFT) instead.
 */
template <typename Vector, typename T, typename U> static inline auto kitamasa(const Vector& A, const T& C, const U& k) {
	using num = typename Vector::value_type;
	int N = int(C.size());
	auto multiply = [&](auto& a, const bool& e) {
		// multiply (a ^ 2) * x ^ e
		std::vector<num> res(N << 1, num(0));
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				res[i + j + e] += a[i] * a[j];
			}
		}
		for (int i = (N << 1) - 1; i >= N; i--) {
			for (int j = 0; j < N; j++) {
				res[i - 1 - j] += res[i] * C[j];
			}
		}
		res.resize(N);
		a = std::move(res);
	};
	std::vector<num> pol(N); pol[0] = num(1);
	for (int i = std::__lg(k); i >= 0; i--) multiply(pol, (k >> i) & 1);
	num ans = num(0);
	for (int i = 0; i < N; i++) ans += pol[i] * A[i];
	return ans;
}
