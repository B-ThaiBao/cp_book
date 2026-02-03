/**
 * BERLEKAMP MASSEY (BM) !!
 *
 * * Finds the shortest linear recurrence for a given sequence A.
 *
 * * Returns a vector C of length L such that:
 *    A[i] = C[0]*A[i-1] + C[1]*A[i-2] + ... + C[L-1]*A[i-L]
 *
 * * Usage: auto C = berlekamp_massey(A);
 *
 * * Requirements:
 *   * value_type must support +, -, *, / and equality comparison
 *   * arithmetic must be over a field (e.g. modulo prime)
 *
 * * Complexity:
 *   * Time : O(N^2) in the worst case
 *   * Space: O(N)
 */
template <typename Vector> static inline Vector berlekamp_massey(const Vector& A) {
	using num = typename Vector::value_type;
	int N = int(A.size()), L = 0;
	std::vector<num> B, C, tmp;
	B.reserve(N + 1), C.reserve(N + 1), tmp.reserve(N + 1);
	B.push_back(num(1)), C.push_back(num(1));

	num b = num(1);
	for (int i = 0, m = 1; i < N; i++, m++) {
		num d = A[i];
		for (int j = 1; j <= L; j++) d += C[j] * A[i - j];
		if (d == num(0)) continue;
		num coef = d / b;

		int M = int(B.size()), P = int(C.size());
		if ((L << 1) <= i) {
			tmp = C;
			if (P < M + m) C.resize(M + m, num(0));
			for (int k = 0; k < M; k++) C[k + m] -= coef * B[k];
			L = i + 1 - L, B = std::move(tmp), b = d, m = 0;
		} else {
			if (P < M + m) C.resize(M + m, num(0));
			for (int k = 0; k < M; k++) C[k + m] -= coef * B[k];
		}
	}
	for (int i = 0; i < L; i++) C[i] = -C[i + 1];
	C.resize(L);
	return C;
}
