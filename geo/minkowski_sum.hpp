/**
 * MINKOWSKI SUM !!
 *
 * Return minkowski sum of two convex polygons A and B with the given point
 *
 * NOTE:
 *  * num_points = num_edge_vectors = tot all distinct vectors in all polygons (multiple)
 *  * The result polygon is start based on the principle of given point such as
 *    top leftmost, top rightmost, bottom leftmost, bottom rightmost, ...
 *  * Every vertex of minkowski sum is sum of vertices of A and B and vice versa
**/
template <typename T, typename P> static inline auto minkowski_sum(const T& A, const P& B, int i, int j) {
	int a = 0, b = 0;
	int N = int(A.size()), M = int(B.size());
#ifdef _GLIBCXX_DEBUG
	assert(0 <= i && i < N);
	assert(0 <= j && j < M);
#endif
	std::vector<decltype(A[i] + B[j])> res; res.reserve(N + M);
	while (a < N || b < M) {
		res.push_back(A[i] + B[j]);
		int ni = i + 1; if (ni == N) ni = 0;
		int nj = j + 1; if (nj == M) nj = 0;
		auto c = cross(A[ni] - A[i], B[nj] - B[j]);
		if (c >= 0) i = ni, ++a;
		if (c <= 0) j = nj, ++b;
	}
	return res;
}
