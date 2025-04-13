/**
 * CONVEX_INCLUSION !!
 *
 * Return {inclusion, idx} of checking whether p in pts by binary search based on root (rt):
 *  * inclusion: -1 (inside), 0 (on 1 of edges), 1 (outside)
 *  * idx:(pts[rt], pts[idx], pts[idx + 1]) is the chosen triangle after binary search and if
 *    point is on the edge of polygon --> it's on the segment [pts[idx], pts[idx + 1])
 *
 * NOTE: num of points >= 3 and time complexity is log(N)
**/
template <typename Vector, typename Point>
static inline std::array<int, 2> convex_inclusion(const Point& p, const int& rt, const Vector& pts) {
	int N = int(pts.size());
#ifdef _GLIBCXX_DEBUG
	assert(N >= 3);
	assert(0 <= rt && rt < N);
#endif
	int mi = rt + 1;
	int ma = mi + N - 2;
	auto idx_norm = [&](const int& v) -> int { return v >= N ? v - N : v; };
	auto a = cross3(pts[rt], pts[idx_norm(mi)], p);
	if (a < 0) return {1, rt};
	int norm_ma = idx_norm(ma);
	auto b = cross3(pts[rt], pts[norm_ma], p);
	if (b > 0) return {1, norm_ma};
	while (mi + 1 < ma) {
		int md = (mi + ma) >> 1;
		if (cross3(pts[rt], pts[idx_norm(md)], p) >= 0) mi = md;
		else ma = md;
	}
	int norm_mi = idx_norm(mi);
	norm_ma = idx_norm(ma);
	auto c = cross3(pts[norm_mi], pts[norm_ma], p);
	if (c < 0) return {1, norm_mi};
	if (c == 0) return {0, norm_mi};
	if (mi == rt + 1 && a == 0) return {0, rt};
	if (ma == rt + N - 1 && b == 0) return {0, norm_ma};
	return {-1, norm_mi};
}
