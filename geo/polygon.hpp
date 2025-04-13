/**
 * Return doubled oriented area (2S)
 * Oriented means that it's positive if polygon vertices are listed in
 * counter-clockwise (ccw) order and negative otherwise.
**/
template <typename T, typename Vector> static inline T area(const Vector& pts) {
	int N = int(pts.size());
	if (N == 0) return 0;
	T res = 0;
	for (int i = 0; i < N - 1; ++i) res += cross(pts[i], pts[i + 1]);
	res += cross(pts[N - 1], pts[0]);
	return res;
}

/**
 * Return centroid of a (possibly non-convex) polygon
 * Assuming that the coordinates are listed in a clockwise orcounterclockwise fashion
 * NOTE: Centroid often known as the "center of gravity" or "center of mass".
**/
template <typename T, typename Vector> static inline pair_point<T> centroid(const Vector& pts) {
#ifdef _GLIBCXX_DEBUG
	assert(!pts.empty());
#endif
	int N = int(pts.size());
	T tot = T(0);
	for (int i = 0; i < N - 1; ++i) tot += cross(pts[i], pts[i + 1]);
	tot += cross(pts[N - 1], pts[0]);
	T scale = T(3) * tot;
	pair_point<T> res(0, 0);
	for (int i = 0; i < N - 1; i++) {
		res = res + (pts[i] + pts[i + 1]) * cross(pts[i], pts[i + 1]);
	}
	res = res + (pts[N - 1] + pts[0]) * cross(pts[N - 1], pts[0]);
	return res / scale;
}
