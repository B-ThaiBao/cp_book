/**
 * CONVEX HULL !!!
 *
 * Usage: auto ch = convex_hull(pts, order);
 *  * pts  : vector of pair_point<T>
 *  * order: order of points sorted by indicies (use std::sort or std::stable_sort, etc ...)
 *
 * NOTE: Returns 2 vectors:
 *  * ch[0]: the counterclockwise order of polygon (include colinear points and no same points)
 *  * ch[1]: is the type of each point based on their indicies
 *   - ch[1][i] == i: This is vertex of polygon (most important one)
 *   - ch[1][i] < N && ch[1][i] != i: This is point that same with another vertex of polygon
 *   - ch[1][i] == N: This is point that only in edges of polygon (colinear point)
 *   - ch[1][i] == N + 1: This is point that not on the polygon
**/
template <typename Vector, typename Order>
static inline std::array<std::vector<int>, 2> convex_hull(const Vector& pts, const Order& order) {
	int N = int(pts.size());
	if (N == 0) return {};
	std::array<std::vector<int>, 2> res;
	res[0].reserve(N + 1); res[1].resize(N, 0);
	res[1][order[0]] = N + 1;
	for (int i = 1; i < N; i++) {
		int prv = order[i - 1], cur = order[i];
		if (pts[cur] == pts[prv]) {
			res[1][cur] = res[1][prv] == N + 1 ? prv : res[1][prv];
		} else {
			res[1][cur] = N + 1;
		}
	}
	// Build upper convex hull
	for (int i = N - 1; i >= 0; --i) {
		int p = order[i];
		if (res[1][p] != N + 1) continue;
		while (int(res[0].size()) > 1) {
			auto c = cross3(pts[res[0].end()[-2]], pts[res[0].back()], pts[p]);
			if (c < 0) {
				res[1][res[0].back()] = N + 1;
				res[0].pop_back();
			} else if (c == 0) {
				// This point is on edge of the polygon but not its vertex
				res[1][res[0].back()] = N;
				break;
			} else break;
		}
		res[1][p] = p;
		res[0].push_back(p);
	}
	res[1][res[0].back()] = N + 1;
	res[1][res[0][0]] = N + 1;
	if (int(res[0].size()) > 1) res[0].pop_back();
	// Build lower convex hull
	auto st = int(res[0].size());
	for (int i = 0; i < N; ++i) {
		int p = order[i];
		if (res[1][p] != N + 1) continue;
		while (int(res[0].size()) > st + 1) {
			auto c = cross3(pts[res[0].end()[-2]], pts[res[0].back()], pts[p]);
			if (c < 0) {
				res[1][res[0].back()] = N + 1;
				res[0].pop_back();
			} else if (c == 0) {
				// This point is on edge of the polygon but not its vertex
				res[1][res[0].back()] = N;
				break;
			} else break;
		}
		res[1][p] = p;
		res[0].push_back(p);
	}
	if (int(res[0].size()) > st + 1) res[0].pop_back();
	if (int(res[0].size()) == 2 && res[0][0] == res[0][1]) res[0].pop_back();
	for (int i = 0; i < N; i++) {
		auto& r = res[1][i];
		if (r < N && r != i) r = res[1][r];
	}
	return res;
}
