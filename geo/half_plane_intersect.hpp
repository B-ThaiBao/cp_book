/**
 * HALF PLANE INTERSECT !!!
 *
 * Usage: Return std::vector means that all edge indexs of half plane intersect in ccw order
 * means that unbounded region can also return some edges if neccesary. To check unbounded region:
 *  * sz(hp) < 3 || cross(dir(front), dir(back)) >= 0: unbounded region
 *  * Otherwise, it is a bounded region
 *
 * Properties:
 *  * Visibility region of abitary polygon.
 *  * Find biggest radius of circle inside convex polygon by binary search and move half plane by r.
 *  * Solve 2D LP: find feasible region --> binary search based on targent vector of LP objective.
 *
 * NOTE: Please sort the edge indexs in ccw order before call half_plane_intersect(pls, order)
 *       Type of line must hold up to x ^ 4 of point coordinates
**/
template <typename Vector, typename Order>
static inline std::vector<int> half_plane_intersect(const Vector& pls, const Order& order) {
	std::deque<int> q;
	auto is_not_good = [&](const auto& u, const auto& v, const auto& t) -> bool {
		// t contains the intersection point of u and v
		auto du = dir(u), dv = dir(v), dt = dir(t);
		auto s = cross(du, dv);
		if (s == 0) return false;
		auto p = cross(dt, u[0] * s + du * cross(v[0] - u[0], dv) - t[0] * s);
		return s > 0 ? p < 0 : p > 0;
	};
	for (const auto& o : order) {
		while (int(q.size()) > 1 && is_not_good(pls[q.end()[-2]], pls[q.back()], pls[o])) q.pop_back();
		while (int(q.size()) > 1 && is_not_good(pls[q.front()], pls[q[1]], pls[o])) q.pop_front();
		// Now we need to handle with special case: parallel half planes
		if (int(q.size()) > 0) {
			auto du = dir(pls[q.back()]), dv = dir(pls[o]);
			if (cross(du, dv) == 0) {
				if (dot(du, dv) > 0) {
					// Same direction
					if (cross(du, pls[o][0] - pls[q.back()][0]) > 0) q.back() = o;
				} else {
					// Opposite direction
					auto c = cross(du, pls[o][0] - pls[q.back()][0]);
					if (c <= 0) q.push_back(o);
					else q.pop_back();
				}
			} else q.push_back(o);
		} else q.push_back(o);
	}
	while (int(q.size()) > 2 && is_not_good(pls[q.end()[-2]], pls[q.back()], pls[q.front()])) q.pop_back();
	while (int(q.size()) > 2 && is_not_good(pls[q.front()], pls[q[1]], pls[q.back()])) q.pop_front();
	return {q.begin(), q.end()};
}
