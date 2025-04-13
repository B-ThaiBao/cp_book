/**
 * TREE DIAMETER !!
 *
 * Usage:
 *   * Return pair of two nodes {x, y} that are farthest from each other (diameter).
 *   * After bfs() called, we receive the rooted tree at x with dist, pe
 *
 * Properties of diameter:
 *   * The farthest node from each node is one end of the diameter
 *   * The height of each component is smaller than the distance to the closest end of the diameter
**/
template <typename T> struct tree_diameter {
	std::vector<T> dist;
	std::vector<int> pe;
	std::vector<int> q;

	tree_diameter() = default;
	tree_diameter(const int& N) : dist(N, -1), pe(N, -1), q(N) {}

	template <typename Graph>
	std::pair<int, int> bfs(const Graph& g, const int& s) {
#ifdef _GLIBXX_DEBUG
		assert(dist[s] == -1);
#endif
		auto do_bfs = [&](const int& u) -> int {
			dist[u] = T(0);
			pe[u] = -1;
			q[0] = u;
			int ma = u;
			for (int b = 0, e = 1; b < e; ++b) {
				const int& v = q[b];
				for (const auto& id : g.adj[v]) {
					if (g.is_ignore(id)) continue;
					if (id == pe[v]) continue;
					int to = g(v, id);
					dist[to] = dist[v] + g.edges[id].cost;
					if (dist[to] > dist[ma]) ma = to;
					pe[to] = id;
					q[e++] = to;
				}
			}
			return ma;
		};
		int x = do_bfs(s);
		int y = do_bfs(x);
		return {x, y};
	}
};
