/**
 * EULER TOUR EDGE !!!
 *
 * Usage: euler_tour_edge et(N); et.dfs(g, u); (can be multiple u)
 * After dfs() called, we can get the euler tour edge from euler and loc
 * (euler: the euler tour, loc: the first location of each node in the
 * euler tour). Euler tour is visited edges twice from root node and push
 * the nodes on that edges into vector euler.
 *
 * NOTE: Can apply lca with rmq: lca(u, v) = euler[rmq[loc[u], loc[v] + 1)]
**/
struct euler_tour_edge {
	std::vector<int> loc;
	std::vector<int> euler;

	template <typename Graph> inline void do_dfs(const Graph& g, const int& u) {
		loc[u] = int(euler.size());
		euler.push_back(u);
		for (const auto& id : g.adj[u]) {
			if (g.is_ignore(id)) continue;
			auto nxt = g(u, id);
			if (loc[nxt] != - 1) continue;
			do_dfs(g, nxt);
			euler.push_back(u);
		}
	}
	template <typename Graph> inline void dfs(const Graph& g, const int& u) {
		do_dfs(g, u);
		euler.push_back(-1);
	}

	inline void assign(const int& N, int v = -1) {
		loc.assign(N, v);
		euler.reserve(N << 1);
	}

	euler_tour_edge() = default;
	euler_tour_edge(const int &N) { this->assign(N); }
};
