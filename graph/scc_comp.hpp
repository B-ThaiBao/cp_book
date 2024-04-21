template <typename T> std::vector<int> scc_comp(const T& g, int& cnt) {
	std::vector<int> low(g.V, - 1), dfn(g.V, - 1), c(g.V, - 1);
	std::vector<int> stk; stk.reserve(g.V);
	int time = 0; cnt = 0;

	auto dfs = [&](auto&& dfs, const int& u) -> void {
		dfn[u] = low[u] = time ++;
		stk.push_back(u);
		for (const auto& id : g.adj[u]) {
			if (g.is_ignore(id)) continue;
			const auto& e = g.edges[id];
			if (dfn[e.to] == - 1) {
				dfs(dfs, e.to);
				low[u] = std::min(low[u], low[e.to]);
			} else if (c[e.to] == - 1) {
				low[u] = std::min(low[u], dfn[e.to]);
			}
		}

		if (dfn[u] == low[u]) {
			int v;
			do {
				v = stk.back();
				stk.pop_back();
				c[v] = cnt;
			} while (v != u);
			++ cnt;
		}
	};

	for (int i = 0; i < g.V; i ++) {
		if (dfn[i] == - 1) dfs(dfs, i);
	}
	for (auto& x : c) x = cnt - 1 - x;
	return c;
	// c[i] <= c[j] for every edge i -> j
	// c[i] == c[j]: i and j are in same_comp
	// NOTE: Consider c[i] is condensation graph, the topo_sort of
	// that graph is: 0 1 2 ... (cnt - 1)
}