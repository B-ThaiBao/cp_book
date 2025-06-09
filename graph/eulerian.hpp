template <typename G> static inline std::pair<int, std::vector<int>> find_eulerian_path(const G& g) {
	std::vector<std::array<int, 2>> deg(g.V);
	int num_edges = int(g.edges.size());
	for (int i = 0; i < num_edges; i++) {
		const auto& e  = g.edges[i];
		deg[e.from][0]++, deg[e.to][1]++;
	}
	int rt = -1, o = 0;
	for (int i = 0; i < g.V; i++) {
		if (int(deg[i][0] + deg[i][1]) & 1) {
			o++;
			if (rt == -1 || deg[i][0] - deg[i][1] > deg[rt][0] - deg[rt][1]) {
				rt = i;
			}
		}
	}
	if (o > 2) return {-1, std::vector<int>()};
	if (rt == -1) {
		rt = 0;
		while (rt < g.V && deg[rt][0] + deg[rt][1] == 0) rt++;
		if (rt == g.V) return {0, std::vector<int>()};
	}

	std::vector<bool> vis(num_edges, false);
	std::vector<int> c(g.V, 0);
	std::vector<int> p(g.V, 0);
	std::vector<int> ans(num_edges);
	int stk = 0, idx = num_edges, v = rt;
	while (true) {
		bool seen = false;
		while (c[v] < int(g.adj[v].size())) {
			int x = g.adj[v][c[v]++];
			if (vis[x]) continue;
			vis[x] = true;
			ans[stk++] = x;
			const auto& e = g.edges[x];
			p[v]++;
			v = g(v, e);
			p[v]--;
			seen = true;
			break;
		}
		if (!seen) {
			if (stk == 0) break;
			int x = ans[--stk];
			ans[--idx] = x;
			const auto& e = g.edges[x];
			v = g(v, e);
		}
	}

	int tot = 0; for (const auto& x : p) tot += std::abs(x);
	if (idx != 0 || tot > 2) return {-1, std::vector<int>()};
	return {rt, std::move(ans)};
}
