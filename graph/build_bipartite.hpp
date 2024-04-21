template <typename T> std::vector<int8_t> build_bipartite(const T& g) {
	const auto& N = g.V;
	std::vector<int8_t> side(N, 2);
	std::vector<int> q; q.reserve(N);
	for (int u = 0; u < N; ++ u) {
		if (side[u] == 2) {
			q.push_back(u);
			side[u] = 0;
			for (int beg = 0; beg < int(q.size()); beg ++) {
				const auto& node = q[beg];
				for (const auto& id : g.adj[node]) {
					if (g.is_ignore(id)) continue;
					auto v = g(node, id);
					if (side[v] == 2) {
						side[v] = side[node] ^ 1;
						q.push_back(v);
					} else if (side[node] == side[v]) {
						// Graph doesn't satisfied. We return empty vector
						return std::vector<int8_t>();
					}
				}
			}
		}
	}
	return side;
}