template <typename T> std::vector<int> topo_sort(const T& g) {
	// TODO: Try to sort node in the topo_order based on BFS
	// If the graph is not DAG, return a empty vector of nodes
	std::vector<int> deg(g.V, 0);
	for (int i = 0; i < int(g.edges.size()); i ++) {
		if (g.is_ignore(i)) continue;
		++ deg[g.edges[i].to];
	}
	std::vector<int> order; order.reserve(g.V);
	for (int i = 0; i < g.V; i ++) {
		if (deg[i] == 0) order.push_back(i);
	}
	for (int beg = 0; beg < int(order.size()); beg ++) {
		int i = order[beg];

		for (const auto& id : g.adj[i]) {
			if (g.is_ignore(id)) continue;
			const auto& e = g.edges[id];
			if (-- deg[e.to] == 0) {
				order.push_back(e.to);
			}
		}
	}

	if ((int) order.size() != g.V) return std::vector<int>();
	return order;
}