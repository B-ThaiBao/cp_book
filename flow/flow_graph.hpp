template <typename flow_t> struct flow_graph {
	static constexpr flow_t eps = flow_t(1e-9);
	struct edge_t {
		int from;
		int to;
		flow_t cap;
		flow_t flow;
	};
	int V;
	std::vector<std::vector<int>> adj;
	std::vector<edge_t> edges;

	flow_graph() {}
	flow_graph(const int& N) : V(N), adj(N) {}
	flow_graph(const int& N, const int& M) : V(N), adj(N) {
		edges.reserve(M << 1); // Avoid relocation
	}

	int add_edge(const int& frm, const int& to, const flow_t& forw, const flow_t& back = 0) {
		int e = int(edges.size());
		edges.emplace_back(edge_t{frm, to, forw, 0});
		edges.emplace_back(edge_t{to, frm, back, 0});
		adj[frm].push_back(e);
		adj[to].push_back(e ^ 1);
		return e;
	}
	void clear_flow() {
		for (auto& e : edges) e.flow = 0;
	}
};

template <typename flow_t, typename cost_t> struct cost_flow_graph {
	static constexpr flow_t eps = flow_t(1e-9);
	struct edge_t {
		int from;
		int to;
		flow_t cap;
		flow_t flow;
		cost_t cost;
	};
	int V;
	std::vector<std::vector<int>> adj;
	std::vector<edge_t> edges;

	cost_flow_graph() {}
	cost_flow_graph(const int& N) : V(N), adj(N) {}
	cost_flow_graph(const int& N, const int& M) : V(N), adj(N) {
		edges.reserve(M << 1); // Avoid relocation
	}

	int add_edge(const int& frm, const int& to, const flow_t& forw, const flow_t& back, const cost_t& cost) {
		int e = int(edges.size());
		edges.emplace_back(edge_t{frm, to, forw, 0, cost});
		edges.emplace_back(edge_t{to, frm, back, 0, - cost});
		adj[frm].push_back(e);
		adj[to].push_back(e ^ 1);
		return e;
	}
	void clear_flow() {
		for (auto& e : edges) e.flow = 0;
	}
};
