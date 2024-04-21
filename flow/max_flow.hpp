// Time: O(V ^ 2 * E)
template <typename flow_t> struct dinic_max_flow {
	std::vector<int> ptr;
	// Assume that network layer is a tree with root is sink
	// depth[sink] = 0 and other node can be reached by back_edges
	std::vector<int> depth;
	std::vector<int> q; // Instead of std::queue

	inline void init(const int& N) {
		ptr.resize(N);
		depth.resize(N);
		q.resize(N);
	}

	dinic_max_flow() {}
	dinic_max_flow(const int& N) { init(N); }

	template <typename Network>
	bool bfs(const Network& g, const int& s, const int& t) {
		std::fill(depth.begin(), depth.end(), - 1); q[0] = t;
		depth[t] = 0;
		for (int beg = 0, end = 1; beg < end; ++ beg) {
			int i = q[beg];
			for (const int& id : g.adj[i]) {
				const auto& e = g.edges[id];
				const auto& back = g.edges[id ^ 1];
				if (back.cap - back.flow > g.eps && depth[e.to] == - 1) {
					depth[e.to] = depth[i] + 1;
					if (e.to == s) {
						// We have the path from s to t
						return true;
					}
					q[end ++] = e.to;
				}
			}
		}
		return false;
	}

	template <typename Network>
	flow_t dfs(Network& g, const int& v, const int& t, const flow_t& w) {
		if (v == t) return w;
		for (; ptr[v] >= 0; -- ptr[v]) {
			const int id = g.adj[v][ptr[v]];
			const auto& e = g.edges[id];
			if (e.cap - e.flow > g.eps && depth[e.to] == depth[v] - 1) {
				flow_t min_flow = dfs(g, e.to, t, std::min<flow_t>(e.cap - e.flow, w));
				if (min_flow > g.eps) {
					g.edges[id].flow += min_flow;
					g.edges[id ^ 1].flow -= min_flow;
					return min_flow;
				}
			}
		}
		return 0;
	}

	template <typename Network>
	flow_t max_flow(Network& g, const int& s, const int& t) {
		flow_t tot_flow = 0;
		while (bfs(g, s, t)) {
			for (int i = 0; i < g.V; ++ i) {
				ptr[i] = int(g.adj[i].size()) - 1;
			}
			flow_t cur_flow = 0;
			// Dinic process
			while (true) {
				flow_t cur_dinic_flow = dfs(g, s, t, std::numeric_limits<flow_t>::max());
				if (cur_dinic_flow <= g.eps) break;
				cur_flow += cur_dinic_flow;
			}
			if (cur_flow <= g.eps) break;
			tot_flow += cur_flow;
		}
		return tot_flow;
	}

	// NOTE: Assume that we already call max_flow() before use this function
	template <typename Network>
	std::vector<bool> min_cut(const Network& g) {
		std::vector<bool> res(g.V);
		for (int i = 0; i < g.V; ++ i) {
			res[i] = (depth[i] != - 1);
		}
		return res;
		// Return a side of each node
		// res[u] = 0 if u on the left partition (can reachable from source)
		// res[u] = 1 on the contrast case
		// For edge e, if res[e.from] != res[e.to] --> e is the cut
	}
};

// Don't know why it works, just copied from: https://loj.ac/s/1488185
// flow_t = double may be give wrong flow (never check it !!!)
// Time: O(V ^ 2 * sqrt(E))
// Just use when dinic gives TLE
template <typename flow_t> struct hlpp_max_flow {
	static constexpr flow_t INF_FLOW = std::numeric_limits<flow_t>::max() / 2;
	std::vector<int> nxt;
	std::vector<int> lst;
	std::vector<flow_t> excess;
	std::vector<int> arc;
	std::vector<int> nxt_gap;
	std::vector<int> prv_gap;
	std::vector<int> depth;
	std::vector<int> q;

	int max_depth;
	int max_gap;
	int cnt;

	inline void init(const int& N) {
		nxt.assign(N, 0);
		lst.assign(N, 0);
		excess.assign(N, 0);
		arc.assign(N, 0);
		nxt_gap.assign(N << 1, 0);
		prv_gap.assign(N << 1, 0);
		depth.assign(N, 0);
		q.assign(N, 0);
		max_depth = max_gap = cnt = 0;
	}

	hlpp_max_flow() {}
	hlpp_max_flow(const int& N) { init(N); }

	inline void push_lst(const int& h, const int& u) {
		nxt[u] = lst[h]; lst[h] = u;
	}
	template <typename Network>
	inline void update_depth(const Network& g, const int& v, int nh) {
		if (depth[v] != g.V) {
			nxt_gap[prv_gap[v]] = nxt_gap[v];
			prv_gap[nxt_gap[v]] = prv_gap[v];
		}
		depth[v] = nh;
		if (nh == g.V) return;
		max_gap = std::max(max_gap, nh);
		if (excess[v] > g.eps) {
			max_depth = std::max(max_depth, nh);
			push_lst(nh, v);
		}
		nh += g.V;
		nxt_gap[v] = nxt_gap[nh];
		prv_gap[v] = nh;
		nxt_gap[nh] = v;
		prv_gap[nxt_gap[v]] = v;
	}
	template <typename Network>
	inline void global_relabel(const Network& g, const int& s, const int& t) {
		cnt = 0;
		std::fill(depth.begin(), depth.end(), g.V);
		std::fill(lst.begin(), lst.end(), - 1);
		std::iota(nxt_gap.begin(), nxt_gap.end(), 0);
		std::iota(prv_gap.begin(), prv_gap.end(), 0);
		depth[t] = 0; q[0] = t;
		for (int i = 0, j = 1; i != j; ++ i) {
			const auto& u = q[i];
			for (const auto& id : g.adj[u]) {
				const auto& e = g.edges[id];
				const auto& back = g.edges[id ^ 1];
				if (depth[e.to] == g.V && back.cap - back.flow > g.eps) {
					update_depth(g, e.to, depth[u] + 1);
					q[j ++] = e.to;
				}
			}
			max_depth = max_gap = depth[u];
		}
	}
	template <typename Network>
	inline void push(Network& g, const int& u, const int& id) {
		auto& e = g.edges[id];
		flow_t x = std::min<flow_t>(excess[u], e.cap - e.flow);
		if (x > g.eps) {
			if (excess[e.to] <= g.eps && excess[e.to] >= 0) push_lst(depth[e.to], e.to);
			e.flow += x;
			g.edges[id ^ 1].flow -= x;
			excess[u] -= x;
			excess[e.to] += x;
		}
	}
	template <typename Network>
	inline void discharge(Network& g, const int& u) {
		int nh = g.V;
		for (int i = arc[u]; i != int(g.adj[u].size()); ++ i) {
			const auto& e = g.edges[g.adj[u][i]];
			if (e.cap - e.flow > g.eps) {
				if (depth[u] == depth[e.to] + 1) {
					push(g, u, g.adj[u][i]);
					if (excess[u] <= g.eps) {
						arc[u] = i; return;
					}
				} else {
					nh = std::min(nh, depth[e.to] + 1);
				}
			}
		}

		for (int i = 0; i != arc[u]; ++ i) {
			const auto& e = g.edges[g.adj[u][i]];
			if (e.cap - e.flow > g.eps) {
				if (depth[u] == depth[e.to] + 1) {
					push(g, u, g.adj[u][i]);
					if (excess[u] <= 0) {
						arc[u] = i; return;
					}
				} else {
					nh = std::min(nh, depth[e.to] + 1);
				}
			}
		}
		++ cnt;
		if (nxt_gap[nxt_gap[depth[u] + g.V]] != depth[u] + g.V) {
			update_depth(g, u, nh);
		} else {
			const int& prv_depth = depth[u];
			for (int h = prv_depth; h <= max_gap; ++ h) {
				for (int i = nxt_gap[h + g.V]; i < g.V; i = nxt_gap[i]) {
					depth[i] = g.V;
				}
				nxt_gap[h + g.V] = prv_gap[h + g.V] = h + g.V;
			}
			max_gap = prv_depth - 1;
		}
	}
	template <typename Network>
	flow_t max_flow(Network& g, const int& s, const int& t) {
		excess[s] = INF_FLOW; excess[t] = - INF_FLOW;
		global_relabel(g, s, t);
		for (const auto& id : g.adj[s]) {
			push(g, s, id);
		}
		for (; max_depth >= 0; -- max_depth) {
			while (lst[max_depth] != - 1) {
				int u = lst[max_depth];
				lst[max_depth] = nxt[u];
				if (depth[u] == max_depth) {
					discharge(g, u);
					if (cnt > (g.V << 2)) {
						global_relabel(g, s, t);
					}
				}
			}
		}
		return excess[t] + INF_FLOW;
	}

	template <typename Network>
	std::vector<bool> min_cut(const Network& g) {
		std::vector<bool> res(g.V);
		for (int i = 0; i < g.V; ++ i) {
			res[i] = (depth[i] != g.V);
		}
		return res;
		// Return a side of each node
		// res[u] = 0 if u on the left partition (can reachable from source)
		// res[u] = 1 on the contrast case
		// For edge e, if res[e.from] != res[e.to] --> e is the cut
	}
};