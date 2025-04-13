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

// Don't know why it works, just copied from:
//   https://github.com/koosaga/olympiad/blob/master/Library/codes/combinatorial_optimization/flow.cpp
// Tested: https://loj.ac/s/2083847
// flow_t = double may be give wrong flow (never check it !!!)
// Time: O(V ^ 2 * sqrt(E))
// Just use when dinic gives TLE
template <typename flow_t> struct hlpp_max_flow {
	std::vector<int> lst;
	std::vector<flow_t> excess;
	std::vector<int> arc;
	std::vector<int> nxt_gap;
	std::vector<int> prv_gap;
	std::vector<int> depth;
	std::vector<int> active;
	std::vector<int> q;

	int max_depth;
	int cnt;
	int cut;

	inline void init(const int& N) {
		active.assign(N << 1, 0);
		lst.assign(N, 0);
		excess.assign(N, 0);
		arc.assign(N, 0);
		nxt_gap.assign(3 * N, 0);
		prv_gap.assign(3 * N, 0);
		depth.assign(N, N << 1);
		q.assign(N, 0);
		max_depth = cut = cnt = 0;
	}

	hlpp_max_flow() {}
	hlpp_max_flow(const int& N) { this->init(N); }

	template <typename Network>
	inline void update_depth(const Network& g, int v, int nh) {
#ifdef _GLIBCXX_DEBUG
		assert(0 <= v && v < g.V);
#endif
		prv_gap[nxt_gap[prv_gap[v]] = nxt_gap[v]] = prv_gap[v];
		depth[v] = nh;
		if (excess[v] > g.eps) {
			lst[v] = active[nh]; active[nh] = v;
			max_depth = std::max(max_depth, nh);
		}
		if (nh < g.V) cut = std::max(cut, nh + 1);
		nxt_gap[v] = nxt_gap[prv_gap[v] = nh += g.V];
		prv_gap[nxt_gap[nxt_gap[nh] = v]] = v;
	}
	template <typename Network>
	inline void global_relabel(const Network& g, const int& s, const int& t) {
#ifdef _GLIBCXX_DEBUG
		assert(0 <= s && s < g.V);
		assert(0 <= t && t < g.V);
#endif
		std::fill(depth.begin(), depth.end(), g.V << 1);
		std::fill(active.begin(), active.end(), -1);
		std::iota(nxt_gap.begin(), nxt_gap.end(), 0);
		std::iota(prv_gap.begin(), prv_gap.end(), 0);
		max_depth = cnt = cut = 0;
		depth[t] = 0; depth[s] = g.V;

		{ // Starting with t
			q[0] = t;
			for (int i = 0, j = 1; i < j; ++i) {
				const auto& v = q[i];
				for (const auto& id : g.adj[v]) {
					const auto& e = g.edges[id];
					const auto& back_e = g.edges[id ^ 1];
					if (depth[e.to] == (g.V << 1) && back_e.cap - back_e.flow > g.eps) {
						q[j++] = e.to;
						this->update_depth(g, e.to, depth[v] + 1);
					}
				}
			}
		}

		{ // Starting with s
			q[0] = s;
			for (int i = 0, j = 1; i < j; ++i) {
				const auto& v = q[i];
				for (const auto& id : g.adj[v]) {
					const auto& e = g.edges[id];
					const auto& back_e = g.edges[id ^ 1];
					if (depth[e.to] == (g.V << 1) && back_e.cap - back_e.flow > g.eps) {
						q[j++] = e.to;
						this->update_depth(g, e.to, depth[v] + 1);
					}
				}
			}
		}
	}
	template <typename Network>
	inline void push(Network& g, const int& u, const int& id, const bool& z) {
#ifdef _GLIBCXX_DEBUG
		assert(0 <= u && u < g.V);
#endif
		auto& e = g.edges[id];
		flow_t x = std::min<flow_t>(excess[u], e.cap - e.flow);
		if (x > g.eps) {
			if (z && 0 <= excess[e.to] && excess[e.to] <= g.eps) {
				lst[e.to] = active[depth[e.to]];
				active[depth[e.to]] = e.to;
			}
			e.flow += x;
			g.edges[id ^ 1].flow -= x;
			excess[u] -= x;
			excess[e.to] += x;
		}
	}
	template <typename Network>
	inline void discharge(Network& g, const int& v) {
#ifdef _GLIBCXX_DEBUG
		assert(0 <= v && v < g.V);
#endif
		int h = g.V << 1, k = depth[v];
 
		for(int j = 0; j < int(g.adj[v].size()); j++){
			const auto& e = g.edges[g.adj[v][arc[v]]];
			if (e.cap - e.flow > g.eps) {
				if (k == depth[e.to] + 1) {
					this->push(g, v, g.adj[v][arc[v]], true);
					if (excess[v] <= g.eps) return;
				} else {
					h = std::min(h, depth[e.to] + 1);
				}
			}
			if (++arc[v] >= int(g.adj[v].size())) arc[v] = 0;
		}
		if (k < g.V && nxt_gap[k + g.V] == prv_gap[k + g.V]) {
			for (int j = k; j < cut; j++) {
				while (nxt_gap[j + g.V] < g.V) {
					this->update_depth(g, nxt_gap[j + g.V], g.V);
				}
			}
			cut = k;
		} else {
			this->update_depth(g, v, h), ++cnt;
		}
	}
	template <typename Network>
	flow_t max_flow(Network& g, const int& s, const int& t) {
#ifdef _GLIBCXX_DEBUG
		assert(0 <= s && s < g.V);
		assert(0 <= t && t < g.V);
#endif
		for (const auto& id : g.adj[s]) {
			excess[s] = g.edges[id].cap - g.edges[id].flow;
			this->push(g, s, id, false);
		}
		global_relabel(g, s, t);
		for (; max_depth > 0; --max_depth) {
			while (active[max_depth] != -1) {
				int u = active[max_depth];
				active[max_depth] = lst[u];
				if (u != s && depth[u] == max_depth) {
					discharge(g, u);
					if (cnt > (g.V << 2)) {
						this->global_relabel(g, s, t);
					}
				}
			}
		}
		return excess[t];
	}

	template <typename Network>
	std::vector<bool> min_cut(const Network& g) {
		std::vector<bool> res(g.V);
		for (int i = 0; i < g.V; ++ i) {
			res[i] = (depth[i] < g.V);
		}
		return res;
		// Return a side of each node
		// res[u] = 0 if u on the left partition (can reachable from source)
		// res[u] = 1 on the contrast case
		// For edge e, if res[e.from] != res[e.to] --> e is the cut
	}
};
