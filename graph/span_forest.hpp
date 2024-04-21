struct dfs_span_forest {
	// TODO: Try to make graph -> span forest by using 1 DFS-traversal
	// Can easily detect back_edge and non-vital node
	std::vector<int> par; // Not need but use too much in helper function
	std::vector<int> pe;
	std::vector<int> order;
	std::vector<int> idx;
	std::vector<int> end;
	std::vector<int> sz;
	std::vector<int> root;
	std::vector<int> depth;
	std::vector<int> min_depth;

	void init(const int& N) {
		par.assign(N, - 1);
		pe.assign(N, - 1);
		order.reserve(N); // Avoid relocation
		idx.assign(N, - 1);
		end.assign(N, - 1);
		sz.assign(N, 0);
		root.assign(N, - 1);
		depth.assign(N, - 1);
		min_depth.assign(N, - 1);
	}
	void clear() {
		par.clear();
		pe.clear(); order.clear();
		idx.clear(); end.clear();
		sz.clear(); root.clear();
		depth.clear();
	}

	dfs_span_forest() {}
	dfs_span_forest(const int& N) { this -> init(N); }

	template <typename G> void do_dfs(const G& g, const int& v) {
		idx[v] = (int) order.size();
		order.push_back(v);
		sz[v] = 1;
		min_depth[v] = depth[v];
		for (auto& id : g.adj[v]) {
			if (id == pe[v] || g.is_ignore(id)) continue;
			const auto& e = g.edges[id];
			int to = e.from ^ e.to ^ v;
			if (root[to] == root[v]) {
				// TODO: Link v with the highest node to tree (using 1 back_edge)
				min_depth[v] = std::min(min_depth[v], depth[to]);
				continue;
			}
			depth[to] = depth[v] + 1;
			pe[to] = id; // This is upper edge in span tree
			par[to] = v;
			root[to] = (root[v] != - 1 ? root[v] : to);
			do_dfs(g, to);
			sz[v] += sz[to];
			min_depth[v] = std::min(min_depth[v], min_depth[to]);
		}
		end[v] = int(order.size()) - 1;
	}
	template <typename G> void do_dfs_from(const G& g, const int& v) {
		depth[v] = 0; root[v] = v; pe[v] = - 1; par[v] = - 1;
		do_dfs(g, v);
	}
	template <typename G> void dfs(const G& g, const int& v, const bool& clear_order = false) {
		if (pe.empty()) init(g.V);
		else if (clear_order) order.clear(); // Just keep order in this comp
		do_dfs_from(g, v);
	}
	template <typename G> void dfs(const G& g) {
		if (pe.empty()) init(g.V);
		for (int v = 0; v < g.V; v ++) {
			if (depth[v] == - 1) do_dfs_from(g, v);
		}
		// assert(int(order.size()) == g.V);
	}
	template <typename G> void dfs(const G& g, std::vector<int>& s) {
		if (pe.empty()) init(g.V);
		for (const auto& v : s) {
			if (depth[v] == - 1) do_dfs_from(g, v);
		}
	}

	template <typename G>
	bool is_back_edge(const G& g, const int& id, const int& u) const {
		// TODO: wheter id is backedge from u to ancestor (not to par[u])
		// This is useful to find the cycle
		int v = g(u, id);
		return id != pe[u] && is_ancestor(v, u);
	}
	template <typename G>
	bool is_span_edge(const G& g, const int& id, const int& u) const {
		return id != pe[u];
	}
	bool is_ancestor(const int& u, const int& v) const {
		// TODO: Check if u is ancestor of v
		return idx[u] <= idx[v] && end[u] <= end[v];
	}

	template <typename G> std::vector<bool> find_cutpoint(const G& g) const {
		std::vector<bool> is_cut(g.V, false);
		for (int i = 0; i < g.V; i ++) {
			if (par[i] != - 1 && min_depth[i] >= depth[par[i]]) {
				is_cut[par[i]] = true;
			}
		}
		std::vector<int> ch(g.V, 0);
		for (int i = 0; i < g.V; i ++) {
			if (par[i] != - 1) ++ ch[par[i]];
		}
		for (int i = 0; i < g.V; i ++) {
			if (par[i] == - 1 && ch[i] < 2) is_cut[i] = false;
		}
		return is_cut;
	}
	template <typename G> std::vector<bool> find_bridge(const G& g) const {
		std::vector<bool> is_bri(g.edges.size(), false);
		for (int i = 0; i < g.V; i ++) {
			if (par[i] != - 1 && min_depth[i] == depth[i]) {
				is_bri[pe[i]] = true;
			}
		}
		return is_bri;
	}
	template <typename G, typename F> void find_cycle(const G& g, const F& f) const {
		// TODO: Find all cycle based on the back edge
		std::vector<int> es; es.reserve(g.edges.size());
		for (const auto& u : order) {
			for (const auto& id : g.adj[u]) {
				if (id == pe[u]) continue;
				int v = g(u, id);
				if (!is_ancestor(v, u)) continue;
				int us = u;
				while (us != v) es.push_back(pe[us]), us = par[us];
				std::reverse(es.begin(), es.end());
				es.push_back(id);
				if (!f(v, es)) return;
				es.clear();
			}
		}
	}

	// NOTE: This returns id (on each vertex) of biconnected components
	// In each component, we doesn't allow any brigde edges
	// ==> All edges in biconnect_compe are brigdes in the original graph 
	template <typename G> std::vector<int> bicon_compe(const G& g, int& cnt) const {
		std::vector<int> vertex_comp(g.V);
		cnt = 0;
		for (const auto& i : order) {
			if (par[i] == - 1 || min_depth[i] == depth[i]) {
				vertex_comp[i] = cnt ++;
			}
			else {
				vertex_comp[i] = vertex_comp[par[i]];
			}
		}
		return vertex_comp;
	}

	// NOTE: This returns id (on each edges) of biconnected components
	// In each component, we doesn't allow any cutpoint
	// If each component shrinked into a node and attach with shared vertices is cutpoint
	// More details: https://tanujkhattar.wordpress.com/2016/01/10/the-bridge-tree-of-a-graph/
	template <typename G> std::vector<int> bicon_compv(const G& g, int& cnt) const {
		std::vector<int> vertex_comp(g.V);
		cnt = 0;
		for (const auto& i : order) {
			if (par[i] == - 1){
				vertex_comp[i] = - 1;
				continue;
			}
			if (min_depth[i] >= depth[par[i]]) {
				vertex_comp[i] = cnt ++;
			}
			else {
				vertex_comp[i] = vertex_comp[par[i]];
			}
		}
		std::vector<int> res(g.edges.size(), - 1);
		for (int i = 0; i < int(g.edges.size()); i ++) {
			if (g.is_ignore(i)) continue;
			int x = g.edges[i].from;
			int y = g.edges[i].to;
			// TODO: choose deeper node
			int z = (depth[x] > depth[y] ? x : y);
			res[i] = vertex_comp[z];
		}
		return res;
	}
};

struct bfs_span_forest {
	std::vector<int> par;
	std::vector<int> pe;
	std::vector<int> order;
	std::vector<int> idx;
	std::vector<int> root;
	std::vector<int> depth;

	void init(const int& N) {
		par.assign(N, - 1); pe.assign(N, - 1);
		order.reserve(N); idx.assign(N, - 1);
		root.assign(N, - 1); depth.assign(N, - 1);
	}
	void clear() {
		par.clear(); pe.clear(); order.clear();
		root.clear(); depth.clear();
	}

	bfs_span_forest() {}
	bfs_span_forest(const int& N) { this -> init(N); }

	template <typename G>
	void bfs(const G& g, const int& s, const bool& clear_order = false) {
		if (clear_order) order.clear();
		if (pe.empty()) init(g.V);
		std::vector<int> q; q.reserve(g.V);
		depth[s] = 0;
		root[s] = s;
		par[s] = pe[s] = - 1;
		q.push_back(s);
		order.push_back(s);
		int beg = 0;
		while (beg < int(q.size())) {
			auto u = q[beg ++];
			idx[u] = int(order.size());
			order.push_back(u);
			for (const auto& id : g.adj[u]) {
				if (g.is_ignore(id)) continue;
				const auto& e = g.edges[id];
				int v = e.from ^ e.to ^ u;
				if (depth[v] != - 1) continue; // Already
				depth[v] = depth[u] + 1;
				par[v] = u;
				pe[v] = id;
				root[v] = (root[u] != - 1) ? root[u] : v;
				q.push_back(v);
			}
		}
	}
	template <typename G> void bfs(const G& g) {
		for (int i = 0; i < g.V; i ++) {
			if (depth[i] == - 1) bfs(g, i);
		}
	}
	template <typename G>
	void bfs(const G& g, const std::vector<int>& roots, const bool& clear_order = false) {
		if (clear_order) order.clear();
		if (pe.empty()) init(g.V);
		std::vector<int> q; q.reserve(g.V); // Instead of std::queue<int>
		for (const auto& s : roots) {
			if (depth[s] != - 1) continue;
			depth[s] = 0;
			root[s] = s;
			par[s] = pe[s] = - 1;
			q.push_back(s);
			order.push_back(s);
		}
		int beg = 0;
		while (beg < int(q.size())) {
			auto u = q[beg ++];
			idx[u] = int(order.size());
			order.push_back(u);
			for (const auto& id : g.adj[u]) {
				if (g.is_ignore(id)) continue;
				const auto& e = g.edges[id];
				int v = e.from ^ e.to ^ u;
				if (depth[v] != - 1) continue; // Already
				depth[v] = depth[u] + 1;
				par[v] = u;
				pe[v] = id;
				root[v] = (root[u] != - 1) ? root[u] : v;
				q.push_back(v);
			}
		}
	}

	template <typename G> bool is_span_edge(const G& g, const int& id, const int& u) const {
		return id == pe[u];
	}
};