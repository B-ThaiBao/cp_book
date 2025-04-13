/**
 * CENTROID TREE !!
 *
 * Properties:
 *   * Height of the centroid tree is O(log(N))
 *   * A vertex belongs to the component of all its ancestors.
 *   * Path A -> B: A -> lca(A, B) and lca(A, B) -> B, which lca(A, B)
 *     is the lowest common ancestor of A and B in the centroid tree.
 *   * Each one of the N * N paths of the original tree is the concatenation
 *     of two paths in a set of O(N * log(N)) paths from a node to all its
 *     ancestors in the centroid decomposition.
 *     --> The idea is that instead of choosing all possible endpoints of the
 *         path in the original tree, we’ll choose all the lowest common ancestors
 *         of two nodes in the centroid decomposition.
 *
 * Summary: Use this when you must
 *   * consider all possible path in tree --> consider all paths go through centroids
 *   * consider lots of queries with different nodes to all others
 *        --> consider that node with all its ancestors
 *
 * NOTE: The centroids has been erased from the original tree is ct.was[u] == ct.iter
 * Usage:
 *   * centroid_tree ct(N);
 *   * Dfs to get size after erase some nodes (was[u] == iter): ct.dfs(g, u)
 *   * Build_tree at node that has sz[u] = -1: ct.build_tree(g, u) if ct.sz[u] == -1
 *   * After build_tree() from all componenents, the result is forest with par[u] and pe[u]
 *   * If want to do centroid decomposition again, set ct.iter++ --> Dfs again
 *   * More advanced: centroid_tree is built by D&C algo, so we can avoid consider all
 *     paths (N * N) by consider all paths go through centroids (N * log(N))
 *     --> Can apply on centroid before erase it with O(A) --> Time O(N * log(N) + A * log(N))
 *     --> Call apply(cur, rt) where cur is the centroid and rt is the root of original subtree
**/
struct centroid_tree {
	std::vector<int> sz;
	std::vector<int> pe;
	std::vector<int> par;
	std::vector<int> was;
	int iter = 0;

	template <typename Graph>
	inline int find_centroid(const Graph& g, const int& s) const {
		int u = s, cur_pe = -1;
		while (true) {
			for (const auto& id : g.adj[u]) {
				if (id == cur_pe || g.is_ignore(id)) continue;
				int v = g(u, id);
				if (was[v] == iter) continue;
				if (sz[v] > (sz[s] >> 1)) {
					u = v, cur_pe = id;
					goto nxt;
				}
			}
			break;
		nxt:;
		}
		return u;
	}

	template <typename Graph>
	inline void dfs(const Graph& g, const int& cur, const int& p = -1) {
		sz[cur] = 1;
		for (const auto& id : g.adj[cur]) {
			if (g.is_ignore(id)) continue;
			int nxt = g(cur, id);
			if (nxt == p || was[nxt] == iter) continue;
			dfs(g, nxt, cur);
			sz[cur] += sz[nxt];
		}
	}

	template <typename Graph, typename F>
	inline void do_dfs(const Graph& g, int cur, const int& p, const int& prv, const F& apply) {
		dfs(g, cur, p);
		auto rt = cur; cur = find_centroid(g, cur);
		par[cur] = p;
		pe[cur] = prv;
		// TODO: Apply on centroid before delete it from the original tree
		apply(cur, rt);
		// Now, we delete it
		was[cur] = iter;
		for (const auto& id : g.adj[cur]) {
			if (g.is_ignore(id)) continue;
			int nxt = g(cur, id);
			if (was[nxt] != iter) do_dfs(g, nxt, cur, id, apply);
		}
	}

	template <typename Graph, typename F>
	inline void build_tree(const Graph& g, const int& r, const F& apply) {
		if (was[r] == iter) return;
		do_dfs(g, r, -1, -1, apply);
	}
	template <typename Graph>
	inline void build_tree(const Graph& g, const int& r) {
		if (was[r] == iter) return;
		do_dfs(g, r, -1, -1, [](const auto&, const auto&) {});
	}

	inline void assign(const int& N, int v = -1) {
		sz.assign(N, v);
		pe.assign(N, v);
		par.assign(N, v);
		was.assign(N, v);
		iter = 0;
	}

	centroid_tree() = default;
	centroid_tree(const int& N, int v = -1) { this->assign(N, v); }
};
