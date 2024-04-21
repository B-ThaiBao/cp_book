/**
 * rmq_ancestor: lca is based on rmq + euler tour
 * 
 * Requires: fast rmq only in O(N) / O(1)
 * 
 * Usage:
 *   * Have fast rmq (using range_min_query template)
 *   * Init: init(N)
 *   * Build euler tour: dfs(g, u)
 *   * Build_rmq: build_rmq()
 *   * Find lca (u, v): find_lca(u, v) only in O(1) !!!
 * 
 * TODO: Please build euler tour before build rmq
**/
struct rmq_ancestor {
	std::vector<int> depth;
	std::vector<int> loc;
	std::vector<int> euler;
	range_min_query<int> rmq; // Fast ???

	void init(const int& N) {
		depth.assign(N, - 1);
		loc.assign(N, - 1);
		euler.reserve(N << 1);
	}
	rmq_ancestor() {}
	rmq_ancestor(const int& N) { this -> init(N); }

	void build_rmq() {
		std::vector<int> depths; depths.reserve(euler.size());
		for (const auto& u : euler) {
			depths.push_back(u < 0 ? u : depth[u]);
		}
		rmq.build(std::move(depths));
	}
	template <typename G> void build(const G& g) {
		if (depth.empty()) init(g.V);
		dfs(g); build_rmq();
	}

	template <typename F> void do_dfs_from(const F& g, const int& u) {
		loc[u] = int(euler.size());
		euler.push_back(u);
		for (const auto& id : g.adj[u]) {
			if (g.is_ignore(id)) continue;
			const auto& e = g.edges[id];
			auto nxt = e.from ^ e.to ^ u;
			if (depth[nxt] != - 1) continue;
			depth[nxt] = depth[u] + 1;
			do_dfs_from(g, nxt);
			euler.push_back(u);
		}
	}
	template <typename F> void dfs(const F& g, const int& u) {
		if (depth.empty()) init(g.V);
		depth[u] = 0;
		do_dfs_from(g, u);
		euler.push_back(- 1);
	}
	template <typename F> void dfs(const F& g) {
		if (depth.empty()) init(g.V);
		for (int u = 0; u < g.V; u ++) {
			if (depth[u] != - 1) continue;
			depth[u] = 0;
			do_dfs_from(g, u);
			euler.push_back(- 1);
		}
	}

	inline int find_lca(int u, int v) const {
		u = loc[u]; v = loc[v];
		if (u > v) std::swap(u, v);
		return euler[rmq.range_query(u, v).first];
	}

	// Note: this is the LCA of any two nodes out of three when the third node is the root.
	// It is also the node with the minimum sum of distances to all three nodes
	// (the centroid of the three nodes).
	inline int rooted_lca(const int& a, const int& b, const int& c) const {
		// Return the deepest node among lca(a, b), lca(b, c), and lca(c, a).
		return find_lca(a, b) ^ find_lca(b, c) ^ find_lca(c, a);
	}
};

/**
 * HEAVY LIGHT DECOMPOSITION (HLD) !!!! ????
 * 
 * 
 * Source (for understand HLD means ???):
 * 
 * Codechef: https://www.youtube.com/watch?v=LK-0QNShaPM
 * USACO: https://usaco.guide/plat/hld?lang=cpp
 * CP-Algorithm: https://cp-algorithms.com/graph/hld.html
 * USACO: https://usaco.guide/CPH.pdf#page=174
 * CF: https://codeforces.com/blog/entry/12239
 * 
 * 
 * 
 * Source (for implementation):
 * 
 * https://codeforces.com/blog/entry/22072
 * https://codeforces.com/blog/entry/53170 (Subtree queries)
 * https://github.com/bqi343/cp-notebook/blob/master/Implementations/content/graphs%20(12)/Trees%20(10)/HLD%20(10.3).h
 * https://github.com/nealwu/competitive-programming/blob/master/heavy_light/subtree_heavy_light.cc (Neal)
 * https://github.com/the-tourist/algo/blame/master/graph/hld_forest.cpp (Tourist)
 * https://github.com/ecnerwala/cp-book/blob/master/src/level_ancestor.hpp (Ecnerwala)
 * https://codeforces.com/contest/1844/submission/213393963 (Jiangly)
 * https://github.com/ShahjalalShohag/code-library/blob/master/Data%20Structures/HLD.cpp
 * https://codeforces.com/contest/1843/submission/210740643 (bashsort)
 * https://nyaannyaan.github.io/library/tree/heavy-light-decomposition.hpp
 * https://codeforces.com/contest/1827/submission/205915582 (Snow-flower)
 * 
 * For more details, dfs_ancestor make tree fit into an array for apply data structure
 * on the node or on the path. For this reason, this code provides you some infor in the
 * tree:
 * 
 * depth[v]: depth of node v
 * par[v]: parent of node v
 * sz[v]: size of subtree v
 * root[v]: root or heavy_par of chain that contains v
 * idx[v]: the idx of v in the array after dfs
 * order[v]: the real node at the idx v of the array
 * end[v]: the time end of dfs all subtree v
 * pe[v]: the id of edge up to par[v]
 * 
 * Furthermore, this implementation can also provides some method that help queries and
 * update much more easy (See more details in the code below) ???? !!!!!!!
 * 
 * Usage:
 *   Init with N nodes: dfs_ancestor lca(N);
 * 
 *   Dfs in the components with root node is u: lca.dfs(adj, u);
 * 
 *   Iterate path from u -> v in O(log_2(N)): no need to care about the order
 *     find_lca(u, v, is_with_lca, f) with f(L, R, is_from_u_to_v);
 * 
 *   Iterate path from u -> v in O(log_2(N)): the order from u -> v is vital
 *     Similar but apply on 2 path: u -> lca(u, v) and lca(u, v) -> v;
 * 
 *   Iterate path from subtree with root u in O(1):
 *     for_subtree(u, is_with_root, f) with f(L, R)
 * 
 *   Iterate u -> kth ancestor of u in O(log_2(K)): direction is always go up
 *     find_ancestor(u, k, is_with_kth_node, f) with f(L, R)
 *   Beware that if k > depth[u] return - 1 but all node will be updated
 * 
 *   Iterate kth node from u -> v in in O(log_2(K)):
 *     find_node_path(u, v, k, f) with f(L, R, is_from_u_to_v);
 *   Beware that if k > dist(u, v) return - 1 but all node will be updated
 * 
**/
struct dfs_ancestor {
	std::vector<int> depth;
	std::vector<int> par;
	std::vector<int> pe;
	std::vector<int> sz;
	std::vector<int> root;
	std::vector<int> idx;
	std::vector<int> order;
	std::vector<int> end;
	std::vector<int> ch; // Yeah, this is heavy child
	                     // TODO: Try not swap in dfs_chain

	int time = 0;
	void init(const int& V){
		time = 0;
		depth.assign(V, - 1);
		par.assign(V, - 1);
		pe.assign(V, - 1);
		sz.assign(V, 1);
		root.assign(V, - 1);
		idx.assign(V, - 1);
		order.assign(V, - 1);
		end.assign(V, - 1);
		ch.assign(V, - 1);
	}
	dfs_ancestor() {}
	dfs_ancestor(const int& N) { init(N); }

	template <typename G> void do_dfs_from(const G& g, const int& cur) {
		for (const auto& id : g.adj[cur]) {
			if (g.is_ignore(id)) continue;
			const auto& e = g.edges[id];
			int nxt = e.from ^ e.to ^ cur;
			if (nxt == par[cur]) continue;
			depth[nxt] = depth[cur] + 1;
			par[nxt] = cur;
			pe[nxt] = id;
			do_dfs_from(g, nxt);
			sz[cur] += sz[nxt];
			if (ch[cur] == - 1 || sz[nxt] > sz[ch[cur]]) {
				ch[cur] = nxt; // We have a heavy child right now
			}
		}
	}
	template <typename G> void dfs_from(const G& g, const int& cur) {
		order[idx[cur] = time ++] = cur;
		if (ch[cur] != - 1) {
			root[ch[cur]] = root[cur];
			dfs_from(g, ch[cur]);
		}
		for (const auto& id : g.adj[cur]) {
			if (g.is_ignore(id)) continue;
			const auto& e = g.edges[id];
			int nxt = e.from ^ e.to ^ cur;
			if (nxt == par[cur] || nxt == ch[cur]) continue;
			root[nxt] = nxt;
			dfs_from(g, nxt);
		}
		end[cur] = time;
	}
	template <typename G> void dfs(const G& g, const int& u) {
		depth[u] = 0;
		par[u] = pe[u] = - 1;
		do_dfs_from(g, u);
		root[u] = u;
		dfs_from(g, u);
	}
	template <typename G> void dfs(const G& g) {
		for (int u = 0; u < g.V; u ++) {
			if (par[u] != - 1) continue;
			depth[u] = 0;
			par[u] = pe[u] = - 1;
			do_dfs_from(g, u);
			root[u] = u;
			dfs_from(g, u);
		}
	}

	template <typename F>
	int find_lca(int u, int v, const bool& with_lca, F&& f) const {
		bool dir = true;
		while (root[u] != root[v]) {
			if (depth[root[u]] > depth[root[v]]){
				std::swap(u, v), dir = !dir;
			}
			auto& vroot = root[v];
			f(idx[vroot], idx[v], dir);
			v = par[vroot];
		}
		if (depth[u] > depth[v]) {
			std::swap(u, v), dir = !dir;
		}
		// u is ancestor of v
		int idx_u = idx[u] + !with_lca;
		if (idx_u <= idx[v]) f(idx_u, idx[v], dir);
		return u;
	}
	int find_lca(int u, int v) const {
		while (root[u] != root[v]) {
			if (depth[root[u]] > depth[root[v]]) {
				std :: swap(u, v);
			}
			v = par[root[v]];
		}
		return depth[u] < depth[v] ? u : v;
	}
	// Note: this is the LCA of any two nodes out of three when the third node is the root.
	// It is also the node with the minimum sum of distances to all three nodes
	// (the centroid of the three nodes).
	inline int rooted_lca(const int& a, const int& b, const int& c) const {
		// Return the deepest node among lca(a, b), lca(b, c), and lca(c, a).
		return find_lca(a, b) ^ find_lca(b, c) ^ find_lca(c, a);
	}

	template <typename F>
	void subtree(const int& u, const bool& rooted, F&& f) const {
		int idx_u = idx[u] + !rooted;
		if (idx_u <= end[u] - 1) f(idx_u, end[u] - 1);
	}

	template <typename F>
	int find_ancestor(int u, int k, const bool& with_high, F&& f) const {
		while (u != - 1 && k > - 1) {
			const auto& root_u = root[u];
			if (k > depth[u] - depth[root_u]) {
				f(idx[root_u], idx[u]);
				k -= (depth[u] - depth[root_u] + 1);
				u = par[root_u];
			}
			else {
				int idx_u = idx[u];
				int a = idx_u - k + !with_high;
				if (a <= idx_u) f(a, idx_u);
				return order[a - !with_high];
			}
		}
		return - 1;
	}
	int find_ancestor(int u, int k) const {
		while (u != - 1 && k > - 1) {
			const auto& root_u = root[u];
			if (k > depth[u] - depth[root_u]) {
				k -= (depth[u] - depth[root_u] + 1);
				u = par[root_u];
			}
			else {
				return order[idx[u] - k];
			}
		}
		return - 1;
	}

	template <typename F>
	int find_node_path(int u, int v, int k, F&& f) const {
		int c = find_lca(u, v);
		if (k <= depth[u] - depth[c]) {
			return find_ancestor(u, k, true, [&](const int& L, const int& R) {
				f(L, R, false);
			});
		}
		else {
			k = depth[u] + depth[v] - (depth[c] << 1) - k;
			if (k <= 0) {
				find_lca(u, v, true, f); return (k == 0) ? v : - 1;
			}
			else {
				find_ancestor(u, depth[u] - depth[c], true, [&](const int& L, const int& R) {
					f(L, R, false);
				});
				v = find_ancestor(v, k);
				find_ancestor(v, depth[v] - depth[c] - 1, true, [&](const int& L, const int& R) {
					f(L, R, true);
				});
				return v;
			}
		}
	}
	int find_node_path(int u, int v, int k) const {
		int c = find_lca(u, v);
		if (k <= depth[u] - depth[c]) {
			return find_ancestor(u, k);
		}
		else {
			k = depth[u] + depth[v] - (depth[c] << 1) - k;
			if (k < 0) return - 1;
			else {
				return find_ancestor(v, k);
			}
		}
	}
};
