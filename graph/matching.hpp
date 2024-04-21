/**
 * DFS_MATCHING
 * 
 * Assume that we have 2 partitions size N (0 ... N - 1) and size M (0 ... M - 1)
 * (can have same id and some node can ignore). Can find max_match from L to R.
 * 
 * Tested:
 *   * https://www.spoj.com/problems/MATCHING/
 *   * https://www.spoj.com/problems/PT07X/
 *   * https://codeforces.com/contest/722/submission/242448054
 * 
 * Related:
 *   * num min_ver_cover = max_match
 *   * max_independent_set = all_nodes - min_ver_cover
 *     (Nodes in this set is the complement with cover vertex)
 * 
 * Usage:
 *   * Init N, M partitions: init(N, M);
 * 
 *   * Set node u is ignore in left or right partition:
 *         mate[0][i] = - 2 or mate[1][i] = - 2
 * 
 *   * Find augmenting path using DFS from node u: dfs(u)
 * 
 *   * Get max_matching: tot_mat = max_match(g)
 *      g is directed graph from left partition to right 
**/
struct dfs_matching {
	int N, M, iter = 0;
	std::vector<int> was;
	// NOTE: Set to - 2 if it is not in corresponding partition
	std::array<std::vector<int>, 2> mate;

	inline const int& size(const bool& z) const { return (z == 0) ? N : M; }
	inline int& size(const bool& z) { return (z == 0) ? N : M; }

	dfs_matching() {}
	dfs_matching(const int& N, const int& M) { init(N, M); }

	inline void init(const int& N, const int& M) {
		this -> N = N; this -> M = M;
		mate[0].assign(N, - 1);
		mate[1].assign(M, - 1);
		was.assign(N, - 1);
	}

	template <typename G> bool dfs_from(const G& g, const int& v) {
		was[v] = iter;
		for (const int& id : g.adj[v]) {
			if (g.is_ignore(id)) continue;
			const int u = g(v, id);
			if (mate[1][u] == - 1) {
				mate[0][v] = mate[1][u] = id;
				return true;
			}
		}
		for (const int& id : g.adj[v]) {
			if (g.is_ignore(id)) continue;
			const int u = g(v, id);
			if (mate[1][u] == - 2) continue;
			const int nxt = g(u, mate[1][u]);
			if (was[nxt] != iter && dfs_from(g, nxt)) {
				mate[0][v] = mate[1][u] = id;
				return true;
			}
		}
		return false;
	}

	template <typename G> bool dfs(const G& g, const int& v) {
		if (mate[0][v] != - 1) return false;
		return dfs_from(g, v);
	}

	template <typename G> int max_match(const G& g) {
		int res = 0;
		while (true) {
			int add = 0;
			for (int i = 0; i < N; i ++) {
				if (mate[0][i] == - 1 && dfs_from(g, i)) ++ add;
			}
			++ iter;
			if (add == 0) break;
			res += add;
		}
		return res;
	}
};

// Copied from: https://codeforces.com/blog/entry/118098
struct hopcroft_karp_matching {
	int N, M;
	// Try to maintain the layer network in the dinic flow
	std::vector<int> depth;
	std::array<std::vector<int>, 2> mate;

	inline const int& size(const bool& z) const { return (z == 0) ? N : M; }
	inline int& size(const bool& z) { return (z == 0) ? N : M; }

	hopcroft_karp_matching() {}
	hopcroft_karp_matching(const int& N, const int& M) { init(N, M); }

	inline void init(const int& N, const int& M) {
		this -> N = N; this -> M = M;
		depth.resize(N);
		mate[0].assign(N, - 1);
		mate[1].assign(M, - 1);
		ptr.resize(N);
	}

	template <typename G> bool bfs(const G& g) {
		std::fill(depth.begin(), depth.end(), N);
		std::vector<int> q(N);
		int beg = 0, end = 0;
		for (int v = 0; v < N; ++ v) {
			if (mate[0][v] == - 1) {
				// v not match with any node on the right
				depth[v] = 0;
				q[end ++] = v;
			}
		}
		bool rep = false;
		while (beg < end) {
			const auto& v = q[beg ++];
			for (const auto& id : g.adj[v]) {
				if (g.is_ignore(id)) continue;
				const int e = g(v, id);
				if (mate[1][e] == - 2) continue;
				if (mate[1][e] == - 1) {
					rep = true;
					continue;
				}
				const int u = g(e, mate[1][e]);
				if (depth[v] + 1 < depth[u]) {
					depth[u] = depth[v] + 1;
					q[end ++] = u;
				}
			}
		}
		return rep;
	}

	std::vector<int> ptr;
	template <typename G> bool dfs(const G& g, const int& v) {
		for (; ptr[v] < int(g.adj[v].size()); ++ ptr[v]) {
			const int& id = g.adj[v][ptr[v]];
			if (g.is_ignore(id)) continue;
			const int nxt = g(v, id);
			if (mate[1][nxt] == - 2) continue;
			if (mate[1][nxt] == - 1) {
				mate[1][nxt] = mate[0][v] = id;
				return true;
			}
			const int u = g(nxt, mate[1][nxt]);
			if (depth[u] == depth[v] + 1 && dfs(g, u)) {
				mate[1][nxt] = mate[0][v] = id;
				return true;
			}
		}
		return false;
	}

	template <typename G> int max_match(const G& g) {
		int tot_mat = 0;
		while (bfs(g)) {
			std::fill(ptr.begin(), ptr.end(), 0);
			// for (int i = 0; i < N; ++ i) {
			// 	ptr[i] = int(g.adj[i].size()) - 1;
			// }
			for (int v = 0; v < N; ++ v) {
				if (mate[0][v] == - 1 && dfs(g, v)) {
					++ tot_mat;
				}
			}
		}
		return tot_mat;
	}
};