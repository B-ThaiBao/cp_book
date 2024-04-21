/**
 * Copied from KATCL (Thanks)
 * 
 * Author: Johan Sannemo, Simon Lindholm
 * Date: 2016-12-15
 * License: CC0
 * 
 * Description: Finds a minimum vertex cover in a bipartite graph.
 *  The size is the same as the size of a maximum matching, and
 *  the complement is a maximum independent set.
 * 
 * Status: stress-tested
**/
template <typename G, typename Matching>
std::array<std::vector<int>, 2> min_vertex_cover(const G& g, const Matching& mat) {
	// NOTE: Matching must be max_match before here !!!
	std::vector<bool> found(mat.N, true), seen(mat.M, false);
	for (int i = 0; i < mat.N; ++ i) {
		if (mat.mate[1][i] >= 0) {
			found[g(i, mat.mate[1][i])] = false;
		}
	}
	std::vector<int> q; q.reserve(mat.N);
	std::array<std::vector<int>, 2> res;
	res[0].reserve(mat.N); res[1].reserve(mat.M);

	for (int i = 0; i < mat.N; ++ i) if (mat.mate[0][i] != - 2 && found[i]) q.push_back(i);
	// Try to DFS from all free node from the left
	while (!q.empty()) {
		int i = q.back(); q.pop_back();
		found[i] = true;
		for (const auto& id : g.adj[i]) {
			const int e = g(i, id);
			if (mat.mate[1][e] == - 2) continue;
			if (!seen[e] && mat.mate[1][e] != - 1) {
				seen[e] = true;
				q.push_back(g(e, mat.mate[1][e]));
			}
		}
	}
	// All vertex not found in the left and seen in the right are cover all edges
	for (int i = 0; i < mat.N; ++ i) if (mat.mate[0][i] != - 2 && !found[i]) res[0].push_back(i);
	for (int i = 0; i < mat.M; ++ i) if (mat.mate[1][i] != - 2 && seen[i]) res[1].push_back(i);
	return res;
}