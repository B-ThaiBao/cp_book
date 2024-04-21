// NOTE: In mst, all weighted edges on the path from node u to v is the minimum compare to another path in the real graph
//       In other word, the min of max weight from u to v is the unique path in the min span tree (mst)
// See:  https://codeforces.com/blog/entry/100666?#comment-893646
template <typename G, typename T, typename Comp = std::less<>>
std::vector<int> build_span_tree(const G& g, T& sum_cost, const Comp& comp = Comp()) {
	std::vector<int> order(g.edges.size());
	std::iota(order.begin(), order.end(), 0);
	std::sort(order.begin(), order.end(), [&](const int& a, const int& b) {
		return comp(g.edges[a].cost, g.edges[b].cost);
	});
	disjoint_set dsu(g.V);
	std::vector<int> res; res.reserve(order.size());
	for (const auto& id : order) {
		if (g.is_ignore(id)) continue;
		const auto& e = g.edges[id];
		if (dsu.merge(e.from, e.to)) {
			res.push_back(id); sum_cost += e.cost;
		}
	}
	return res;
}
