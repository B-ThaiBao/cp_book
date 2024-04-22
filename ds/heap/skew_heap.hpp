template <typename skew_heap_node> struct skew_heap_node_base {
	std::array<skew_heap_node*, 2> c{nullptr, nullptr};
	// NOTE: Useful to cut from node to it child if only know child
	skew_heap_node* par = nullptr;

	skew_heap_node* derived_this() {
		return static_cast<skew_heap_node*>(this);
	}
	const skew_heap_node* derived_this() const {
		return static_cast<const skew_heap_node*>(this);
	}

	inline void downdate() { derived_this() -> do_downdate(); }

	template <typename Comp>
	friend skew_heap_node* merge(skew_heap_node* x, skew_heap_node* y, const Comp& comp) {
		// comp(x, y): return true if x is root and false if y is root of heap
		if (x == nullptr) return y;
		if (y == nullptr) return x;
		if (!comp(x, y)) std::swap(x, y);
		// Now x must be the root of heap
		x -> downdate();
		x -> c[1] = merge(x -> c[1], y, comp);
		if (x -> c[1] != nullptr) x -> c[1] -> par = x;
		std::swap(x -> c[0], x -> c[1]);
		return x;
	}

	friend skew_heap_node* find_root(skew_heap_node* v) {
		if (v == nullptr) return v;
		while (v -> par != nullptr) v = v -> par;
		return v;
	}

	template <typename Comp>
	friend skew_heap_node* pop(skew_heap_node* x, const Comp& comp) {
		if (x == nullptr) return x;
		return merge(x -> c[0], x -> c[1], comp);
	}
};

struct skew_heap_node : public skew_heap_node_base<skew_heap_node> {
	

	void apply() {
		// Apply lazy propagation from par node
		
	}

	void do_downdate() {
		// // Push everything else except flip ....
		// if (add != 0){
		// 	if (c[0] != nullptr){
		// 		c[0] -> apply(add);
		// 	}
		// 	if (c[1] != nullptr){
		// 		c[1] -> apply(add);
		// 	}
		// 	add = 0;
		// }
	}
};

template <typename T, typename G>
std::vector<T> dijkstra_skew_heap(const G& g, const int &s, std::vector<int>& pe) {
	struct alignas(32) dijkstra_skew_heap_node {
		T key= T();
		int p = - 1, L = - 1, R = - 1;
	};
	std::vector<dijkstra_skew_heap_node> nodes(g.V);
	int root = s;

	auto merge = [&](auto&& merge, int x, int y) -> int {
		if (x == - 1) return y;
		if (y == - 1) return x;
		if (nodes[x].key > nodes[y].key) std::swap(x, y);
		nodes[x].R = merge(merge, nodes[x].R, y);
		nodes[nodes[x].R].p = x;
		std::swap(nodes[x].L, nodes[x].R);
		return x;
	};

	static constexpr T INF_COST = std::numeric_limits<T>::max();
	std::vector<T> dist(g.V, INF_COST);
	pe.assign(g.V, - 1); dist[s] = 0;
	while (root != - 1) {
		auto cst = nodes[root].key;
		int cur = root;
		root = merge(merge, nodes[root].L, nodes[root].R);
		if (dist[cur] != cst) continue;
		for (const int& id : g.adj[cur]) {
			if (g.is_ignore(id)) continue;
			const auto& e = g.edges[id];
			int to = e.from ^ e.to ^ cur;
			if (dist[to] == INF_COST) {
				dist[to] = dist[cur] + e.cost;
				pe[to] = id;
				nodes[to].key = dist[to];
				root = merge(merge, root, to);
			} else if (dist[cur] + e.cost < dist[to]) {
				dist[to] = dist[cur] + e.cost;
				pe[to] = id;
				nodes[to].key = dist[to];
				int par = nodes[to].p;
				if (par == - 1) continue;
				if (nodes[par].key <= dist[to]) continue;
				(nodes[par].L == to ? nodes[par].L : nodes[par].R) = - 1;
				nodes[to].p = - 1;
				root = merge(merge, root, to);
			}
		}
	}
	return dist;
	// NOTE: returns numeric_limits<T>::max() if there's no path
}
