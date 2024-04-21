template <typename key_t, typename val_t> struct radix_heap {
	// Idea: https://www.giaithuatlaptrinh.com/?tag=radix-heap
	// Code was copied from: https://people.ksp.sk/~kuko/ds/mat/radix-heap.pdf
	// Just apply for dijkstra algorithm
	static_assert(std::is_integral<key_t>::value, "only integers supported");
	static constexpr int32_t bit = sizeof(key_t) << 3;
	std::array<std::vector<std::pair<key_t, val_t>>, bit> data;
	size_t sz = size_t();
	key_t last = key_t();

	inline bool empty() const { return sz == 0; }
	inline size_t size() const { return sz; }
	static inline int32_t log_2(const int32_t& a) {
		return a ? 31 - __builtin_clz(a) : - 1;
	}
	static inline int32_t log_2(const int64_t& a) {
		return a ? 63 - __builtin_clzll(a) : - 1;
	}
	inline void push(const std::pair<key_t, val_t>& v) {
		assert(v.first >= last);
		++ sz;
		data[log_2(v.first ^ last) + 1].emplace_back(v);
	}
	template <typename... T>
	inline void emplace(const key_t& key, const T&... args) {
		assert(key >= last);
		++ sz;
		data[log_2(key ^ last) + 1].emplace_back(key, args...);
	}

	inline std::pair<key_t, val_t> pop() {
		if (data[0].empty()) {
			int idx = 1;
			while (data[idx].empty()) ++ idx;
			last = std::get<0>(*std::min_element(data[idx].begin(), data[idx].end()));
			for (const auto& p : data[idx]) {
				data[log_2(p.first ^ last) + 1].emplace_back(p);
			}
			data[idx].clear();
		}
		-- sz;
		auto res = data[0].back();
		data[0].pop_back();
		return res;
	}
	inline std::pair<key_t, val_t> top() {
		if (data[0].empty()) {
			int idx = 1;
			while (data[idx].empty()) ++ idx;
			last = std::get<0>(*std::min_element(data[idx].begin(), data[idx].end()));
			for (const auto& p : data[idx]) {
				data[log_2(p.first ^ last) + 1].emplace_back(p);
			}
			data[idx].clear();
		}
		auto res = data[0].back();
		return res;
	}
};

template <typename T, typename G>
std::vector<T> dijkstra_radix_heap(const G& g, const int &s, std::vector<int>& pe) {
	std::vector<T> dist(g.V, std::numeric_limits<T>::max());
	radix_heap<T, int> H;
	pe.assign(g.V, - 1);
	dist[s] = 0; H.emplace(dist[s], s);
	while (!H.empty()) {
		auto t = H.pop();
		if (dist[t.second] != t.first) continue;
		for (const int& id : g.adj[t.second]) {
			if (g.is_ignore(id)) continue;
			const auto& e = g.edges[id];
			int to = e.from ^ e.to ^ t.second;
			if (dist[t.second] + e.cost < dist[to]) {
				dist[to] = dist[t.second] + e.cost;
				pe[to] = id;
				H.emplace(dist[to], to);
			}
		}
	}
	return dist;
	// NOTE: returns numeric_limits<T>::max() if there's no path
}