/**
 * GRAPH !!!
 * 
 * This graph template store the edge by using indicies so that
 * can apply for multiple edge with same node. For this reason,
 * lots of algorithm on the graph template can be written much
 * more easy and exact based on generic properties.
**/
template <typename T> struct graph {
	struct edge_base {
		int from;
		int to;
	};
	int V;
	std::vector<std::vector<int>> adj;
	std::vector<T> edges;

	graph() {}
	graph(const int& N) : V(N), adj(N) {}
	graph(const int& N, const int& M) : V(N), adj(N) {
		edges.reserve(M); // Avoid relocation
	}
};

template <typename T> struct undigraph : public graph<T> {
	using graph<T>::graph;
	using graph<T>::V;
	using graph<T>::adj;
	using graph<T>::edges;

	int add_edge(const T& c) {
		int e = int(edges.size());
		adj[c.from].push_back(e);
		adj[c.to].push_back(e);
		edges.emplace_back(c);
		return e;
	}
	template <typename... Args>
	int add_edge(const int& frm, const int& to, const Args&... args) {
		int e = int(edges.size());
		adj[frm].push_back(e);
		adj[to].push_back(e);
		edges.emplace_back(T{frm, to, args...});
		return e;
	}

	inline int operator () (const int& frm, const int& i) const {
		return frm ^ edges[i].from ^ edges[i].to;
	}
	inline int operator () (const int& frm, const T& e) const {
		return frm ^ e.from ^ e.to;
	}

	static constexpr bool is_ignore(const int& id) { return false; }
};

template <typename T> struct digraph : public graph<T> {
	using graph<T>::graph;
	using graph<T>::V;
	using graph<T>::adj;
	using graph<T>::edges;

	int add_edge(const T& c) {
		int e = int(edges.size());
		adj[c.from].push_back(e);
		edges.emplace_back(c);
		return e;
	}
	template <typename... Args>
	int add_edge(const int& frm, const int& to, const Args&... args) {
		int e = int(edges.size());
		adj[frm].push_back(e);
		edges.emplace_back(T{frm, to, args...});
		return e;
	}

	inline int operator () (const int& frm, const int& i) const {
		return frm ^ edges[i].from ^ edges[i].to;
	}
	inline int operator () (const int& frm, const T& e) const {
		return frm ^ e.from ^ e.to;
	}

	static constexpr bool is_ignore(const int& id) { return false; }
};

template <typename T, typename F = std::function<bool(int)>>
struct condi_undigraph : public graph<T> {
	using graph<T>::graph;
	using graph<T>::V;
	using graph<T>::adj;
	using graph<T>::edges;
	F ignore = nullptr;

	int add_edge(const T& c) {
		int e = int(edges.size());
		adj[c.from].push_back(e);
		adj[c.to].push_back(e);
		edges.emplace_back(c);
		return e;
	}
	template <typename... Args>
	int add_edge(const int& frm, const int& to, const Args&... args) {
		int e = int(edges.size());
		adj[frm].push_back(e);
		adj[to].push_back(e);
		edges.emplace_back(T{frm, to, args...});
		return e;
	}

	inline int operator () (const int& frm, const int& i) const {
		return frm ^ edges[i].from ^ edges[i].to;
	}
	inline int operator () (const int& frm, const T& e) const {
		return frm ^ e.from ^ e.to;
	}
	
	template <typename U> inline void set_ignore(const U& f) { ignore = f; }
	template <typename U> inline void clear_ignore() { ignore = nullptr; }

	constexpr bool is_ignore(const int& id) const {
		if (ignore == nullptr) return false;
		return ignore(id);
	}
};

template <typename T, typename F = std::function<bool(int)>>
struct condi_digraph : public graph<T> {
	using graph<T>::graph;
	using graph<T>::V;
	using graph<T>::adj;
	using graph<T>::edges;
	F ignore = nullptr;

	int add_edge(const T& c) {
		int e = int(edges.size());
		adj[c.from].push_back(e);
		adj[c.to].push_back(e);
		edges.emplace_back(c);
		return e;
	}
	template <typename... Args>
	int add_edge(const int& frm, const int& to, const Args&... args) {
		int e = int(edges.size());
		adj[frm].push_back(e);
		edges.emplace_back(T{frm, to, args...});
		return e;
	}

	inline int operator () (const int& frm, const int& i) const {
		return frm ^ edges[i].from ^ edges[i].to;
	}
	inline int operator () (const int& frm, const T& e) const {
		return frm ^ e.from ^ e.to;
	}

	template <typename U> inline void set_ignore(const U& f) { ignore = f; }
	template <typename U> inline void clear_ignore() { ignore = nullptr; }

	constexpr bool is_ignore(const int& id) const {
		if (ignore == nullptr) return false;
		return ignore(id);
	}
};