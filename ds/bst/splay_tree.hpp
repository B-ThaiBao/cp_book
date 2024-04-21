/**
 * SPLAY_TREE !!!
 * 
 * Mostly inspired here:
 *   * https://codeforces.com/contest/1495/submission/109606482 (Tourist)
 *   * https://codeforces.com/contest/899/submission/44463457 (Tourist)
 * 
 * Application: Link-cut-tree !!!
 * 
 * Details:
 *   * Self-adjusting tree: after access node A, A will become root of tree
 *   * Beware: If don't want to change the structure of the tree. Don't use
 *   * function with node input.
 * 
 *   * Before do_operation on node or link to outer node, we must do_downdate_down() to
 *     maintain lazy propagation in exact subtree
 *   * When you update_up() for node, assume that we already do_downdate everything down
 *   * When child nodes receive information when do_downdate_down() from par node, we must
 *     apply_lazy() and can get query from all subtree (with root is child node)
 * 
 * Usage:
 *   * Make a struct splay_tree_node : public splay_tree_node_base<splay_tree_node>
 *   * based on CRTP method, which implements:
 *       * void do_downdate(): do_downdate_down something when lazy propagation.
 *       * void update(): update_up information when go_down and go_up again
 * 
 *   * Build tree: build()
 *   * Find: find(), find_first(), find_last(), find_implicit(), find_pos(), find_root()
 *   * Split: split(), split_implicit()
 *   * Merge: merge()
 *   * Combine: combine()
 *   * More and more .... 
 *   
 *   * To apply for single node u in the tree:
 *       * splay(u)
 *       * do_operation_on_node(u)
**/
template <typename splay_tree_node> struct splay_tree_node_base {
	std::array<splay_tree_node*, 2> c{nullptr, nullptr};
	splay_tree_node* par = nullptr;
	bool flip = false;
	int32_t sz = 1;

	splay_tree_node* derived_this() {
		return static_cast<splay_tree_node*>(this);
	}
	const splay_tree_node* derived_this() const {
		return static_cast<const splay_tree_node*>(this);
	}
	void do_downdate() { derived_this() -> do_downdate(); }
	void update() { derived_this() -> update(); }
	void reverse() {
		do_downdate();
		flip ^= true;
		std::swap(c[0], c[1]);
		update();
	}
	void downdate() {
		if (flip) {
			if (c[0] != nullptr){
				c[0] -> reverse();
			}
			if (c[1] != nullptr){
				c[1] -> reverse();
			}
			flip = false;
		}
		do_downdate();
	}

	friend bool is_root(splay_tree_node* v) {
		if (v == nullptr) return false;
		return (v -> par == nullptr || (v -> par -> c[0] != v && v -> par -> c[1] != v));
	}
	friend void rotate(splay_tree_node* v) {
		// TODO: Try to rotate this node up to par node and par node will go down
		auto u = v -> par;
		// assert(u != nullptr);
		u -> downdate(); v -> downdate();
		v -> par = u -> par;
		if (v -> par != nullptr) {
			if (v -> par -> c[0] == u) {
				v -> par -> c[0] = v;
			}
			if (v -> par -> c[1] == u) {
				v -> par -> c[1] = v;
			}
		}
		if (v == u -> c[0]) {
			u -> c[0] = v -> c[1];
			v -> c[1] = u;
		} else {
			u -> c[1] = v -> c[0];
			v -> c[0] = u;
		}
		u -> update(); v -> update();
	}
	friend void splay(splay_tree_node* v) {
		if (v == nullptr) return;
		while (!is_root(v)) {
			auto u = v -> par;
			if (!is_root(u)) {
				// NOTE: Zig-zig or zig-zag case
				if ((u -> c[0] == v) ^ (u -> par -> c[0] == u)) {
					rotate(v);
				} else {
					rotate(u);
				}
			}
			rotate(v);
		}
	}

	template <typename F>
	friend std::pair<splay_tree_node*, int> find(splay_tree_node* v, const F& f) {
		if (v == nullptr) return {nullptr, 0};
		splay(v); // For sure, v is root of tree
		int dir;
		while (true) {
			v -> downdate();
			dir = f(v);
			if (dir == 0) break;
			auto u = (dir == - 1 ? v -> c[0] : v -> c[1]);
			if (u == nullptr) break;
			v = u;
		}
		splay(v);
		return {v, dir};
	}
	friend splay_tree_node* find_first(splay_tree_node* v) {
		return find(v, [&](splay_tree_node*) { return - 1; }).first;
	}
	friend splay_tree_node* find_last(splay_tree_node* v) {
		return find(v, [&](splay_tree_node*) { return 1; }).first;
	}
	friend splay_tree_node* find_implicit(splay_tree_node* v, int k) {
		auto p = find(v, [&](splay_tree_node* u) {
			if (u -> c[0] != nullptr){
				if (u -> c[0] -> sz > k) return - 1; // go_left
				k -= (u -> c[0] -> sz);
			}
			if (k == 0) return 0; // Found
			-- k; return 1; // go_right
		});
		return (p.second == 0 ? p.first : nullptr);
	}
	friend int32_t find_pos(splay_tree_node* v) {
		splay(v);
		return (v -> c[0] != nullptr ? v -> c[0] -> sz : 0);
	}
	friend splay_tree_node* find_root(splay_tree_node* v) {
		splay(v); return v;
	}

	template <typename F>
	friend std::pair<splay_tree_node*, splay_tree_node*> split(splay_tree_node* v, const F& is_right) {
		if (v == nullptr) return {nullptr, nullptr};
		auto p = find(v, [&](splay_tree_node* u) { return is_right(u) ? - 1 : 1; } );
		v = p.first; v -> downdate();
		if (p.second == - 1) {
			auto u = v -> c[0];
			if (u == nullptr) return {nullptr, v};
			v -> c[0] = nullptr;
			u -> par = v -> par;
			u = find_last(u);
			v -> par = u; // NOTE: Something confuse here !!!
			v -> update();
			return {u, v};
		} else {
			auto u = v -> c[1];
			if (u == nullptr) return {v, nullptr};
			v -> c[1] = nullptr;
			v -> update();
			return {v, u};
		}
	}
	template <typename F>
	friend std::pair<splay_tree_node*, splay_tree_node*> split_implicit(splay_tree_node* v, int k, const F& is_right) {
		// TODO: Try to split a tree at most k node in left subtree
		return split(v, [&](splay_tree_node* u) {
			int hold_left = (u -> c[0] != nullptr ? u -> c[0] -> sz : 0) + 1;
			if (k < hold_left && is_right(u, k)) {
				return true;
			} else {
				k -= hold_left; return false;
			}
		});
	}
	friend std::pair<splay_tree_node*, splay_tree_node*> split_implicit(splay_tree_node* v, int k) {
		// TODO: Try to split a tree k node in left subtree
		return split(v, [&](splay_tree_node* u) {
			int hold_left = (u -> c[0] != nullptr ? u -> c[0] -> sz : 0) + 1;
			if (k < hold_left) {
				return true;
			} else {
				k -= hold_left; return false;
			}
		});
	}

	friend splay_tree_node* merge(splay_tree_node* v, splay_tree_node* u) {
		if (v == nullptr) return u;
		if (u == nullptr) return v;
		v = find_last(v);
		// assert(v -> c[1] == nullptr);
		splay(u);
		v -> downdate();
		v -> c[1] = u;
		v -> update();
		return v;
	}

	template <typename F>
	friend int count_left(splay_tree_node* v, const F& is_right) {
		if (v == nullptr) return 0;
		auto p = find(v, [&](splay_tree_node* u) { return is_right(u) ? - 1 : 1; });
		auto u = p.first;
		return (u -> c[0] != nullptr ? u -> c[0] -> sz : 0) + (p.second == 1);
	}
	template <typename F>
	friend splay_tree_node* insert(splay_tree_node* r, splay_tree_node* v, const F& is_right) {
		auto p = split(r, is_right);
		return merge(p.first, merge(v, p.second));
	}

	friend splay_tree_node* erase(splay_tree_node* v) {
		splay(v);
		v -> downdate();
		auto x = v -> c[0];
		auto y = v -> c[1];
		v -> c[0] = v -> c[1] = nullptr;
		auto z = merge(x, y);
		if (z != nullptr) z -> par = v -> par;
		v -> par = nullptr;
		v -> downdate(); v -> update();
		return z;
	}

	friend splay_tree_node* next(splay_tree_node* v) {
		splay(v);
		v -> downdate();
		if (v -> c[1] == nullptr) return nullptr;
		v = v -> c[1];
		while (v -> c[0] != nullptr) {
			v -> downdate();
			v = v -> c[0];
		}
		splay(v);
		return v;
	}
	friend splay_tree_node* previous(splay_tree_node* v) {
		splay(v);
		v -> downdate();
		if (v -> c[0] == nullptr) return nullptr;
		v = v -> c[0];
		while (v -> c[1] != nullptr) {
			v -> downdate();
			v = v -> c[1];
		}
		splay(v);
		return v;
	}
	friend int32_t size(splay_tree_node* v) {
		splay(v);
		return (v != nullptr ? v -> sz : 0);
	}

	template <typename F>
	friend splay_tree_node* split_and_merge(splay_tree_node* v, const int& L, const F& do_node) {
		auto p = split_implicit(v, L);
		do_node(p.first, p.second);
		p.first = merge(p.first, p.second);
		return p.first;
	}
	template <typename F>
	friend splay_tree_node* split_and_merge(splay_tree_node* v, const int& L, const int& R, const F& do_node) {
		auto p = split_implicit(v, L);
		auto u = split_implicit(p.second, R - L + 1);
		do_node(p.first, u.first, u.second);
		u.first = merge(u.first, u.second);
		p.first = merge(p.first, u.first);
		return p.first;
	}

	// NOTE: Try to travel all node and get id for each node
	template <typename F>
	friend void dfs_implicit(splay_tree_node* v, const F& f, const int& sz = 0) {
		if (v == nullptr) return;
		int cur = sz + (v -> c[0] != nullptr ? v -> c[0] -> sz : 0);
		f(cur, v);
		dfs_implicit(v -> c[0], f, sz);
		dfs_implicit(v -> c[1], f, cur + 1);
	}

	friend splay_tree_node* build_dfs(std::vector<splay_tree_node*>& nodes, const int& L, const int& R) {
		if (R < L) return nullptr;
		int mid = (L + R) >> 1;
		nodes[mid] -> c[0] = build_dfs(nodes, L, mid - 1);
		nodes[mid] -> c[1] = build_dfs(nodes, mid + 1, R);
		nodes[mid] -> update();
		return nodes[mid];
	}
	friend splay_tree_node* build_dfs(std::vector<splay_tree_node>& nodes, const int& L, const int& R) {
		if (R < L) return nullptr;
		int mid = (L + R) >> 1;
		nodes[mid].c[0] = build_dfs(nodes, L, mid - 1);
		nodes[mid].c[1] = build_dfs(nodes, mid + 1, R);
		nodes[mid].update();
		return &nodes[mid];
	}
	template <typename Vec> friend int build(Vec& nodes) {
		int N = int(nodes.size());
		build_dfs(nodes, 0, N - 1);
		return (N - 1) >> 1;
	}
};

struct splay_tree_node : public splay_tree_node_base<splay_tree_node> {
	

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

	void update() {
		sz = 1;
		if (c[0] != nullptr){
			c[0] -> par = this;
			sz += c[0] -> sz;
			// Take more information from c[0] ....
			
		}
		if (c[1] != nullptr){
			c[1] -> par = this;
			sz += c[1] -> sz;
			// Take more information from c[1] ....
			
		}
	}
};