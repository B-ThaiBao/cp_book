/**
 * LINK_CUT_TREE (LCT) !!!
 * 
 * Idea is from:
 *   * https://codeforces.com/blog/entry/80383
 *   * https://codeforces.com/blog/entry/67637 (Subtree with LCT but can use ETT)
 *   * https://codeforces.com/blog/entry/80145 (lazy prop with LCT but can use top_tree)
 * 
 * Code was copied from:
 *   * https://codeforces.com/contest/1495/submission/109606482
 * 
 * More details:
 *   * LCT is dynamic HLD with all nodes in same comp connected by heavy/light edges.
 *   * Trivially, LCT is bottom-up tree --> just go up only
 * 
 * Usage:
 *   * Link_root u, v: v must be root of comp and attach v to child of u
 *   * Link u, v: u, v are abitrary nodes --> u will be root of comp
 *       * BEWARE: If false, v will be become root of old tree
 * 
 *   * Cut_root u, v: v must be the root of comp and delete edge (u -> v)
 *   * Cut_root u: delete u out of chain to par of u
 *   * Cut u, v: u, v are abitrary nodes --> delete a link if exist
 * 
 *   * Change u to root of that tree: make_lct_root(u)
 * 
 *   * lca, is_ancestor is no change
 * 
 *   * For rooted tree, apply PATH from node u to root is easy:
 *       * expose(u);
 *       * do_operation_on_node(u);
 *       * u -> downdate();
 *       * u -> update();
 * 
 *   * To apply for single node u in the tree:
 *       * splay(u)
 *       * do_operation_on_node(u)
 * 
 * WARNING:
 *   * For easy to use, this tree just apply for operation on PATH, alternative ways:
 *       * For subtree: euler_tour_tree
 *       * For subtree + path: top_tree (https://codeforces.com/blog/entry/103726)
**/
template <typename link_cut_node> struct link_cut_node_base{
	std::array<link_cut_node*, 2> c{nullptr, nullptr};
	link_cut_node* par = nullptr;
	bool flip = false;
	int32_t sz = 1;

	link_cut_node* derived_this() {
		return static_cast<link_cut_node*>(this);
	}
	const link_cut_node* derived_this() const {
		return static_cast<const link_cut_node*>(this);
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

	friend bool is_root(link_cut_node* v) {
		if (v == nullptr) return false;
		return (v -> par == nullptr || (v -> par -> c[0] != v && v -> par -> c[1] != v));
	}
	friend void rotate(link_cut_node* v) {
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
	friend void splay(link_cut_node* v) {
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
	friend std::pair<link_cut_node*, int> find(link_cut_node* v, const F& f) {
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
	friend link_cut_node* find_first(link_cut_node* v) {
		return find(v, [&](link_cut_node*) { return - 1; }).first;
	}
	friend link_cut_node* find_last(link_cut_node* v) {
		return find(v, [&](link_cut_node*) { return 1; }).first;
	}
	friend link_cut_node* find_implicit(link_cut_node* v, int k) {
		auto p = find(v, [&](link_cut_node* u) {
			if (u -> c[0] != nullptr){
				if (u -> c[0] -> sz > k) return - 1; // go_left
				k -= (u -> c[0] -> sz);
			}
			if (k == 0) return 0; // Found
			-- k; return 1; // go_right
		});
		return (p.second == 0 ? p.first : nullptr);
	}
	friend int32_t find_pos(link_cut_node* v) {
		splay(v);
		return (v -> c[0] != nullptr ? v -> c[0] -> sz : 0);
	}
	friend link_cut_node* find_root(link_cut_node* v) {
		splay(v); return v;
	}

	template <typename F>
	friend std::pair<link_cut_node*, link_cut_node*> split(link_cut_node* v, const F& is_right) {
		if (v == nullptr) return {nullptr, nullptr};
		auto p = find(v, [&](link_cut_node* u) { return is_right(u) ? - 1 : 1; } );
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
	friend std::pair<link_cut_node*, link_cut_node*> split_implicit(link_cut_node* v, int k, const F& is_right) {
		// TODO: Try to split a tree at most k node in left subtree
		return split(v, [&](link_cut_node* u) {
			int hold_left = (u -> c[0] != nullptr ? u -> c[0] -> sz : 0) + 1;
			if (k < hold_left && is_right(u, k)) {
				return true;
			} else {
				k -= hold_left; return false;
			}
		});
	}
	friend std::pair<link_cut_node*, link_cut_node*> split_implicit(link_cut_node* v, int k) {
		// TODO: Try to split a tree k node in left subtree
		return split(v, [&](link_cut_node* u) {
			int hold_left = (u -> c[0] != nullptr ? u -> c[0] -> sz : 0) + 1;
			if (k < hold_left) {
				return true;
			} else {
				k -= hold_left; return false;
			}
		});
	}

	friend link_cut_node* merge(link_cut_node* v, link_cut_node* u) {
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
	friend int count_left(link_cut_node* v, const F& is_right) {
		if (v == nullptr) return 0;
		auto p = find(v, [&](link_cut_node* u) { return is_right(u) ? - 1 : 1; });
		auto u = p.first;
		return (u -> c[0] != nullptr ? u -> c[0] -> sz : 0) + (p.second == 1);
	}
	template <typename F>
	friend link_cut_node* insert(link_cut_node* r, link_cut_node* v, const F& is_right) {
		auto p = split(r, is_right);
		return merge(p.first, merge(v, p.second));
	}

	friend link_cut_node* erase(link_cut_node* v) {
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

	friend link_cut_node* next(link_cut_node* v) {
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
	friend link_cut_node* previous(link_cut_node* v) {
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
	friend int32_t size(link_cut_node* v) {
		splay(v);
		return (v != nullptr ? v -> sz : 0);
	}

	template <typename F>
	friend link_cut_node* split_and_merge(link_cut_node* v, const int& L, const F& do_node) {
		auto p = split_implicit(v, L);
		do_node(p.first, p.second);
		p.first = merge(p.first, p.second);
		return p.first;
	}
	template <typename F>
	friend link_cut_node* split_and_merge(link_cut_node* v, const int& L, const int& R, const F& do_node) {
		auto p = split_implicit(v, L);
		auto u = split_implicit(p.second, R - L + 1);
		do_node(p.first, u.first, u.second);
		u.first = merge(u.first, u.second);
		p.first = merge(p.first, u.first);
		return p.first;
	}

	// NOTE: Try to travel all node and get id for each node
	template <typename F>
	friend void dfs_implicit(link_cut_node* v, const F& f, const int& sz = 0) {
		if (v == nullptr) return;
		int cur = sz + (v -> c[0] != nullptr ? v -> c[0] -> sz : 0);
		f(cur, v);
		dfs_implicit(v -> c[0], f, sz);
		dfs_implicit(v -> c[1], f, cur + 1);
	}

	friend link_cut_node* build_dfs(std::vector<link_cut_node*>& nodes, const int& L, const int& R) {
		if (R < L) return nullptr;
		int mid = (L + R) >> 1;
		nodes[mid] -> c[0] = build_dfs(nodes, L, mid - 1);
		nodes[mid] -> c[1] = build_dfs(nodes, mid + 1, R);
		nodes[mid] -> update();
		return nodes[mid];
	}
	friend link_cut_node* build_dfs(std::vector<link_cut_node>& nodes, const int& L, const int& R) {
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
	friend void expose(link_cut_node* v) {
		link_cut_node* r = nullptr;
		auto u = v;
		while (u != nullptr) {
			splay(u);
			u -> downdate();
			u -> c[1] = r;
			u -> update();
			r = u; u = u -> par;
		}
		splay(v);
		// assert(v -> par == nullptr);
		// NOTE: After this, v become root of that splay_tree with no right subtree and no parent
	}

	friend link_cut_node* find_lct_root(link_cut_node* v) {
		expose(v);
		// NOTE: After expose, v will have only left subtree with root is the leftmost one
		return find_first(v);
	}
	friend void make_lct_root(link_cut_node* v) {
		if (v == nullptr) return;
		// NOTE: After expose(v): the heavy path will be v to root (similar to one chain of HLD)
		// So we just reverse the array in this chain and v become the root of chain
		expose(v); v -> reverse();
	}

	friend bool link(link_cut_node* u, link_cut_node* v) {
		if (u == v) return false;
		make_lct_root(v);
		expose(u);
		if (v -> par != nullptr) return false;
		v -> par = u;
		return true;
	}
	friend bool link_root(link_cut_node* u, link_cut_node* v) {
		if (u == v) return false;
		splay(v);
		if (v -> par != nullptr || v -> c[0] != nullptr) {
			// v is not a root of tree, the condition is not met
			return false;
		}
		expose(u);
		if (v -> par != nullptr) return false;
		v -> par = u; return true;
	}

	friend bool cut(link_cut_node* u, link_cut_node* v) {
		if (u == v) return false;
		expose(u);
		splay(v);
		if (v -> par != u) {
			std::swap(u, v);
			expose(u);
			splay(v);
			if (v -> par != u) return false;
		}
		v -> par = nullptr; return true;
	}
	friend bool cut_root(link_cut_node* u, link_cut_node* v) {
		if (u == v) return false;
		expose(u);
		splay(v);
		if (v -> par != u) {
			return false; // v is not root
		}
		v -> par = nullptr; return true;
	}
	friend bool cut_root(link_cut_node* v) {
		expose(v);
		v -> downdate();
		if (v -> c[0] == nullptr) {
			// NOTE: v is root --> no delete anything
			return false;
		}
		v -> c[0] -> par = nullptr;
		v -> c[0] = nullptr;
		v -> update();
		return true;
	}

	friend bool same_comp(link_cut_node* v, link_cut_node* u) {
		expose(v); expose(u);
		return v -> par != nullptr;
	}
	friend link_cut_node* find_lca(link_cut_node* v, link_cut_node* u) {
		if (u == v) return u;
		expose(v); expose(u);
		if (v -> par == nullptr) return - 1;
		splay(v);
		if (v -> par == nullptr) return v;
		return v -> par;
	}
	friend bool is_ancestor(link_cut_node* v, link_cut_node* u) {
		if (u == v) return true;
		expose(u); splay(v);
		return v -> par == nullptr && u -> par != nullptr;
	}
};

struct link_cut_node : public link_cut_node_base<link_cut_node> {
	

	void apply() {
		// Apply lazy propagation from par node
		
	}

	void do_downdate() {
		// // do_downdate everything else except flip ....
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