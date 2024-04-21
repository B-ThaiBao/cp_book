/**
 * TREAP !!!
 * 
 * Mostly inspired here:
 *   * https://codeforces.com/contest/1737/submission/175043607 (Tourist)
 *   * https://codeforces.com/contest/1737/submission/175037024 (Ecnerwala)
 *   * https://codeforces.com/contest/1852/submission/215285692 (Ecnerwala)
 *   * https://usaco.guide/adv/treaps?lang=cpp (Benq)
 *   * https://github.com/OpenGenus/cosmos/pull/2407/files
 *   * https://www.youtube.com/watch?v=6x0UlIBLRsc (Theory and code)
 *   * https://www.youtube.com/watch?v=erKlLEXLKyY&t=3468s (Explain well !!!)
 * 
 * Details:
 *   * Before do_operation on node or link to outer node, we must do_downdate_down() to
 *     maintain lazy propagation in exact subtree
 *   * When you update_up() for node, assume that we already do_downdate everything down
 *   * When child nodes receive information when do_downdate_down() from par node, we must
 *     apply_lazy() and can get query from all subtree (with root is child node)
 * 
 *  Usage:
 *   * Make a struct treap_node : public treap_node_base<treap_node> based on CRTP method,
 *     which implements:
 *      * void do_downdate(): do_downdate_down something when lazy propagation.
 *      * void update(): update_up information when go_down and go_up again
 * 
 *   * Build treap: build_by_stack() or build_by_heap()
 *   * Find: find(), find_first(), find_last(), find_implicit(), find_pos(), find_root()
 *   * Split: split(), split_implicit()
 *   * Merge: merge()
 *   * Combine: combine()
 *   * More and more .... 
**/
std::mt19937 treap_random(std::chrono::steady_clock::now().time_since_epoch().count());
template <typename treap_node> struct treap_node_base {
	std::array<treap_node*, 2> c{nullptr, nullptr};
	treap_node* par = nullptr;
	bool flip = false;
	int32_t sz = 1;
	std::mt19937::result_type priority = treap_random();

	treap_node* derived_this() {
		return static_cast<treap_node*>(this);
	}
	const treap_node* derived_this() const {
		return static_cast<const treap_node*>(this);
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

	// NOTE: diretion :: - 1: go left, 0: found, 1: go right
	template <typename F>
	friend std::pair<treap_node*, int32_t> find(treap_node* v, const F& where) {
		if (v == nullptr) return {nullptr, 0};
		int direct;
		while (true){
			v -> downdate();
			direct = where(v);
			if (direct == 0){
				// Found
				break;
			}
			auto u = (direct == - 1 ? v -> c[0] : v -> c[1]);
			if (u == nullptr) break;
			v = u;
		}
		return {v, direct};
	}
	friend treap_node* find_first(treap_node* v) {
		// TODO: Return the leftmost node in tree (first element in array)
		return find(v, [&](const auto& v){ return - 1; }).first;
	}
	friend treap_node* find_last(treap_node* v) {
		// TODO: Return the leftmost node in tree (first element in array)
		return find(v, [&](const auto& v){ return 1; }).first;
	}
	friend treap_node* find_implicit(treap_node* v, int k) {
		auto p = find(v, [&](treap_node* u) -> int{
			if (u -> c[0] != nullptr){
				if (u -> c[0] -> sz > k) return - 1; // go_left
				k -= (u -> c[0] -> sz);
			}
			if (k == 0) return 0; // Found
			-- k; return 1; // go_right
		});
		return (p.second == 0 ? p.first : nullptr);
	}
	friend int find_pos(treap_node* v) {
		int k = (v -> c[0] != nullptr ? v -> c[0] -> sz : 0);
		while (v -> par != nullptr){
			if (v == v -> par -> c[1]){
				++ k;
				if (v -> par -> c[0] != nullptr){
					k += v -> par -> c[0] -> sz;
				}
			}
			v = v -> par; // Go_to par node
		}
		return k;
	}
	friend treap_node* find_root(treap_node* v) {
		while (v -> par != nullptr){
			v = v -> par;
		}
		return v;
	}

	template <typename F>
	friend std::pair<treap_node*, treap_node*> split(treap_node* v, const F& is_right) {
		if (v == nullptr) return {nullptr, nullptr};
		v -> downdate();
		if (is_right(v)){
			std::pair<treap_node*, treap_node*> p = split(v -> c[0], is_right);
			if (p.first != nullptr){
				p.first -> par = nullptr;
			}
			v -> c[0] = p.second;
			v -> update();
			return {p.first, v};
		}
		else{
			std::pair<treap_node*, treap_node*> p = split(v -> c[1], is_right);
			v -> c[1] = p.first;
			if (p.second != nullptr){
				p.second -> par = nullptr;
			}
			v -> update();
			return {v, p.second};
		}
	}
	template <typename F>
	friend std::pair<treap_node*, treap_node*> split_implicit(treap_node* v, const int& k, const F& is_right) {
		// TODO: Try to split a tree k node in left subtree
		if (v == nullptr) return {nullptr, nullptr};
		v -> downdate();
		int hold_left = (v -> c[0] != nullptr ? v -> c[0] -> sz : 0) + 1;
		if (k < hold_left && is_right(v, k)){
			std::pair<treap_node*, treap_node*> p = split_implicit(v -> c[0], k, is_right);
			if (p.first != nullptr){
				p.first -> par = nullptr;
			}
			v -> c[0] = p.second;
			v -> update();
			return {p.first, v};
		}
		else{
			std::pair<treap_node*, treap_node*> p = split_implicit(v -> c[1], k - hold_left, is_right);
			v -> c[1] = p.first;
			if (p.second != nullptr){
				p.second -> par = nullptr;
			}
			v -> update();
			return {v, p.second};
		}
	}
	friend std::pair<treap_node*, treap_node*> split_implicit(treap_node* v, const int& k) {
		// TODO: Split a tree with k-node in left subtree
		return split_implicit(v, k, [](treap_node*, const int& x) { return true; });
	}

	friend treap_node* merge(treap_node* v, treap_node* u) {
		if (v == nullptr) return u;
		if (u == nullptr) return v;
		if (v -> priority > u -> priority){
			v -> downdate();
			v -> c[1] = merge(v -> c[1], u);
			v -> update();
			return v;
		}
		else{
			u -> downdate();
			u -> c[0] = merge(v, u -> c[0]);
			u -> update();
			return u;
		}
	}

	// NOTE: This doesn't allow combine_implicit and must have key member to combine
	friend treap_node* combine(treap_node* u, treap_node* v) {
		if (u == nullptr) return v;
		if (v == nullptr) return u;
		if (u -> priority < v -> priority) std::swap(u, v);
		auto N = split(v, [&](const auto& x) -> bool{
			return x -> key > u -> key;
		});
		u -> c[0] = combine(u -> c[0], N.first);
		u -> c[1] = combine(u -> c[1], N.second);
		return u;
	}

	template <typename F>
	friend int count_left(treap_node* v, const F& is_right) {
		if (v == nullptr) return 0;
		v -> downdate();
		if (is_right(v)) return count_left(v -> c[0], is_right);
		return (v -> c[0] != nullptr ? v -> c[0] -> sz : 0) + count_left(v -> c[1], is_right);
	}
	template <typename F>
	friend treap_node* insert(treap_node* u, treap_node* v, const F& is_right) {
		std::pair<treap_node*, treap_node*> p = split(v, is_right);
		return merge(p.first, merge(u, p.second));
	}

	friend treap_node* erase(treap_node* v) {
		// TODO: Try to erase v out of the tree --> Return new root
		v -> downdate();
		treap_node* x = v -> c[0];
		treap_node* y = v -> c[1];
		treap_node* p = v -> par;
		v -> c[0] = v -> c[1] = v -> par = nullptr;
		v -> downdate();
		v -> update(); // now v might be reusable...
		treap_node* z = merge(x, y);
		if (p == nullptr){
			if (z != nullptr){
				z -> par = nullptr;
			}
			return z;
		}
		if (p -> c[0] == v){
			p -> c[0] = z;
		}
		if (p -> c[1] == v){
			p -> c[1] = z;
		}
		while (true){
			p -> downdate();
			p -> update();
			if (p -> par == nullptr) break;
			p = p -> par;
		}
		return p;
	}

	friend treap_node* next(treap_node* v){
		// TODO: Return next element in the array of this node
		if (v -> c[1] == nullptr){
			while (v -> par != nullptr && v -> par -> c[1] == v){
				v = v -> par;
			}
			return v -> par;
		}
		v -> do_downdate();
		v = v -> c[1];
		while (v -> c[0] != nullptr){
			v -> do_downdate();
			v = v -> c[0];
		}
		return v;
	}
	friend treap_node* previous(treap_node* v){
		// TODO: Return previous element in the array of this node
		if (v -> c[0] == nullptr){
			while (v -> par != nullptr && v -> par -> c[0] == v){
				v = v -> par;
			}
			return v -> par;
		}
		v -> do_downdate();
		v = v -> c[0];
		while (v -> c[1] != nullptr){
			v -> do_downdate();
			v = v -> c[1];
		}
		return v;
	}
	friend int size(treap_node* v){
		return (v != nullptr ? v -> sz : 0);
	}

	template <typename F>
	friend treap_node* split_and_merge(treap_node* v, const int& L, const F& do_node) {
		auto p = split_implicit(v, L);
		do_node(p.first, p.second);
		p.first = merge(p.first, p.second);
		return p.first;
	}
	template <typename F>
	friend treap_node* split_and_merge(treap_node* v, const int& L, const int& R, const F& do_node) {
		auto p = split_implicit(v, L);
		auto u = split_implicit(p.second, R - L + 1);
		do_node(p.first, u.first, u.second);
		u.first = merge(u.first, u.second);
		p.first = merge(p.first, u.first);
		return p.first;
	}

	friend void heapify(treap_node* v) {
		if (v == nullptr) return;
		while (true) {
			auto ma = v;
			if (v -> c[0] != nullptr && v -> c[0] -> priority < ma -> priority) {
				ma = v -> c[0];
			}
			if (v -> c[1] != nullptr && v -> c[1] -> priority < ma -> priority) {
				ma = v -> c[1];
			}
			if (ma != v) {
				// TODO: We can swap priority to make the tree better
				std::swap(ma -> priority, v -> priority);
				v = ma;
			}
			else break;
		}
	}
	friend treap_node* build_dfs(std::vector<treap_node*>& A, const int& L, const int& R) {
		if (R < L) return nullptr;
		int mid = (L + R) >> 1;
		A[mid] -> c[0] = build_dfs(A, L, mid - 1);
		A[mid] -> c[1] = build_dfs(A, mid + 1, R);
		heapify(A[mid]);
		A[mid] -> update();
		return A[mid];
	}
	friend treap_node* build_dfs(std::vector<treap_node>& A, const int& L, const int& R) {
		if (R < L) return nullptr;
		int mid = (L + R) >> 1;
		A[mid].c[0] = build_dfs(A, L, mid - 1);
		A[mid].c[1] = build_dfs(A, mid + 1, R);
		heapify(&A[mid]);
		A[mid].update();
		return &A[mid];
	}
	template <typename Vec> friend int build_by_heap(Vec& A) {
		int N = int(A.size());
		if (N == 0) return - 1;
		build_dfs(A, 0, N - 1);
		return (N - 1) >> 1;
	}

	friend void build_stack_dfs(treap_node* v) {
		if (v -> c[0] != nullptr) build_stack_dfs(v -> c[0]);
		if (v -> c[1] != nullptr) build_stack_dfs(v -> c[1]);
		v -> update();
	}
	friend int build_by_stack(std::vector<treap_node*>& nodes) {
		int N = int(nodes.size());
		if (N == 0) return - 1;
		std::vector<int> stk; stk.reserve(N);
		std::vector<int> par_idx(N);
		for (int i = 0; i < N; ++ i){
			while (!stk.empty() && nodes[stk.back()] -> priority < nodes[i] -> priority){
				stk.pop_back();
			}
			if (!stk.empty()){
				nodes[i] -> par = nodes[stk.back()];
				par_idx[i] = stk.back();
			}
			stk.do_downdate_back(i);
		}
		stk.clear();
		for (int i = N - 1; i >= 0; -- i){
			while (!stk.empty() && nodes[stk.back()] -> priority < nodes[i] -> priority){
				stk.pop_back();
			}
			if (!stk.empty()){
				if (nodes[i] -> par == nullptr || nodes[i] -> par -> priority > nodes[stk.back()] -> priority){
					nodes[i] -> par = nodes[stk.back()];
					par_idx[i] = stk.back();
				}
			}
			stk.do_downdate_back(i);
		}

		int id = 0;
		for (int i = 0; i < N; ++ i){
			if (nodes[i] -> priority > nodes[id] -> priority) id = i;
			if (nodes[i] -> par != nullptr){
				// Update nodes[i] -> c[0 ... 1]
				if (par_idx[i] < i){
					nodes[i] -> par -> c[1] = nodes[i];
				}
				else{
					nodes[i] -> par -> c[0] = nodes[i];
				}
			}
		}
		build_stack_dfs(nodes[id]);
		return id;
	}
	friend int build_by_stack(std::vector<treap_node>& nodes) {
		int N = int(nodes.size());
		if (N == 0) return - 1;
		std::vector<int> stk; stk.reserve(N);
		std::vector<int> par_idx(N);
		for (int i = 0; i < N; ++ i){
			while (!stk.empty() && nodes[stk.back()].priority < nodes[i].priority){
				stk.pop_back();
			}
			if (!stk.empty()){
				nodes[i].par = &nodes[stk.back()];
				par_idx[i] = stk.back();
			}
			stk.do_downdate_back(i);
		}
		stk.clear();
		for (int i = N - 1; i >= 0; -- i){
			while (!stk.empty() && nodes[stk.back()].priority < nodes[i].priority){
				stk.pop_back();
			}
			if (!stk.empty()){
				if (nodes[i].par == nullptr || nodes[i].par -> priority > nodes[stk.back()].priority){
					nodes[i].par = &nodes[stk.back()];
					par_idx[i] = stk.back();
				}
			}
			stk.do_downdate_back(i);
		}

		int id = 0;
		for (int i = 0; i < N; ++ i){
			if (nodes[i].priority > nodes[id].priority) id = i;
			if (nodes[i].par != nullptr){
				// Update nodes[i] -> c[0 ... 1]
				if (par_idx[i] < i){
					nodes[i].par -> c[1] = &nodes[i];
				}
				else{
					nodes[i].par -> c[0] = &nodes[i];
				}
			}
		}
		build_stack_dfs(&nodes[id]);
		return id;
	}

	// NOTE: Try to travel all node and get id for each node
	template <typename F>
	friend void dfs_implicit(treap_node* v, const F& f, const int& sz = 0) {
		if (v == nullptr) return;
		int cur = sz + size(v -> c[0]);
		f(cur, v);
		dfs_implicit(v -> c[0], f, sz);
		dfs_implicit(v -> c[1], f, cur + 1);
	}
};

struct treap_node : public treap_node_base<treap_node> {
	

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