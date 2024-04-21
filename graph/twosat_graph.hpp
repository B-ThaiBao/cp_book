/**
 * TWOSAT (2 - SAT) !!!
 * 
 * Tutorial: https://codeforces.com/blog/entry/92977 (WELL !!!)
 * Code implementation: KATCL and more ...
 * 
 * Description: Solve sat prob in linear time
 * Negated variables are represented by bit-inversions (~ x).
 * Usage:
 *   * twosat_graph sat(num of var)
 *   * sat.either(a, b): a or b
 *   * sat.implies(a, b): a -> b
 *   * sat.same(a, b): a <-> b
 *   * sat.single(a): a = true
 *   * sat.diff(a, b): a xor b
 *   * sat.at_most_one(vec_num): <= 1 element in vec_num = true
 *  * Main idea of at_most_one: https://discuss.codechef.com/t/oh1dcare-editorial/99787
 * 
 * Status: stress-tested
 *   * At_most_one: https://www.codechef.com/viewsolution/1040190016
 *   * Normal: See on KATCL
 */
struct twosat_graph {
	int V = 0;
	std::vector<std::vector<int>> adj;

	twosat_graph() {}
	twosat_graph(const int& N) : V(N), adj(N << 1) {}

	inline void resize(const int& N) { V = N; adj.resize(N << 1); }
	inline void reserve(const int& N) { adj.reserve(N << 1); }
	inline void clear() { adj.clear(); V = 0; }

	inline int add_var() {
		adj.emplace_back(); adj.emplace_back();
		return V ++;
	}

	void either(int f, int j) {
		f = std::max(2 * f, - 1 - 2 * f);
		j = std::max(2 * j, - 1 - 2 * j);
		adj[f ^ 1].push_back(j);
		adj[j ^ 1].push_back(f);
	}
	void implies(int f, int j) {
		f = std::max(2 * f, - 1 - 2 * f);
		j = std::max(2 * j, - 1 - 2 * j);
		adj[f].push_back(j);
		adj[j ^ 1].push_back(f ^ 1);
	}
	void single(int f) {
		f = std::max(2 * f, - 1 - 2 * f);
		adj[f ^ 1].push_back(f);
		// either(f, f);
	}

	// NOTE: The value of them must be same to make the clause true
	// Similar to <-> operation
	void same(const int& f, const int& j) {
		implies(f, j); implies(~ f, ~ j);
	}

	// NOTE: The value of them must be different to make the clause true
	// Similar to xor operation
	void diff(const int& f, const int& j) {
		implies(f, ~ j); implies(~ f, j);
	}

	// NOTE: The new var will be added and returned
	// This var = true if (a == true or b == true)
	int add_or_var(const int &a, const int& b) {
		int res = add_var();
		implies(a, res); implies(b, res);
		return res;
	}

	// NOTE: The new var will be added and returned
	// This var = false if (a == false or b == false)
	int add_and_var(const int &a, const int& b) {
		int res = add_var();
		implies(res, a); implies(res, b);
		return res;
	}

	int at_most_one(const int& a, const int& b) {
		either(~ a, ~ b);
		return add_or_var(a, b);
	}

	template <typename Vec> int at_most_one(const Vec& vars) {
		if (vars.empty()) return - 1;
		int aux = vars[0];
		for (auto it = vars.begin() + 1; it != vars.end(); ++ it) {
			aux = at_most_one(aux, *it);
		}
		return aux;
	}

	std::vector<bool> solve() const {
		std::vector<int> c(V << 1, - 1), dfn(V << 1, - 1), low(V << 1, - 1);
		std::vector<int> stk; stk.reserve(V << 1);
		int now = 0, cnt_comp = 0;
		auto dfs = [&](auto&& dfs, const int& u) -> void {
			stk.push_back(u);
			dfn[u] = low[u] = now ++;
			for (const auto& v : adj[u]) {
				if (dfn[v] == - 1) {
					dfs(dfs, v);
					low[u] = std::min(low[u], low[v]);
				} else if (c[v] == - 1) {
					low[u] = std::min(low[u], dfn[v]);
				}
			}
			if (dfn[u] == low[u]) {
				int v;
				do {
					v = stk.back(); stk.pop_back();
					c[v] = cnt_comp;
				} while (v != u);
				++ cnt_comp;
			}
		};

		for (int i = 0; i < (V << 1); ++ i) {
			if (dfn[i] == - 1) dfs(dfs, i);
		}
		std::vector<bool> res(V);
		for (int i = 0; i < V; ++ i) {
			if (c[i << 1] == c[(i << 1) + 1]) return std::vector<bool>();
			res[i] = c[i << 1] > c[(i << 1) + 1];
		}
		return res;
	}

	std::vector<bool> operator () () const { return solve(); }
};