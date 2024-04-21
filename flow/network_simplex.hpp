/**
 * NETWORK SIMPLEX !!!
 * 
 * See: https://codeforces.com/blog/entry/94190
 * 
 * Tested: https://loj.ac/s/1990723
 * 
 * Details: * Network simplex is a method to solve by simplex method.
 *          * The result is based on min_cost-optimization
 *          * Supports edge lower bounds, negative costs and negative cost cycles
 * 
 * Runtime:
 *   * Expected runtime: O(VE) for positive costs, O(E²) for negative costs too
 *   * With positive costs, it really fast. Change to positive costs to solve (i.e by increasing)
 * 
 * Usage:
 *   * Init: network_simplex<flow_t, cost_t> simplex(num_node, num_edge)
 *       flow_t: should be large enough to hold node supply/excess and sum of all capacities
 *       cost_t: should be large enough to hold costs and potentials (usually >=64 bits)
 * 
 *   * Set supply for each node: simplex.nodes[u].supply = x;
 * 
 *   * Get excess for every vertex: excess(u) = flow(out of u) - flow(into u)
 *        auto excess = simplex.excess(); // std::vector<flow_t>
 * 
 *   * Solve the prob based on min_cost_optimization : simplex.simplex();
 * 
 *   * Check the prob can be solved: must call simplex() before check
 *        auto solvable = simplex.min_cost();
 *      Besides that, you can check sum_of_all_supply = 0 if you want
 * 
 *   * Get min_cost result: Iterate by yourself
 *        cost_t tot_cost = 0; for (edges : simplex) tot_cost += e_cost * e_flow;
 * 
 *   * Solve min_cost_max_flow with source s, sink t:
 *        simplex.nodes[s].supply += INF; simplex.nodes[t].supply -= INF;
 *        simplex.simplex(); // Solve it
 *        max_flow = simplex.min_cost_flow();
 *        Get min_cost: no change
 *      NOTE: The excess at a supply/source node u will be in the range [0,supply[u]].
 *            The excess at a demand/sink   node u will be in the range [supply[u],0].
**/
template <typename flow_t, typename cost_t> struct network_simplex {
	struct node_t {
		int par; // par node
		int pe; // par edge
		flow_t supply;
		cost_t pi; // potential of node
	};
	enum arc_state : int8_t {
		STATE_UPPER = - 1,
		STATE_TREE = 0,
		STATE_LOWER = 1
	};
	struct edge_t {
		int from;
		int to;
		flow_t lower;
		flow_t upper;
		cost_t cost;
		flow_t flow = 0;
		arc_state state = STATE_LOWER;
	};
	struct linked_list {
		int L, N;
		std::vector<int> next, prev;

		linked_list(const int& L = 0, const int& N = 0) { assign(L, N); }

		inline void clear() {
			std::iota(std::begin(next) + N, std::end(next), N);
			std::iota(std::begin(prev) + N, std::end(prev), N);
		}
		inline void assign(const int& L, const int& N) {
			this -> L = L; this -> N = N;
			next.resize(N + L), prev.resize(N + L), clear();
		}

		inline int rep(const int& l) const { return l + N; }
		inline int head(const int& l) const { return next[rep(l)]; }
		inline int tail(const int& l) const { return prev[rep(l)]; }
		inline void meet(const int& u, const int& v) { next[u] = v, prev[v] = u; }
		inline void meet(const int& u, const int& v, const int& w) { meet(u, v), meet(v, w); }

		inline void push_front(const int& l, const int& n) { meet(rep(l), n, head(l)); }
		inline void push_back(const int& l, const int& n) { meet(tail(l), n, rep(l)); }
		inline void erase(const int& n) { meet(prev[n], next[n]); }
	};

	int V = 0, E = 0;
	std::vector<node_t> nodes;
	std::vector<edge_t> edges;
	linked_list ch;

	int next_arc = 0, block_size = 0;
	std::vector<int> bfs, perm;

	network_simplex() {}
	network_simplex(const int& N) : V(N), nodes(N + 1) {}
	network_simplex(const int& N, const int& M) : V(N), nodes(N + 1) {
		edges.reserve(N + M); // Avoid allocation
	}

	inline int add_edge(const int& u, const int& v, const flow_t& lo, const flow_t& hi, const cost_t& cost) {
		return edges.emplace_back(edge_t{u, v, lo, hi, cost, STATE_LOWER}), E ++;
	}
	inline int add_node() { return nodes.emplace_back(), V ++; }

	inline cost_t reduced_cost(const int& e) const {
		return edges[e].cost + nodes[edges[e].from].pi - nodes[edges[e].to].pi;
	}
	inline cost_t signed_reduced_cost(const int& e) const { return edges[e].state * reduced_cost(e); }

	int select_pivot_edge() {
		// lemon-like block search, check block_size edges and pick the best one
		cost_t minimum = 0;
		int in_arc = -1;
		int count = block_size, seen_edges = E + V;
		for (int& e = next_arc; seen_edges -- > 0; e = e + 1 == E + V ? 0 : e + 1) {
			int x = perm[e];
			if (minimum > signed_reduced_cost(x)) {
				minimum = signed_reduced_cost(x);
				in_arc = x;
			}
			if (-- count == 0 && minimum < 0) {
				break;
			} else if (count == 0) {
				count = block_size;
			}
		}
		return in_arc;
	}
	void pivot(int in_arc) {
		// Find join node (lca of u_in and v_in)
		int u_in = edges[in_arc].from;
		int v_in = edges[in_arc].to;
		int a = u_in, b = v_in;
		while (a != b) {
			a = nodes[a].par == -1 ? v_in : nodes[a].par;
			b = nodes[b].par == -1 ? u_in : nodes[b].par;
		}
		int join = a;

		// Orient edge so that we add flow to source->target
		int source = edges[in_arc].state == STATE_LOWER ? u_in : v_in;
		int target = edges[in_arc].state == STATE_LOWER ? v_in : u_in;

		enum out_arc_side : int8_t {
			SAME_EDGE = 0,
			SOURCE_SIDE = 1,
			TARGET_SIDE = 2
		};
		flow_t flow_delta = edges[in_arc].upper;
		out_arc_side side = SAME_EDGE;
		int u_out = - 1;

		// Go up the cycle from source to the join node
		for (int u = source; u != join && flow_delta; u = nodes[u].par) {
			int e = nodes[u].pe;
			bool edge_down = u == edges[e].to;
			flow_t d = edge_down ? edges[e].upper - edges[e].flow : edges[e].flow;
			if (flow_delta > d) {
				flow_delta = d;
				u_out = u;
				side = SOURCE_SIDE;
			}
		}

		// Go up the cycle from target to the join node
		for (int u = target; u != join && (flow_delta || side != TARGET_SIDE); u = nodes[u].par) {
			int e = nodes[u].pe;
			bool edge_up = u == edges[e].from;
			flow_t d = edge_up ? edges[e].upper - edges[e].flow : edges[e].flow;
			if (flow_delta >= d) {
				flow_delta = d;
				u_out = u;
				side = TARGET_SIDE;
			}
		}

		// Augment along the cycle
		if (flow_delta) {
			auto delta = edges[in_arc].state * flow_delta;
			edges[in_arc].flow += delta;
			for (int u = edges[in_arc].from; u != join; u = nodes[u].par) {
				int e = nodes[u].pe;
				edges[e].flow += u == edges[e].from ? - delta : + delta;
			}
			for (int u = edges[in_arc].to; u != join; u = nodes[u].par) {
				int e = nodes[u].pe;
				edges[e].flow += u == edges[e].from ? + delta : - delta;
			}
		}

		if (side == SAME_EDGE) {
			edges[in_arc].state = arc_state(- edges[in_arc].state);
			return;
		}

		// Replace out_arc with in_arc in the spanning tree
		int out_arc = nodes[u_out].pe;
		edges[in_arc].state = STATE_TREE;
		edges[out_arc].state = edges[out_arc].flow ? STATE_UPPER : STATE_LOWER;

		// Put u_in on the same side as u_out
		u_in = side == SOURCE_SIDE ? source : target;
		v_in = side == SOURCE_SIDE ? target : source;

		// Walk up from u_in to u_out, then fix par/pred/child pointers backwards
		int i = 0, S = 0;
		for (int u = u_in; u != u_out; u = nodes[u].par) {
			bfs[S ++] = u;
		}
		for (i = S - 1; i >= 0; i--) {
			int u = bfs[i], p = nodes[u].par;
			ch.erase(p);
			ch.push_back(u, p);
			nodes[p].par = u;
			nodes[p].pe = nodes[u].pe;
		}
		ch.erase(u_in);
		ch.push_back(v_in, u_in);
		nodes[u_in].par = v_in;
		nodes[u_in].pe = in_arc;

		// Adjust potentials in the subtree of u_in (pi_delta is not 0)
		cost_t current_pi = reduced_cost(in_arc);
		cost_t pi_delta = u_in == edges[in_arc].from ? - current_pi : + current_pi;

		bfs[0] = u_in;
		for (i = 0, S = 1; i < S; i++) {
			int u = bfs[i];
			nodes[u].pi += pi_delta;
			for (int v = ch.head(u); v != ch.rep(u); v = ch.next[v]) {
				bfs[S ++] = v;
			}
		}
	}
	void simplex() {
		// Remove non-zero lower bounds and compute artif_cost as sum of all costs
		cost_t artif_cost = 1;
		for (int e = 0; e < E; ++ e) {
			const int& u = edges[e].from;
			const int& v = edges[e].to;
			edges[e].flow = 0;
			edges[e].state = STATE_LOWER;
			edges[e].upper -= edges[e].lower;
			nodes[u].supply -= edges[e].lower;
			nodes[v].supply += edges[e].lower;
			artif_cost += edges[e].cost < 0 ? - edges[e].cost : edges[e].cost;
		}

		edges.resize(E + V);
		bfs.resize(V + 1);
		ch.assign(V + 1, V + 1);

		// Add root <-> node artificial edges with initial supply for feasible flow
		int root = V;
		nodes[root] = node_t{- 1, - 1, 0, 0};

		for (int u = 0, e = E; u < V; ++ u, ++ e) {
			nodes[u].par = root, nodes[u].pe = e;
			ch.push_back(root, u);
			auto supply = nodes[u].supply;
			if (supply >= 0) {
				nodes[u].pi = - artif_cost;
				edges[e] = edge_t{u, root, 0, supply, artif_cost, supply, STATE_TREE};
			} else {
				nodes[u].pi = artif_cost;
				edges[e] = edge_t{root, u, 0, - supply, artif_cost, - supply, STATE_TREE};
			}
		}

		// We want to, hopefully, find a pivot edge in O(sqrt(E))
		// This should be < E to check different sets of edges in each select()
		// Otherwise we are vulnerable to simplex killers
		block_size = std::max(int(std::ceil(std::sqrt(E + V))), std::min(5, V + 1));
		next_arc = 0;

		// Random permutation of the edges; helps with wide graphs and killer test cases
		static std::mt19937 rng(random_device{}());
		perm.resize(E + V);
		std::iota(std::begin(perm), std::end(perm), 0);
		std::shuffle(std::begin(perm), std::end(perm), rng);

		// Pivot until we're done
		int in_arc = select_pivot_edge();
		while (in_arc != -1) {
			pivot(in_arc);
			in_arc = select_pivot_edge();
		}

		// Restore flows and supplies
		for (int e = 0; e < E; ++ e) {
			const int& u = edges[e].from;
			const int& v = edges[e].to;
			edges[e].flow += edges[e].lower;
			edges[e].upper += edges[e].lower;
			nodes[u].supply += edges[e].lower;
			nodes[v].supply -= edges[e].lower;
		}
	}

	bool min_cost() const {
		// TODO: Ensure that simplex() must be called and not change the edges.size()
		for (int e = E; e < E + V; ++ e) {
			if (edges[e].flow != 0) {
				// This prob can not be solved (INFEASIBLE)
				return false;
			}
		}
		return true;
	}
	flow_t min_cost_flow() const {
		flow_t tot_flow = 0;
		for (int e = E; e < E + V; ++ e) {
			if (edges[e].to == V) {
				tot_flow += edges[e].upper - edges[e].flow;
			}
		}
		return tot_flow;
	}
	std::vector<flow_t> excess() const {
		std::vector<flow_t> res(V);
		for (int e = 0; e < E; ++ e) {
			const int& u = edges[e].from;
			const int& v = edges[e].to;
			res[u] += edges[e].flow;
			res[v] -= edges[e].flow;
		}
		return res;
	}

	template <typename F>
	inline void for_each_edge(const int& L, const int& R, const F& f) {
		for (int e = L; e < R; ++ e) {
			f(edges[e]);
		}
	}
	template <typename F>
	inline void for_each_node(const int& L, const int& R, const F& f) {
		for (int e = L; e < R; ++ e) {
			f(nodes[e]);
		}
	}
};
