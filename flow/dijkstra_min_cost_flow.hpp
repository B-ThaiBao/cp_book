#include <ext/pb_ds/priority_queue.hpp>

/**
 * NOTE: This doesn't support neg-cost edges, adjust it !!!
 * (by using spfa_potential or make it INF if know the start point and end point)
 * Besides that, this doesn't allow neg-cost cycle
 * 
 * Mostly inspired here:
 *   * https://www.geeksforgeeks.org/johnsons-algorithm/ (Idea)
 *   * https://codeforces.com/blog/entry/95823 (Theory)
 *   * https://jeffe.cs.illinois.edu/teaching/algorithms/notes/G-mincostflow.pdf (Paper)
 *   * https://codeforces.com/blog/entry/105658 (Theory)
 *   * https://codeforces.com/blog/entry/105658 (Theory)
 *   * https://usaco.guide/adv/min-cost-flow?lang=cpp (Theory)
 *   * https://codeforces.com/contest/1913/submission/237747038 (Implementation)
 *   * https://codeforces.com/contest/1913/submission/237769118 (Implementation)
 *   * https://codeforces.com/contest/1913/submission/237808589 (Implementation)
 *   * https://codeforces.com/contest/1913/submission/237768315 (Cancel cycle !!!)
 * 
 * Tested on: https://open.kattis.com/problems/mincostmaxflow
**/
template <typename flow_t, typename cost_t>
struct dijkstra_min_cost_flow {
	static constexpr cost_t INF_COST = std::numeric_limits<cost_t>::max() >> 1;
	static constexpr flow_t INF_FLOW = std::numeric_limits<flow_t>::max() >> 1;

	std::vector<cost_t> pi, dist;
	std::vector<int> pe;
	__gnu_pbds::priority_queue<std::pair<cost_t, int>> q;
	std::vector<typename decltype(q)::point_iterator> its;

	dijkstra_min_cost_flow() {}
	dijkstra_min_cost_flow(const int& N) : pi(N), dist(N), pe(N), its(N) {}

	// TODO: To calculate potential in the case that cost is neg
	// TRUE if using spfa, FALSE if no neg cost or can set pi by yourself
	bool using_spfa = false;
	template <typename G> void spfa_potential(const G& g, const int& s) {
		using_spfa = true;
		std::vector<int> deg(g.V, 0);
		for (int i = 0; i < g.V; i ++) {
			for (const auto& id : g.adj[i]) {
				const auto& e = g.edges[id];
				if (e.cap - e.flow > g.eps) {
					++ deg[e.to];
				}
			}
		}
		std::vector<int> q; q.reserve(g.V);
		for (int i = 0; i < g.V; i ++) {
			if (deg[i] == 0) q.push_back(i);
		}
		for (int beg = 0; beg < int(q.size()); beg ++) {
			for (const auto& id : g.adj[q[beg]]) {
				const auto& e = g.edges[id];
				if (e.cap - e.flow > g.eps) {
					-- deg[e.to];
					if (deg[e.to] == 0) {
						q.push_back(e.to);
					}
				}
			}
		}
		std::fill(pi.begin(), pi.end(), INF_COST);
		pi[s] = 0;
		if (int(q.size()) == g.V) {
			// TODO: This graph is DAG. We not need spfa to calculate potential
			// Just use topo_sort + DP to calculate potential in O(V + E)
			for (const auto& v : q) {
				if (pi[v] < INF_COST) {
					for (const auto& id : g.adj[v]) {
						const auto& e = g.edges[id];
						if (e.cap - e.flow > g.eps) {
							if (pi[v] + e.cost < pi[e.to]) {
								pi[e.to] = pi[v] + e.cost;
								pe[e.to] = id;
							}
						}
					}
				}
			}
		}
		else {
			// TODO: This graph have non-neg cycle. Using spfa
			q.assign(1, s);
			std::vector<bool> in_queue(g.V, false);
			in_queue[s] = true;
			for (int b = 0; b < int(q.size()); b ++) {
				int i = q[b];
				in_queue[b] = false;
				for (const auto& id : g.adj[i]) {
					const auto& e = g.edges[id];
					if (e.cap - e.flow > g.eps && pi[i] + e.cost < pi[e.to]) {
						pi[e.to] = pi[i] + e.cost;
						pe[e.to] = id;
						if (!in_queue[e.to]) {
							q.push_back(e.to);
							in_queue[e.to] = true;
						}
					}
				}
			}
		}
	}

	template <typename G> void dijkstra(const G& g, const int& s) {
		std::fill(dist.begin(), dist.end(), INF_COST);
		dist[s] = 0;
		std::fill(its.begin(), its.end(), q.end());
		its[s] = q.push({pi[s], s});
		while (!q.empty()) {
			auto i = q.top().second; q.pop();
			for (const auto& id : g.adj[i]) {
				const auto& e = g.edges[id];
				int j = e.to;
				cost_t nd = dist[i] + e.cost;
				if (e.cap - e.flow > g.eps && nd < dist[j]) {
					dist[j] = dist[i] + e.cost;
					pe[j] = id;
					if (its[j] == q.end()) {
						its[j] = q.push({pi[j] - dist[j], j});
					}
					else {
						q.modify(its[j], {pi[j] - dist[j], j});
					}
				}
			}
		}
		std::swap(dist, pi);
	}

	template <typename G>
	std::pair<flow_t, cost_t> min_cost_max_flow(G& g, const int& s, const int& t, const cost_t& max_cost = INF_COST - 1) {
		flow_t tot_flow = 0; cost_t tot_cost = 0;
		if (using_spfa) {
			// TODO: check that if t is reachable from s by using directed path
			assert(pi[t] <= max_cost);
			auto cur_flow = std::numeric_limits<flow_t>::max();
			int cur = t;
			while (cur != s) {
				const auto &e = g.edges[pe[cur]];
				cur_flow = std::min(cur_flow, e.cap - e.flow);
				cur = e.from;
			}
			cur = t;
			while (cur != s) {
				auto& e = g.edges[pe[cur]];
				e.flow += cur_flow;
				auto& back = g.edges[pe[cur] ^ 1];
				back.flow -= cur_flow;
				cur = e.from;
			}
			tot_flow += cur_flow;
			tot_cost += cur_flow * pi[t];
			using_spfa = false;
		}
		while (true) {
			dijkstra(g, s);
			if (pi[t] > max_cost) break;
			auto cur_flow = std::numeric_limits<flow_t>::max();
			int cur = t;
			while (cur != s) {
				const auto &e = g.edges[pe[cur]];
				cur_flow = std::min(cur_flow, e.cap - e.flow);
				cur = e.from;
			}
			cur = t;
			while (cur != s) {
				auto& e = g.edges[pe[cur]];
				e.flow += cur_flow;
				auto& back = g.edges[pe[cur] ^ 1];
				back.flow -= cur_flow;
				cur = e.from;
			}
			tot_flow += cur_flow;
			tot_cost += cur_flow * pi[t];
		}
		return std::make_pair(tot_flow, tot_cost);
	}

	template <typename G>
	std::pair<flow_t, cost_t> min_cost_flow(G& g, const int& s, const int& t, const flow_t& tot_goal, const cost_t& max_cost = INF_COST - 1) {
		flow_t tot_flow = 0; cost_t tot_cost = 0;
		if (using_spfa && tot_flow < tot_goal) {
			// TODO: check that if t is reachable from s by using directed path
			assert(pi[t] <= max_cost);
			auto cur_flow = tot_goal - tot_flow;
			int cur = t;
			while (cur != s) {
				const auto &e = g.edges[pe[cur]];
				cur_flow = std::min(cur_flow, e.cap - e.flow);
				cur = e.from;
			}
			cur = t;
			while (cur != s) {
				auto& e = g.edges[pe[cur]];
				e.flow += cur_flow;
				auto& back = g.edges[pe[cur] ^ 1];
				back.flow -= cur_flow;
				cur = e.from;
			}
			tot_flow += cur_flow;
			tot_cost += cur_flow * pi[t];
			using_spfa = false;
		}
		while (tot_flow < tot_goal) {
			dijkstra(g, s);
			if (pi[t] > max_cost) break;
			auto cur_flow = tot_goal - tot_flow;
			int cur = t;
			while (cur != s) {
				const auto &e = g.edges[pe[cur]];
				cur_flow = std::min(cur_flow, e.cap - e.flow);
				cur = e.from;
			}
			cur = t;
			while (cur != s) {
				auto& e = g.edges[pe[cur]];
				e.flow += cur_flow;
				auto& back = g.edges[pe[cur] ^ 1];
				back.flow -= cur_flow;
				cur = e.from;
			}
			tot_flow += cur_flow;
			tot_cost += cur_flow * pi[t];
		}
		return std::make_pair(tot_flow, tot_cost);
	}
};

template <typename flow_t, typename cost_t>
struct spfa_min_cost_flow {
	static constexpr cost_t INF_COST = std::numeric_limits<cost_t>::max() >> 1;
	static constexpr flow_t INF_FLOW = std::numeric_limits<flow_t>::max() >> 1;

	std::vector<cost_t> d;
	std::vector<int> q;  // Instead of std::queue
	std::vector<bool> in_queue;
	std::vector<int> par;

	spfa_min_cost_flow() {}
	spfa_min_cost_flow(const int& N) : d(N), in_queue(N), par(N) {}

	template <typename graph_t>
	bool spfa(graph_t& g, const int& s, const int& k) {
		std::fill(d.begin(), d.end(), INF_COST);
		q.clear();
		q.push_back(s);
		d[s] = 0;
		in_queue[s] = true;
		int beg = 0;
		bool found = false;
		while (beg < int(q.size())) {
			int i = q[beg ++];
			if (i == k) { found = true; } // Exist on path from s to k

			in_queue[i] = false;
			for (const auto& id : g.adj[i]) {
				const auto &e = g.edges[id];
				if (e.cap - e.flow > g.eps && d[i] + e.cost < d[e.to]) {
					d[e.to] = d[i] + e.cost;
					par[e.to] = id;
					if (!in_queue[e.to]) {
						q.push_back(e.to);
						in_queue[e.to] = true;
					}
				}
			}
		}
		return found;
	}

	template <typename graph_t>
	std::pair<flow_t, cost_t> min_cost(graph_t& g, const int& s, const int& t, const flow_t& tot_goal) {
		flow_t flow = 0;
		cost_t cost = 0;
		while (flow < tot_goal && spfa(g, s, t)) {
			auto cur_flow = tot_goal - flow;
			int v = t;
			while (v != s) {
				const auto& e = g.edges[par[v]];
				cur_flow = std::min(cur_flow, e.cap - e.flow);
				v = e.from;
			}

			v = t;
			while (v != s) {
				auto& e = g.edges[par[v]];
				e.flow += cur_flow;
				auto& back = g.edges[par[v] ^ 1];
				back.flow -= cur_flow;
				v = e.from;
			}

			flow += cur_flow;
			cost += cur_flow * d[t];
		}
		return std::make_pair(flow, cost);
	}

	template <typename graph_t>
	std::pair<flow_t, cost_t> min_cost_max_flow(graph_t& g, const int& s, const int& t) {
		flow_t flow = 0;
		cost_t cost = 0;
		while (spfa(g, s, t)) {
			auto cur_flow = std::numeric_limits<flow_t>::max();
			int v = t;
			while (v != s) {
				const auto& e = g.edges[par[v]];
				cur_flow = std::min(cur_flow, e.cap - e.flow);
				v = e.from;
			}

			v = t;
			while (v != s) {
				auto& e = g.edges[par[v]];
				e.flow += cur_flow;
				auto& back = g.edges[par[v] ^ 1];
				back.flow -= cur_flow;
				v = e.from;
			}

			flow += cur_flow;
			cost += cur_flow * d[t];
		}
		return std::make_pair(flow, cost);
	}
};