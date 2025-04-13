/**
 * WAVELET_TREE !!!
 *
 * Usage:
 *  - Create struct that contains bit_vector and 2 child nodes:
 *     + Example:
 *         struct wavelet_tree_node {
 *             std::array<int, 2> ch{-1, -1};
 *             bit_vector<T> bit;
 *             // Store more information like sum, max, min, .. for your needs
 *         };
 *         std::vector<wavelet_tree_node> nodes; nodes.reserve(N * log(A));
 *
 *     + NOTE: bit_vector<T> (bit[i] = 1 if val[i] is on the left node and 0 if val[i] is on the right node).
 *             Besides that, it also supports pref(i) is the prefix sum of bits in range [0, i).
 *             So that, it can change to satifiy your needs:
 *        * Without any update (change): use bit_vector<T> is good enough.
 *        * With update (change element): use binary_index_tree<int> or seg_tree instead.
 *        * Without delete element: use splay_tree or treap instead.
 *
 *  - Build wavelet_tree: build recursively and use std::stable_partition to split the array into two parts
 *     + Example:
 *         int build(auto begin, auto end, int mi, int ma) {
 *             if (begin >= end || mi >= ma) return -1;
 *             int idx = int(nodes.size());
 *             nodes.emplace_back();
 *             int md = (mi + ma) >> 1;
 *
 *             // Build prefix sum
 *             for (auto it = begin; it != end; it++) {
 *                 if (*it < md) // This is on the left  --> bit[i] = 1;
 *                 else          // This is on the right --> bit[i] = 0;
 *             }
 *
 *             if (ma - mi == 1) return idx;
 *
 *             // Split the array into two parts
 *             auto mit = std::stable_partition(begin, end, [&](auto x) { return x < md; });
 *
 *             nodes[idx].ch[0] = build(begin, mit, mi, md);
 *             nodes[idx].ch[1] = build(mit, end, md, ma);
 *             return idx;
 *         }
 *
 *  - For range [L, R):
 *     --> On the left node:  [pref(L), pref(R))
 *     --> On the right node: [L - pref(L), R - pref(R))
**/
struct bit_vector {
	int N;
	std::vector<std::pair<uint32_t, uint32_t>> data;

	inline void resize(const int& N_) {
		N = N_;
		// int sz = (N + 63) >> 5;
		data.resize((N >> 5) + 1, {0, 0});
	}
	inline size_t size() const { return size_t(N); }

	bit_vector() {}
	bit_vector(const int& N_) { this->resize(N_); }

	inline void set(const int& k, bool value = true) {
#ifdef _GLIBCXX_DEBUG
		assert(0 <= k && k < N);
#endif
		data[k >> 5].first &= ~(uint32_t(1) << (k & 31));
		data[k >> 5].first |= uint32_t(value) << (k & 31);
	}
	inline bool get(const int& k) const {
#ifdef _GLIBCXX_DEBUG
		assert(0 <= k && k < N);
#endif
		return (bool((data[k >> 5].first >> (k & 31)) & 1));
	}
	inline bool operator[](const int& k) const { return get(k); }

	inline void build() {
		// data[0].second = 0;
		for (int i = 0; i + 1 < int(data.size()); i++) {
			data[i + 1].second = data[i].second + __builtin_popcount(data[i].first);
		}
	}
	inline void clear() { data.clear(); }

	inline int prefix(const int& k) const {
#ifdef _GLIBCXX_DEBUG
		assert(0 <= k && k <= N);
#endif
		return data[k >> 5].second + __builtin_popcount(data[k >> 5].first & ((uint32_t(1) << (k & 31)) - 1));
	}
	inline int prefix(const int& k, const bool& val) const {
		return val ? prefix(k) : k - prefix(k);
	}
};
