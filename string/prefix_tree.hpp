/**
 * PREFIX_TREE: Code was made by thaibao with lots of sources
 *
 * Usage:
 *   * Init the tree and call make_root: prefix_tree p(num_childs); p.make_root();
 *   * Insert: p.insert(string);
 *   * Erase: p.erase(string);
 *   * Find: p.find(string);
 *   * Count prefix: p.count_prefix(string, is_full);
 *   * Get next node: p[node][char];
 *
 * NOTE: You can store more information for each node by outside std::vector that you like.
**/
struct prefix_tree_node : std::vector<int> {
	int num_starts = 0; // number of strings that has the prefix is the string from root to this node
	int num_ends = 0;   // number of strings that fit with string from root to this node

	prefix_tree_node() {};
	prefix_tree_node(const int& num_chs) : std::vector<int>(num_chs, -1) {}
};

struct prefix_tree : std::vector<prefix_tree_node> {
	int num_childs = 0;

	prefix_tree() = default;
	prefix_tree(const int& num_chs) : num_childs(num_chs) {}

	// NOTE: Make sure to call make_root() before do_ops with prefix_tree
	inline void make_root() { this->emplace_back(this->num_childs); }
	static constexpr int root() { return 0; }
	inline int next(const int& p, const int& c) {
		if (this->at(p).at(c) < 0) {
			this->at(p).at(c) = int(this->size());
			this->emplace_back(this->num_childs);
		}
		return this->at(p).at(c);
	}
	template <typename String> inline int build(const String& S, const int& d) {
		int r = 0;
		for (const auto& ch : S) {
			this->at(r).num_starts += d;
			r = this->next(r, ch);
		}
		this->at(r).num_starts += d;
		return r;
	}

	template <typename String> int insert(const String& S) {
		int r = build(S, 1);
		this->at(r).num_ends++;
		return r;
	}
	template <typename String> int erase(const String& S) {
		int r = build(S, -1);
		this->at(r).num_ends--;
		return r;
	}
	template <typename String> int find(const String& S) const {
		int r = 0;
		for (const auto& ch : S) {
			r = this->at(r).at(ch);
			if (r < 0) break;
		}
		return r;
	}

	// Given a string, how many words in the trie are prefixes of the string?
	template <typename String> int count_prefix(const String& S, bool full = true) const {
		int r = 0, c = 0;
		for (const auto& ch : S) {
			c += this->at(r).num_ends;
			r = this->at(r).at(ch);
			if (r < 0) break;
		}
		if (full && r >= 0) c += this->at(r).num_ends;
		return c;
	}
};
