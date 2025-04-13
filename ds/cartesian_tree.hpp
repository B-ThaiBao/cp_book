/**
 * CARTESIAN TREE !!
 *
 * Return par of nodes in cartesian tree by principles:
 *   * Leftmost min tree : std::less<>
 *   * Rightmost min tree: std::less_equal<>
 *   * Leftmost max tree : std::greater<>
 *   * Rightmost max tree: std::greater_equal<>
**/
template <typename Vector, typename Comp = std::less<>>
static inline std::vector<int> build_cartesian_tree(const Vector& A, Comp comp = Comp()) {
	int N = int(A.size());
	// We maintain monotonic stack that top of stack is the par of cur idx
	std::vector<int> stk; stk.reserve(N);
	std::vector<int> par(N, -1);
	for (int i = 0; i < N; ++i) {
		int lst = -1;
		while (!stk.empty() && comp(A[i], A[stk.back()])) {
			lst = stk.back(); stk.pop_back();
		}
		if (lst != -1) par[lst] = i;
		if (!stk.empty()) par[i] = stk.back();
		stk.push_back(i);
	}
	return par;
}
