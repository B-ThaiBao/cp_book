/**
 * COMPRESSOR:
 * Provide machine to compress the number fix into the array (number can be smallest as possible)
 *
 * Example: Compress: {1, 3, 6, 3} -> {0, 1, 2, 1}
 * For compress an array, you just use compress (it may be faster than this compressor)
 *
 * NOTE: Please sort before compress() (this function doesn't sort directly)
**/
template <typename T> struct compressor : public std::vector<T> {
	using std::vector<T>::vector;

	inline void compress() {
		this->erase(std::unique(this->begin(), this->end()), this->end());
	}
	inline int operator () (const T &x) const {
		return int(std::lower_bound(this->begin(), this->end(), x) - this->begin());
	}
};

// NOTE: Using when you want to compress the numbers but don't want to change their order in O(nlogn)
// Source: https://codeforces.com/blog/entry/84164?#comment-716682
template <typename Vector> std::vector<int> compress(const Vector &A, int nxt = 0) {
	if (A.empty()) return {};
	int N = int(A.size());
	std::vector<int> res(N);
	std::vector<std::pair<typename std::decay<decltype(A[0])>::type, int>> pairs(N);
	for (int i = 0; i < N; ++i) pairs[i] = {A[i], i};
	std::sort(pairs.begin(), pairs.end());
	for (int i = 0; i < N; ++ i) {
		if (i > 0 && pairs[i - 1].first != pairs[i].first) ++nxt;
		res[pairs[i].second] = nxt;
	}
	return res;
}
