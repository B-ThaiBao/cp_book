/**
 * COMPRESSOR:
 * Provide machine to compress the number fix into the array (number can be smallest as possible)
 * 
 * Example: Compress: {1, 3, 6, 3} -> {0, 1, 2, 1}
 * For compress an array, you just use compress (it may be faster than this compressor)
**/
template <typename T> struct compressor {
	bool already_compress = false;
	std::vector<T> data;

	compressor(){}
	compressor(const int& N){
		// Beware of this one, may be it can be relocation if number of elements > N
		data.reserve(N); // Avoid relocation
	}

	inline void compress(){
		std::sort(data.begin(), data.end());
		data.erase(std::unique(data.begin(), data.end()), data.end());
		already_compress = true;
	}
	inline void insert(const T& x) { data.emplace_back(x); }

	inline int operator () (const T& x) {
		if (!already_compress) compress();
		return std::lower_bound(data.begin(), data.end(), x) - data.begin();
	}
	inline T& operator [] (const int& x) {
		if (!already_compress) compress();
		return data[x];
	}
	inline T& find_value(const int& x) {
		if (!already_compress) compress();
		return data[x];
	}

	inline int size() {
		if (!already_compress) compress();
		return static_cast<int>(data.size());
	}
	inline int max_element(){
		return size() - 1;
	}
}; // compressor

// NOTE: Using when you want to compress the numbers but don't want to change their order in O(nlogn)
// Source: https://codeforces.com/blog/entry/84164?#comment-716682
template <typename Vector> std::vector<int> compress(const Vector &A) {
	if (A.empty()) return {};
	int N = int(A.size());
	std::vector<int> res(N);
	std::vector<std::pair<decltype(A[0]), int>> pairs(N);
	for (int i = 0; i < N; ++ i){
		pairs[i] = {A[i], i};
	}
	std::sort(pairs.begin(), pairs.end());
	int nxt = 0;
	for (int i = 0; i < N; ++ i){
		if (i > 0 && pairs[i - 1].first != pairs[i].first) ++ nxt;
		res[pairs[i].second] = nxt;
	}
	return res;
}