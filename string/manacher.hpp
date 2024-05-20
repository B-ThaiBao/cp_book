/**
 * Mostly inspired here: https://vnoi.info/wiki/algo/string/manacher.md
 *                       https://cp-algorithms.com/string/manacher.html
 *
 * Manacher(S): provides great algorithm for palindrome string
 *
 * For more details, we begin with radius of palindrome string (array name res):
 *   * Odd radius: radius centered at i-indexed ==> len of palindrome: 2 * res[i] + 1
 *   * Even radius: radius centered at space between point i and i + 1 ==> len of palindrome: 2 * res[i]
 *
 * For example, let consider this string s = "abaa" -> res = {0, 0, 1, 0, 0, 1, 0}
 *
 * Supports:
 *   * build(): input a vector or string to apply manacher algorithm (without restricted character)
 *   * [{x, r}]: x is index and r is remainder of type radius when modulo by 2
 *      * r = 0: return even radius between index x and x + 1 (similar to [2 * x + 1])
 *      * r = 1: return odd radius at x-indexed (similar to [2 * x])
**/
struct manacher : public std::vector<int> {
	template <typename String> inline void build(const String &S) {
		int N = int(S.size());
		if (N == 0) return;
		this->resize((N << 1) - 1);
		int l = -1, r = -1;
		for (int z = 0; z < ((N << 1) - 1); ++z) {
			int i = (z + 1) >> 1;
			int j = z >> 1;
			int p = (i >= r) ? 0 : std::min<int>(r - i, this->at(((l + r) << 1) - z));
			while (j + p + 1 < N && i - p - 1 >= 0){
				if (!(S[j + p + 1] == S[i - p - 1])) break;
				++p;
			}
			if (j + p > r){
				l = i - p; r = j + p;
			}
			this->at(z) = p;
		}
	}

	manacher() {}
	template <typename String> manacher(const String &S) { this->build(S); }

	inline int& operator [] (const std::array<int, 2>& x) {
		return this->at((x[0] << 1) | (x[1] ^ 1));
	}
	inline const int& operator [] (const std::array<int, 2>& x) const {
		return this->at((x[0] << 1) | (x[1] ^ 1));
	}

	// Check whether substring [l, r) is palindrome
	inline bool is_palindrome(const int& l, const int& r) const {
		if (l >= r) return false;
		auto p = l + r - 1;
		if (p & 1) return r - l <= (this->at(p) << 1);
		return r - l <= (this->at(p) << 1) + 1;
	}
};
