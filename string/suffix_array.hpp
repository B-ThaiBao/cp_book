/**
 * SUFFIX_ARRAY using INDUCED_SORT (SAIS) !!!!
 * 
 * Main idea is from:
 *   * https://www.rahmannlab.de/lehre/alsa21/02-3-sais.pdf
 *   * https://web.stanford.edu/class/archive/cs/cs166/cs166.1196/lectures/04/Small04.pdf
 *   * https://zork.net/~st/jottings/sais.html
 *   * https://codeforces.com/blog/entry/91417
 *   * https://codeforces.com/edu/course/2/lesson/2
 * 
 * Copied from:
 *   * https://judge.yosupo.jp/problem/suffixarray
 *   * https://github.com/ecnerwala/cp-book/blob/master/src/suffix_array.hpp
 * 
 * Details:
 *   * Fatest algorithm as opposed to DC3, SAM or Radix_sort and only work in O(N)
 * 
 * Usage:
 *   * 1, Compress your array or your string as min as possible non-neg number
 *   * (e.g by using compressor or change the field of them)
 *   * Example:
 *   *              decrease by 'a'
 *   *   "acbade" ------------------> 0 2 1 0 3 4
 *   *
 *   * 2, Build a suffix array on the compressed array
 * 
 * Supports: prefix_array by reverse string and apply suffix array
 * 
 * Application: LCP array (by RMQ)
 * 
 * How to use it in constest: https://web.stanford.edu/class/cs97si/suffix-array.pdf
**/
struct suffix_array : public std::vector<int32_t> {
	using index_t = int32_t;
	template <typename String>
	static void sais(const int& N, const String& S, index_t* sa, int sigma, index_t* tmp) {
		if (N == 0) {
			sa[0] = 0;
			return;
		} else if (N == 1) {
			sa[0] = 1;
			sa[1] = 0;
			return;
		}

		// Phase 1: Initialize the frequency array, which will let us lookup buckets.
		index_t* freq = tmp; tmp += sigma;
		memset(freq, 0, sizeof(*freq) * sigma);
		for (int i = 0; i < N; ++ i) {
			++ freq[index_t(S[i])];
		}
		auto build_bucket_start = [&]() {
			int cur = 1;
			for (int v = 0; v < sigma; ++ v) {
				tmp[v] = cur;
				cur += freq[v];
			}
		};
		auto build_bucket_end = [&]() {
			int cur = 1;
			for (int v = 0; v < sigma; ++ v) {
				cur += freq[v];
				tmp[v] = cur;
			}
		};

		int num_pieces = 0;

		int first_endpoint = 0;
		// Phase 2: find the right-endpoints of the pieces
		{
			build_bucket_end();

			// Initialize the final endpoint out-of-band this way so that we don't try to look up tmp[-1].
			// This doesn't count towards num_pieces.
			sa[0] = N;

			index_t c0 = S[N - 1], c1 = - 1; bool isS = false;
			for (int i = N - 2; i >= 0; -- i) {
				c1 = c0;
				c0 = S[i];
				if (c0 < c1) {
					isS = true;
				} else if (c0 > c1 && isS) {
					isS = false;
					// insert i+1
					sa[first_endpoint = -- tmp[c1]] = i + 1;
					++ num_pieces;
				}
			}
		}

		// If num_pieces <= 1, we don't need to actually run the recursion, it's just sorted automatically
		// Otherwise, we're going to rebucket
		if (num_pieces > 1) {
			// Remove the first endpoint, we don't need to run the IS on this
			sa[first_endpoint] = 0;

			// Run IS for L-type
			{
				build_bucket_start();
				for (int z = 0; z <= N; z++) {
					int v = sa[z];
					if (!v) continue;

					// Leave for the S-round
					if (v < 0) continue;

					// clear out our garbage
					sa[z] = 0;

					-- v;
					index_t c0 = S[v - 1], c1 = S[v];
					sa[tmp[c1] ++] = (c0 < c1) ? ~ v : v;
				}
			}

			index_t* const sa_end = sa + N + 1;

			index_t* pieces = sa_end;
			// Run IS for S-type and compactify
			{
				build_bucket_end();
				for (int z = N; z >= 0; -- z) {
					int v = sa[z];
					if (!v) continue;

					// clear our garbage
					sa[z] = 0;

					if (v > 0) {
						* --pieces = v;
						continue;
					}

					v = ~ v;

					-- v;
					index_t c0 = S[v - 1], c1 = S[v];
					sa[-- tmp[c1]] = (c0 > c1) ? v : ~ v;
				}
			}

			// Compute the lengths of the pieces in preparation for equality
			// comparison, and store them in sa[v/2]. We set the length of the
			// final piece to 0; it compares unequal to everything because of
			// the sentinel.
			{
				int prv_start = N;
				index_t c0 = S[N-1], c1 = - 1; bool isS = false;
				for (int i = N - 2; i >= 0; -- i) {
					c1 = c0;
					c0 = S[i];
					if (c0 < c1) {
						isS = true;
					} else if (c0 > c1 && isS) {
						isS = false;

						// insert i+1
						int v = i + 1;
						sa[v >> 1] = prv_start == N ? 0 : prv_start - v;
						prv_start = v;
					}
				}
			}

			// Compute the alphabet, storing the result into sa[v/2].
			int next_sigma = 0;
			{
				int prv_len = -1, prv_v = 0;
				for (int i = 0; i < num_pieces; ++ i) {
					int v = pieces[i];
					int len = sa[v >> 1];

					bool eq = prv_len == len;
					for (int a = 0; eq && a < len; ++ a) {
						eq = S[v + a] == S[prv_v + a];
					}
					if (!eq) {
						++ next_sigma;
						prv_len = len;
						prv_v = v;
					}

					sa[v >> 1] = next_sigma; // purposely leave this 1 large to check != 0
				}
			}

			if (next_sigma == num_pieces) {
				sa[0] = N;
				memcpy(sa + 1, pieces, sizeof(*sa) * num_pieces);
			} else {
				index_t* next_S = sa_end;

				// Finally, pack the input to the SA
				{
					for (int i = (N - 1) >> 1; i >= 0; -- i) {
						int v = sa[i];
						if (v) *-- next_S = v - 1;
						sa[i] = 0;
					}
				}

				memset(sa, 0, sizeof(*sa) * (num_pieces+1));
				sais<const index_t*>(num_pieces, next_S, sa, next_sigma, tmp);

				{ // Compute the piece start points again and use those to map up the suffix array
					next_S = sa_end;
					index_t c0 = S[N-1], c1 = - 1; bool isS = false;
					for (int i = N-2; i >= 0; -- i) {
						c1 = c0;
						c0 = S[i];
						if (c0 < c1) {
							isS = true;
						} else if (c0 > c1 && isS) {
							isS = false;

							int v = i+1;
							*--next_S = v;
						}
					}
					sa[0] = N;
					for (int i = 1; i <= num_pieces; i++) {
						sa[i] = next_S[sa[i]];
					}
				}
			}

			// zero everything else
			memset(sa+num_pieces+1, 0, sizeof(*sa) * (N - num_pieces));

			{
				// Scatter the finished pieces
				build_bucket_end();
				for (int i = num_pieces; i > 0; -- i) {
					int v = sa[i];
					sa[i] = 0;

					index_t c1 = S[v];
					sa[-- tmp[c1]] = v;
				}
			}
		}

		// Home stretch! Just finish out with the L-type and then S-type
		{
			build_bucket_start();
			for (int z = 0; z <= N; ++ z) {
				int v = sa[z];
				if (v <= 0) continue;
				--v;
				index_t c1 = S[v];
				index_t c0 = v ? S[v - 1] : c1; // if v = 0, we don't want to invert
				sa[tmp[c1]++] = (c0 < c1) ? ~ v : v;
			}
		}

		// This just aggressively overwrites our original scattered pieces with the correct values
		{
			build_bucket_end();
			for (int z = N; z >= 0; -- z) {
				int v = sa[z];
				if (v >= 0) continue;
				sa[z] = v = ~v;
				-- v;
				index_t c1 = S[v];
				index_t c0 = v ? S[v-1] : c1+1;
				sa[-- tmp[c1]] = (c0 > c1) ? v : ~v;
			}
		}
	}

	template <typename String> void build(const String& S) {
		int N = int(S.size());
		(*this).resize(N + 1);
		// for (const auto& s : S) assert(index_t(s) >= 0);
		int sigma = N ? *std::max_element(S.begin(), S.end()) + 1 : 0;
		std::vector<index_t> tmp(std::max(N, sigma << 1));
		sais<String>(N, S, (*this).data(), sigma, tmp.data());
		(*this).erase((*this).begin());
	}
	suffix_array() {}
	template <typename String> suffix_array(const String& s) { build(s); }
};

// NOTE: This supports text search based on suffix array, time O(sigma * N)
// Source: https://codeforces.com/blog/entry/58991?#comment-425858
//         https://www.youtube.com/watch?v=X4kl3-UK_p0
struct burrows_wheeler_transform {
	int N = 0, sigma = 0;
	int last = 0;
	std::vector<int> loc;
	std::vector<int> fm_index;

	template <typename String, typename SuffixArray>
	void build(const String& S, const SuffixArray& sa) {
		N = int(S.size());
		sigma = N ? *std::max_element(S.begin(), S.end()) + 1 : 0;
		loc.assign(sigma + 1, 0);
		last = S.empty() ? - 1 : S.back();

		// N + 1 is the length of string s + '$'
		fm_index.assign(sigma * (N + 1), 0);
		for (int i = 0; i < N; ++ i) {
			++ loc[S[i] + 1];
			for (int j = 0; j < sigma; ++ j) {
				fm_index[(i + 1) * sigma + j] = fm_index[i * sigma + j];
			}
			if (sa[i] == 0) continue;
			++ fm_index[(i + 1) * sigma + S[sa[i] - 1]];
		}
		std::partial_sum(loc.begin(), loc.end(), loc.begin());
	}

	burrows_wheeler_transform() {}
	template <typename String, typename SuffixArray>
	burrows_wheeler_transform(const String& S, const SuffixArray& sa) { build(S, sa); }

	// NOTE: T must be compressed before pass thourgh function
	template <typename String>
	constexpr std::pair<int, int> find(const String& T, int lo = - 1, int hi = - 1) const {
		if (lo == - 1) lo = 0;
		if (hi == - 1) hi = int(T.size());

		int mi = 0, ma = N;
		bool x = false;
		for (int i = hi - 1; mi < ma && i >= lo; -- i) {
			const auto& c = T[i];
			mi = loc[c] + fm_index[mi * sigma + c] + (c == last && x);
			ma = loc[c] + fm_index[ma * sigma + c] + (c == last);
			x = true;
		}
		return std::make_pair(mi, ma);
	}
};

struct suffix_lcp_array {
	int N;
	// NOTE: lcp[i] = lcp(S[sa[i] ... N - 1], S[sa[i] + 1 ... N - 1]) in the order of suffix array
	std::vector<int32_t> lcp;
	// NOTE: rank[i]: suffix [i ... N - 1] is which rank in the order of suffix array
	std::vector<int32_t> rank;
	range_min_query<int> rmq;

	template <typename String, typename SuffixArray>
	void build(const String& s, const SuffixArray& sa, const bool& want_rmq = false) {
		N = int(s.size());
		rank.resize(N); lcp.resize(N - 1);
		for (int i = 0; i < N; ++ i) rank[sa[i]] = i;

		for (int i = 0, k = 0; i < N; ++ i) {
			if (rank[i] != N - 1) {
				int j = sa[rank[i] + 1];
				while (i + k < N && j + k < N && s[i + k] == s[j + k]) ++ k;
				lcp[rank[i]] = k;
				if (k > 0) -- k;
			} else {
				k = 0; continue;
			}
		}
		if (want_rmq) rmq.build(lcp);
	}

	suffix_lcp_array() {}
	template <typename String, typename SuffixArray>
	suffix_lcp_array(const String& s, const SuffixArray& sa, const bool& want_rmq = false) {
		build(s, sa, want_rmq);
	}

	inline int32_t find_lcp(const int& a) const {
		// lcp(a, a + 1) : Near one
		return lcp[rank[a]];
	}
	inline int32_t find_lcp(int a, int b) const {
		// assert(!rmq.data.empty());
		if (a == b) return N - a;
		a = rank[a], b = rank[b];
		if (a > b) std::swap(a, b);
		return rmq.range_query(a, b - 1).second;
	}

	inline int32_t operator () (const int& a) const { return find_lcp(a); }
	inline int32_t operator () (const int& a, const int& b) const { return find_lcp(a, b); }
	inline int32_t operator [] (const int& a) const { return find_lcp(a); }
	inline int32_t operator [] (const std::array<int, 2>& a) const { return find_lcp(a[0], a[1]); }

	uint64_t count_diff() const {
		// Idea: https://cp-algorithms.com/string/suffix-array.html#number-of-different-substrings
		uint64_t res = (uint64_t(N + 1) * N) >> 1;
		for (int i = 0; i < N - 1; ++ i) res -= lcp[i];
		return res;
	}
};