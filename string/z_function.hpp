struct z_function : public std::vector<int> {
	template <typename String> inline void build(const String &S) {
		int N = int(S.size());
		this->assign(N, 0);
		for (int i = 1, l = 0, r = 0; i < N; i++) {
			this->at(i) = (i > r ? 0 : std ::min(r - i + 1, this->at(i - l)));
			while (i + this->at(i) < N && S[this->at(i)] == S[i + this->at(i)])
				this->at(i)++;
			if (i + this->at(i) - 1 > r) {
				l = i;
				r = i + this->at(i) - 1;
			}
		}
	}

	z_function() {}
	template <typename String> z_function(const String &S) { this->build(S); }

	constexpr int compress(int N = -1) const {
		if (N == - 1) N = int(this->size());
		for (int i = 1; i < N; i++) {
			if (N % i == 0 && i + this->at(i) == N) return i;
		}
		return N;
	}

	inline std::vector<int> count_prefix() const {
		int N = int(this->size());
		std::vector<int> res(N + 1);
		for (int i = 1; i < N; i++) {
			if (this->at(i) >= 1) ++res[this->at(i) - 1];
		}
		res.back()++;
		for (int i = N - 1; i >= 0; --i) {
			res[i] += res[i + 1];
		}
		res.pop_back();
		return res;
	}
};
