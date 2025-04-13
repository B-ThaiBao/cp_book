std::vector<bool> eratos_sieve(const int& N) {
	std::vector<bool> p(N, true);
	p[0] = p[1] = false;
	for (int i = 2; i * i < N; i++){
		if (!p[i]) continue;
		for (int j = i * i; j < N; j += i){
			p[j] = false;
		}
	}
	return p;
}

std::vector<int> block_sieve(const int K, const int Q = 17, const int L = 1 << 15) {
	if (K < 3) return {};
	if (K == 3) return {2};
	const int N = K - 1;
	static const int rs[] = {1, 7, 11, 13, 17, 19, 23, 29};
	struct P {
		P(int p_) : p(p_) {}
		int p; int pos[8];
	};
	auto tot_primes = [] (const int M) -> int {
		return M > 60184 ? int(M / (log(M) - 1.1)) : int(std::max(1., M / (log(M) - 1.11)) + 1);
	};

	const int v = int(std::sqrt(N)), vv = int(std::sqrt(v));
	std::vector<bool> isp(v + 1, true);
	for (int i = 2; i <= vv; ++i) if (isp[i]) {
		for (int j = i * i; j <= v; j += i) isp[j] = false;
	}

	const int rsize = tot_primes(N + 30);
	std::vector<int> primes = {2, 3, 5}; int psize = 3;
	primes.resize(rsize);

	std::vector<P> sprimes; size_t pbeg = 0;
	int pt = 1;
	for (int p = 7; p <= v; ++p) {
		if (!isp[p]) continue;
		if (p <= Q) pt *= p, ++pbeg, primes[psize++] = p;
		auto pp = P(p); 
		for (int t = 0; t < 8; ++t) {
			int j = (p <= Q) ? p : p * p;
			while (j % 30 != rs[t]) j += p << 1;
			pp.pos[t] = j / 30;
		}
		sprimes.push_back(pp);
	}

	std::vector<unsigned char> pre(pt, 0xFF);
	for (size_t pi = 0; pi < pbeg; ++pi) {
		auto pp = sprimes[pi]; const int p = pp.p;
		for (int t = 0; t < 8; ++t) {
			const uint8_t m = uint8_t(~(1 << t));
			for (int i = pp.pos[t]; i < pt; i += p) pre[i] &= m;
		}
	}

	const int block_size = (L + pt - 1) / pt * pt;
	std::vector<unsigned char> block(block_size); unsigned char* pblock = block.data();
	const int M = (N + 29) / 30;

	for (int beg = 0; beg < M; beg += block_size, pblock -= block_size) {
		int end = std::min(M, beg + block_size);
		for (int i = beg; i < end; i += pt) {
			copy(pre.begin(), pre.end(), pblock + i);
		}
		if (beg == 0) pblock[0] &= 0xFE;
		for (size_t pi = pbeg; pi < sprimes.size(); ++pi) {
			auto& pp = sprimes[pi];
			const int p = pp.p;
			for (int t = 0; t < 8; ++t) {
				int i = pp.pos[t]; const uint8_t m = uint8_t(~(1 << t));
				for (; i < end; i += p) pblock[i] &= m;
				pp.pos[t] = i;
			}
		}
		for (int i = beg; i < end; ++i) {
			for (int m = pblock[i]; m > 0; m &= m - 1) {
				primes[psize++] = i * 30 + rs[__builtin_ctz(m)];
			}
		}
	}
	assert(psize <= rsize);
	while (psize > 0 && primes[psize - 1] > N) --psize;
	primes.resize(psize);
	return primes;
}

struct linear_sieve : public std::vector<int> {
	std::vector<int> primes;
	inline void build(const int& N) {
		this->assign(N, 0); primes.reserve(N);
		for (int i = 2; i < N; ++i) {
			if (this->at(i) == 0) {
				this->at(i) = i;
				primes.emplace_back(i);
			}
			for (const auto& x : this->primes) {
				if (x > this->at(i) || i * x >= N) break;
				this->at(i * x) = x;
			}
		}
	}

	linear_sieve() {}
	linear_sieve(const int& N) { this->build(N); }
};

template <typename T, typename Q = T> struct range_sieve : public std::vector<bool> {
	T L;
	Q R;
	template <typename M, typename C> inline void build(const M& l, const C& r) {
		L = l, R = r - 1;
		Q sq = Q(std::sqrt(R));
		std::vector<bool> is_primes(sq + 1, true);
		std::vector<Q> primes;
		for (Q p = 2; p <= sq; ++p) {
			if (is_primes[p]) continue;
			primes.push_back(p);
			for (Q j = p * p; j <= sq; j += p) is_primes[j] = true;
		}
		this->assign(R - L + 1, true);
		for (const Q& p: primes) {
			for (Q j = std::max(p * p, Q(L + p - 1) / p * p); j <= R; j += p) {
				this->at(j - L) = false;
			}
		}
		++R;
		if (L == 0) this->at(0) = this->at(1) = false;
		if (L == 1) this->at(0) = false;
	}

	range_sieve() {}
	range_sieve(const int& l, const int& r) { this->build(l, r); }

	template <typename R> inline bool is_prime(const R& x) const { return this->at(x - L); }
};
