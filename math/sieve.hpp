// Sieve nums in range [0, N]
std::vector<bool> eratos_sieve(const int& N) {
	std::vector<bool> is_prime(N + 1, true);
	is_prime[0] = is_prime[1] = false;
	for (int i = 2; i * i <= N; i ++){
		if (!is_prime[i]) continue;
		for (int j = i * i; j <= N; j += i){
			is_prime[j] = false;
		}
	}
	return is_prime;
}

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

// credit: min_25
// takes 0.5s for n = 1e9
std::vector<int> min25_sieve(const int N, const int Q = 17, const int L = 1 << 15) {
	static const int rs[] = {1, 7, 11, 13, 17, 19, 23, 29};
	struct P {
		P(int p) : p(p) {}
		int p; int pos[8];
	};
	auto approx_prime_count = [] (const int N) -> int {
		return N > 60184 ? N / (log(N) - 1.1) : std::max(1., N / (log(N) - 1.11)) + 1;
	};

	const int v = sqrt(N), vv = sqrt(v);
	std::vector<bool> isp(v + 1, true);
	for (int i = 2; i <= vv; ++i) if (isp[i]) {
		for (int j = i * i; j <= v; j += i) isp[j] = false;
	}

	const int rsize = approx_prime_count(N + 30);
	std::vector<int> primes = {2, 3, 5}; int psize = 3;
	primes.resize(rsize);

	std::vector<P> sprimes; size_t pbeg = 0;
	int prod = 1; 
	for (int p = 7; p <= v; ++p) {
		if (!isp[p]) continue;
		if (p <= Q) prod *= p, ++pbeg, primes[psize++] = p;
		auto pp = P(p); 
		for (int t = 0; t < 8; ++t) {
			int j = (p <= Q) ? p : p * p;
			while (j % 30 != rs[t]) j += p << 1;
			pp.pos[t] = j / 30;
		}
		sprimes.push_back(pp);
	}

	std::vector<unsigned char> pre(prod, 0xFF);
	for (size_t pi = 0; pi < pbeg; ++pi) {
		auto pp = sprimes[pi]; const int p = pp.p;
		for (int t = 0; t < 8; ++t) {
			const unsigned char m = ~(1 << t);
			for (int i = pp.pos[t]; i < prod; i += p) pre[i] &= m;
		}
	}

	const int block_size = (L + prod - 1) / prod * prod;
	std::vector<unsigned char> block(block_size); unsigned char* pblock = block.data();
	const int M = (N + 29) / 30;

	for (int beg = 0; beg < M; beg += block_size, pblock -= block_size) {
		int end = std::min(M, beg + block_size);
		for (int i = beg; i < end; i += prod) {
			copy(pre.begin(), pre.end(), pblock + i);
		}
		if (beg == 0) pblock[0] &= 0xFE;
		for (size_t pi = pbeg; pi < sprimes.size(); ++pi) {
			auto& pp = sprimes[pi];
			const int p = pp.p;
			for (int t = 0; t < 8; ++t) {
				int i = pp.pos[t]; const unsigned char m = ~(1 << t);
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
