// Copied from: https://judge.yosupo.jp/submission/189742
// NOTE: This supports for number that fixed into "uint64_t type" (Larger type doesn't allowed)
struct factorizer {
	static std::vector<int> min_prime;
	static std::vector<int> primes;

	static inline void eratos_sieve(const int& N) {
		min_prime.assign(N + 1, 0);
		for (int i = 2; i * i <= N; ++ i) {
			if (min_prime[i] == 0) {
				for (int j = i * i; j <= N; j += i) {
					if (min_prime[j] == 0) {
						min_prime[j] = i;
					}
				}
			}
		}
		primes.clear();
		for (int i = 2; i <= N; ++ i) {
			if (min_prime[i] == 0) {
				min_prime[i] = i;
				primes.emplace_back(i);
			}
		}
	}
	static inline void linear_sieve(const int& N) {
		min_prime.assign(N + 1, 0); primes.clear();
		for (int i = 2; i <= N; ++ i) {
			if (min_prime[i] == 0) {
				min_prime[i] = i;
				primes.emplace_back(i);
			}
			for (const auto& x : primes) {
				if (x > min_prime[i] || i * x > N) break;
				min_prime[i * x] = x;
			}
		}
	}

	template <typename T> constexpr static bool is_prime(const T& a) {
		if (a < 2) return false;
		int N = int(min_prime.size()) - 1;
		if (a <= N) return min_prime[a] == a;
		return miller_rabin::is_prime(a);
	}

	template <typename T> using factorize_t = std::vector<std::pair<T, int>>;
	template <typename T>
	constexpr static factorize_t<T> merge_factor(const factorize_t<T>& A, const factorize_t<T>& B) {
		size_t N = A.size(), M = B.size();
		factorize_t<T> C; C.reserve(N + M);
		size_t i = 0, j = 0;
		// NOTE: We try to use two pointers to merge them in the sorted factors
		while (i < N || j < M) {
			if (i < N && j < M && A[i].first == B[j].first) {
				C.emplace_back(A[i].first, A[i].second + B[j].second);
				++ i, ++ j;
				continue;
			}
			if (j == M || (i < N && A[i].first < B[j].first)) {
				C.emplace_back(A[i ++]);
			} else {
				C.emplace_back(B[j ++]);
			}
		}
		return C;
	}

	struct factorizer_mod {
		uint64_t mod;
		uint64_t mon_t;
		factorizer_mod(const uint64_t& N) : mod(N), mon_t(N) {
			for(int i = 0; i < 5; ++ i) {
			  mon_t *= 2 - mod * mon_t;
			}
		}

		constexpr uint64_t fast_mod(const uint64_t& a, const uint64_t& b, const uint64_t& c) const {
			const __uint128_t d = __uint128_t(a) * b;
			const uint64_t e = c + mod + (d >> 64);
			const uint64_t f = uint64_t(d) * mon_t;
			const uint64_t g = (__uint128_t(f) * mod) >> 64;
			return e - g;
		}

		constexpr uint64_t multiply(const uint64_t& a, const uint64_t& b) const {
			return fast_mod(a, b, 0);
		}
	};

	template <typename T>
	constexpr static factorize_t<T> pollard_rho(const T& N, const T& c) {
		if (N <= 1) return {};
		if ((N & 1) == 0){
			return merge_factor({{2, 1}}, pollard_rho(N >> 1, c));
		}
		if (is_prime(N)) return {{N, 1}};

		auto g = [](const T& Q) -> T {
			const factorizer_mod m(Q);

			constexpr uint64_t C_first = 1;
			constexpr uint64_t C_second = 2;
			constexpr uint64_t M = 512;

			uint64_t Z_first = 1;
			uint64_t Z_second = 2;
			try_find_divisor:
				uint64_t z_first = Z_first;
				uint64_t z_second = Z_second;
				for (size_t k = M;; k *= 2) {
					const uint64_t x_first = z_first + Q;
					const uint64_t x_second = z_second + Q;
					for (size_t j = 0; j < k; j += M) {
						const uint64_t y_first = z_first;
						const uint64_t y_second = z_second;

						uint64_t q_first = 1;
						uint64_t q_second = 2;
						z_first = m.fast_mod(z_first, z_first, C_first);
						z_second = m.fast_mod(z_second, z_second, C_second);
						for (size_t i = 0; i < M; ++ i) {
							const uint64_t t_first = x_first - z_first;
							const uint64_t t_second = x_second - z_second;
							z_first = m.fast_mod(z_first, z_first, C_first);
							z_second = m.fast_mod(z_second, z_second, C_second);
							q_first = m.multiply(q_first, t_first);
							q_second = m.multiply(q_second, t_second);
						}
						q_first = m.multiply(q_first, x_first - z_first);
						q_second = m.multiply(q_second, x_second - z_second);

						const uint64_t q_t = m.multiply(q_first, q_second);
						const uint64_t g_t = std::__gcd<uint64_t>(Q, q_t);
						if(g_t == 1) continue;
						if(g_t != Q) return g_t;

						const uint64_t g_first = std::__gcd<uint64_t>(Q, q_first);
						const uint64_t g_second = std::__gcd<uint64_t>(Q, q_second);

						const uint64_t C = g_first != 1 ? C_first : C_second;
						const uint64_t x = g_first != 1 ? x_first : x_second;
						uint64_t z = g_first != 1 ? y_first : y_second;
						uint64_t G = g_first != 1 ? g_first : g_second;

						if (G == Q) {
							do {
								z = m.fast_mod(z, z, C);
								G = std::__gcd<uint64_t>(Q, x - z);
							} while(G == 1);
						}
						if (G != Q) return G;

						Z_first += 2; Z_second += 2;
						goto try_find_divisor;
					}
				}
		}(N);
		return merge_factor(pollard_rho(g, c + 1), pollard_rho(N / g, c + 1));
	}

	template <typename T>
	constexpr static factorize_t<T> factorize(T x) {
		if (x <= 1) return {};
		int N = int(min_prime.size()) - 1;
		if (x <= N) {
			factorize_t<T> res;
			while (x > 1) {
				if (!res.empty() && res.back().first == min_prime[x]) {
					++ res.back().second;
				} else {
					res.emplace_back(min_prime[x], 1);
				}
				x /= min_prime[x];
			}
			return res;
		} 
		if (x <= int64_t(N) * N) {
			factorize_t<T> res;
			if (!is_prime(x)) {
				for (const auto& p : primes) {
					T t = x / p;
					if (p > t) break;
					if (x == t * p) {
						int cnt = 0;
						while (x % p == 0) {
							x /= p; ++ cnt;
						}
						res.emplace_back(p, cnt);
						if (is_prime(x)) break;
					}
				}
			}
			if (x > 1) res.emplace_back(x, 1);
			return res;
		}
		return pollard_rho(x, T(1));
	}

	template <typename T>
	constexpr static std::vector<T> find_divisor(const factorize_t<T>& factors, bool want_sort = false) {
		int num_divisors = 1;
		for (const auto& p : factors) num_divisors *= (p.second + 1);
		std::vector<T> divisors; divisors.reserve(num_divisors); divisors.emplace_back(1);
		for (const auto& p : factors) {
			size_t sz = divisors.size();
			for (size_t i = 0; i < sz; ++ i) {
				T cur = divisors[i];
				for (int j = 0; j < p.second; ++ j) {
					cur *= p.first;
					divisors.emplace_back(cur);
				}
			}
		}
		if (want_sort) std::sort(divisors.begin(), divisors.end());
		// assert(int(divisors.size()) == num_divisors);
		return divisors;
	}
};

std::vector<int> factorizer::min_prime = {0, 1};
std::vector<int> factorizer::primes;
