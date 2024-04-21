namespace eratosthenes{
	// Sieve of Eratosthenes: 0 <= x <= n
	vector<bool> sieve(const int& n){
		vector<bool> p(n + 1, true);
		p[0] = p[1] = false;
		for (int i = 2; i * i <= n; i ++){
			if (!p[i]) continue;
			for (int j = i * i; j <= n; j += i){
				p[j] = false;
			}
		}
		return p;
	}

	/**
	Range_sieve: l <= x <= r (r - l + 1 <= 1e7):
		- size(): get size of the range from l to r
		- build(l, r): sieve all prime numbers in range l to r and similar for constructor from (l, r)
		- is_prime(x): return true if x is prime number and false for the opposite case
		- count_prime(l, r): count how many prime numbers from l to r
		- all_prime(l, r): return the vector that contains all prime numbers from l to r
	
		NOTE:
			- This sieve is based on segmented sieve so that works fast as possible
			- To call count_prime and all_prime, you need to make sure that sieve has already built by build function or constructor
			- Beware about your range before using this range_sieve
	**/
	template <typename T = int64_t>
	struct range_sieve{
		T l, r;
		vector<bool> prime;
	
		int size(){
			return r - l + 1;
		}
	
		void build(T left, T right){  // depends on segmented sieve
			l = left; r = right;
			T lim = sqrt(r);
			vector<bool> marked(lim + 1, false);
			vector<ll> primes; // contains all prime number in range [2, sqrt(r)]
			for (ll i = 2; i <= lim; i ++){
				if (marked[i]) continue;
				primes.push_back(i);
				for (ll j = i * i; j <= lim; j += i){
					marked[j] = true;
				}
			}
			prime.assign(r - l + 1, true);
			for (ll i : primes){
				for (ll j = max(i * i, (l + i - 1) / i * i); j <= r; j += i){
					prime[j - l] = false;
				}
			}
			if (l == 1) prime[0] = false;
			if (l == 0) prime[0] = prime[1] = false;
		}
	
		// Constructor 1
		range_sieve(){}
	
		// Constructor 2
		range_sieve(T l, T r) : l(l), r(r){ // depends on segmented sieve
			T lim = sqrt(r);
			vector<bool> marked(lim + 1, false);
			vector<T> primes; // contains all prime number in range [2, sqrt(r)]
			for (T i = 2; i <= lim; i ++){
				if (marked[i]) continue;
				primes.push_back(i);
				for (T j = i * i; j <= lim; j += i){
					marked[j] = true;
				}
			}
			prime.assign(r - l + 1, true);
			for (T i : primes){
				for (T j = max(i * i, (l + i - 1) / i * i); j <= r; j += i){
					prime[j - l] = false;
				}
			}
			if (l == 1) prime[0] = false;
			if (l == 0) prime[0] = prime[1] = false;
		}
	
		// Check x is prime number
		inline bool is_prime(T x){
			assert(x >= l && x <= r);
			return prime[x - l];
		}
	
		// How many prime numbers from left to right
		int count_prime(T left, T right){
			assert(left >= l && right <= r);
			int cnt = 0;
			for (T i = left; i <= right; i ++){
				if (is_prime(i)) cnt ++;
			}
			return cnt;
		}
	
		vector<T> all_prime(T left, T right){
			assert(left >= l && right <= r);
			vector<T> all;
			for (T i = left; i <= right; i ++){
				if (is_prime(i)) all.push_back(i);
			}
			return all;
		}
	};

	// Linear sieve: 0 <= x <= n: return the minimum prime factor => usedful for number factorization
	// More memory so using when n <= 1e7
	vector<int> linear_sieve(const int& n){
		vector<int> lp(n + 1, 0);
		vector<int> pr;
		for (int i = 2; i <= n; i ++){
			if (lp[i] == 0){
				lp[i] = i;
				pr.push_back(i);
			}
			for (int j = 0; i * pr[j] <= n; j ++){ // Consider all numbers that divisible by i
				lp[i * pr[j]] = pr[j];
				if (pr[j] == lp[i]) break;
			}
		}
		return lp;
	}

	// Segmented sieve: time in O(nloglogn) and space O(sqrt(n))
	vector<int> segmented_sieve(const int& n){
		vector<pair<int, int>> primes;
		int nsqrt = sqrt(n);
		vector<char> isprime(nsqrt + 2, true);
		for (int i = 3; i < nsqrt; i += 2){
			if (!isprime[i]) continue;
			primes.push_back({i, (i * i - 1) / 2});
			for (int j = i * i; j <= nsqrt; j += 2 * i){
				isprime[j] = false;
			}
		}
	
		int S = nsqrt;
		vector<int> all = {2};
		vector<char> block(S);
		int high = (n - 1) / 2;
		for (int low = 0; low <= high; low += S){
			fill(block.begin(), block.end(), true);
			for (auto &i : primes){
				int p = i.first, idx = i.second;
				for (; idx < S; idx += p){
					block[idx] = false;
				}
				i.second = idx - S;
			}
			if (low == 0) block[0] = false;
			for (int i = 0; i < S && low + i <= high; i ++){
				if (block[i]){
					all.emplace_back((low + i) * 2 + 1);
				}
			}
		}
		return all;
	}

	// credit: min_25
	// takes 0.5s for n = 1e9
	vector<int> sieve_faster(const int N, const int Q = 17, const int L = 1 << 15) {
		static const int rs[] = {1, 7, 11, 13, 17, 19, 23, 29};
		struct P { 
			P(int p) : p(p) {}
			int p; int pos[8];
		};
		auto approx_prime_count = [] (const int N) -> int {
			return N > 60184 ? N / (log(N) - 1.1)
											 : max(1., N / (log(N) - 1.11)) + 1;
		};

		const int v = sqrt(N), vv = sqrt(v);
		vector<bool> isp(v + 1, true);
		for (int i = 2; i <= vv; ++i) if (isp[i]) {
			for (int j = i * i; j <= v; j += i) isp[j] = false;
		}

		const int rsize = approx_prime_count(N + 30);
		vector<int> primes = {2, 3, 5}; int psize = 3;
		primes.resize(rsize);

		vector<P> sprimes; size_t pbeg = 0;
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

		vector<unsigned char> pre(prod, 0xFF);
		for (size_t pi = 0; pi < pbeg; ++pi) {
			auto pp = sprimes[pi]; const int p = pp.p;
			for (int t = 0; t < 8; ++t) {
				const unsigned char m = ~(1 << t);
				for (int i = pp.pos[t]; i < prod; i += p) pre[i] &= m;
			}
		}

		const int block_size = (L + prod - 1) / prod * prod;
		vector<unsigned char> block(block_size); unsigned char* pblock = block.data();
		const int M = (N + 29) / 30;

		for (int beg = 0; beg < M; beg += block_size, pblock -= block_size) {
			int end = min(M, beg + block_size);
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
} // namespace eratosthenes