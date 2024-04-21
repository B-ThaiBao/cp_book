/**
EULER TOTIENT FUNCTION:
	- phieuler(n): return the value of phi(n) in O(sqrt(n))
	- range_phieuler(n): return the value of phi(i) (1 <= i <= n) in O(nloglogn) (this idea based on sieve of Eratosthenes)

	DEFINITION:
		+ phi(n) count number from 1 to n, which are coprime with n (__gcd(n, x) == 1)

	PROPERTIES:
		+ phi(p) = p - 1 if p is a prime number
		+ phi(p ^ k) = p ^ k - p ^ (k - 1) if p is a prime number
		+ phi(a * b) = phi(a) * phi(b) if a and b is coprime (based on Chinese remainder theorem)
		+ General: phi(a * b) = phi(a) * phi(b) * d / phi(d) with d = gcd(a, b)
		==> FORMULA FOR PHI FUNCTION: phi(n) = n * (1 - 1 / p1) * (1 - 1 / p2) * ... * (1 - 1 / pk) with n = p1 ^ a1 * p2 ^ a2 * ... * pk ^ ak
		(p1, p2, ..., pk is a prime number)

	DIVISOR SUM PROPERTY:
		+ sum of all phi(d) = n with n % d == 0

	APPLICATIONS:
		+ a ^ phi(m) % m == 1 (a and m is coprime), this is generic more than Fermat's little theorem
			===> Can use to find INVERSE MODULO
		+ k = n % phi(m) ===> a ^ n % m == a ^ k (this just true with m is a prime number)
		===> Can use to find x ^ n mod m with n is very big value
		===> GENERALIZATION (if m is not a prime number):
			- x ^ n % m == x ^ (k + phi(m)) with x and m is not coprime
			- NOTE: This formula is just true for gcd(x, m) != 1 and n >= log2(m)

	SOURCE FOR THIS ONE: https://cp-algorithms.com/algebra/phi-function.html
**/
namespace phi_function{
	template <typename T>
	T phieuler(T n){
		T res = n;
		for (int i = 2; i * i <= n; i ++){
			if (n % i == 0){
				while (n % i == 0){
					n /= i;
				}
				res -= res / i;
			}
		}
		if (n > 1) res -= res / n;
		return res;
	}

	template <typename T>
	vector<T> range_phieuler(const T& n){
		vector<T> phi(n + 1);
		for (int i = 0; i <= n; ++ i) phi[i] = i;
		for (int i = 2; i <= n; ++ i){
			if (phi[i] == i){
				for (int j = i; j <= n; j += i){
					phi[j] -= phi[j] / i;
				}
			}
		}
		return phi;
	}

	template <typename T, typename Q =  int64_t>
	Q smallest_inverse_phieuler(const T &k){
		auto primes = sieve_faster((T) sqrt(k) + 1);
		int cnt_prime = primes.size();
		auto decomposition = y_combinator([&](auto decomposition, const T& n, const int& start, const T& pre) -> Q{
			if (n == 1) return 1;
			if (primetester.miller_rabin(n + 1) && n + 1 > pre) return n + 1;
			Q min_res = numeric_limits<Q> :: max();
			for (int i = start; i < cnt_prime; i ++){
				assert(primes[i] != 1);
				if (n % (primes[i] - 1)) continue;
				T cnt = 0;
				T tmp = n / (primes[i] - 1);
				Q res = decomposition(tmp, i + 1, primes[i]);
				if (res >= 0) min_res = min(min_res, (Q) primes[i] * res);
				while (tmp % primes[i] == 0){
					cnt ++;
					tmp /= primes[i];
					res = decomposition(tmp, i + 1, primes[i]);
					if (res >= 0) min_res = min(min_res, bin_pow((Q)(primes[i]), cnt + 1) * res);
				}
			}
			return (min_res == numeric_limits<Q> :: max()) ? - 1 : min_res;
		});
		Q res = decomposition(k, 0, - 1);
		return res;
	}
} // namespace euler_totient_function