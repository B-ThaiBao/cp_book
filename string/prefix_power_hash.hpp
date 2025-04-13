/**
 * Mostly inspired here:
 *   * https://codeforces.com/contest/1909/submission/238567745 (Tourist)
 *   * https://codeforces.com/contest/825/submission/209664916 (Propane)
 * 
 * Usage:
 *   * Type T can be hashnum, pairnum<hashnum>, arraynum<hashnum, 3>
 *       -->>> Just set base one time, second must clear it
 *   * Random base: base_power<T>::base = mt();
 *   * Range hash: range_hash(frm, to)
 *   * Find lcp: find_lcp(hs, ss, ht, st)
 *   * Compare: compare(A, B) [- 1, 1]
**/
struct hashnum {
	using value_type = uint64_t;
	static constexpr uint64_t MOD = 18446744069414584321LLU;
	value_type v;

	static constexpr uint64_t modulo(const __uint128_t& y) {
		/**
		 * Copied from:
		 *  https://nyaannyaan.github.io/library/ntt/arbitrary-ntt-mod18446744069414584321.hpp
		**/
		uint64_t l = y & uint64_t(- 1);
		uint64_t m = (y >> 64) & unsigned(-1);
		uint64_t h = uint64_t(y >> 96);
		uint64_t u = h + m + (m ? MOD - (m << 32) : 0);
		uint64_t n = MOD <= l ? l - MOD : l;
		return n - u + (n < u ? MOD : 0);
	}

	constexpr hashnum() : v(0) {}
	template <typename U> hashnum(const U& x) {
		if (x >= 0 && x < MOD) v = static_cast<uint64_t>(x);
		else if (x > - MOD && x < 0) v = static_cast<uint64_t>(x + MOD);
		else if (x > MOD) v = modulo(x);
		else if (x < - MOD) v = (x % MOD) + MOD;
		else if (x == MOD || x == - MOD) v = 0;
	}
	const uint64_t& operator () () const { return v; }
	template <typename U>
	explicit operator U() const { return static_cast<U>(v); }
	constexpr static uint64_t mod() { return MOD; }
	template <typename U> friend U& operator << (U& out, const hashnum& o) { return out << o.v; }
	template <typename U> friend U& operator >> (U& in, hashnum& o) { int64_t n; in >> n; o = hashnum(n); return in; }
	friend void __print(const hashnum& o) { std::cerr << '{' << o.v << '}'; }
	friend std::string to_string(const hashnum& o) { return std::to_string(o.v); }
	friend bool is_zero(const hashnum& o) { return o.v == uint64_t(0); }
	friend const uint64_t& abs(const hashnum& x) { return x.v; }

	// This is comparision (Rare but still happen)
	friend bool operator == (const hashnum& x, const hashnum& y) { return x.v == y.v; }
	friend bool operator != (const hashnum& x, const hashnum& y) { return x.v != y.v; }
	friend bool operator < (const hashnum& x, const hashnum& y) { return x.v < y.v; }
	friend bool operator <= (const hashnum& x, const hashnum& y) { return x.v <= y.v; }
	friend bool operator > (const hashnum& x, const hashnum& y) { return x.v > y.v; }
	friend bool operator >= (const hashnum& x, const hashnum& y) { return x.v >= y.v; }
	hashnum& operator += (const hashnum& o) {
		if (v < MOD - o.v) v = v - (MOD - o.v) + MOD;
		else v -= (MOD - o.v);
		return *this;
	}
	hashnum& operator -= (const hashnum& o) {
		if (v < o.v) v = v - o.v + MOD;
		else v -= o.v;
		return *this;
	}
	hashnum& operator ++ () {
		++ v;
		if (v == MOD) v = 0;
		return *this;
	}
	hashnum& operator -- () {
		if (v == 0) v = MOD;
		-- v;
		return *this;
	}
	hashnum operator ++ (int) { hashnum r(*this); ++ (*this); return r; }
	hashnum operator -- (int) { hashnum r(*this); -- (*this); return r; }
	hashnum operator - () const {
		hashnum res;
		res.v = v ? MOD - v : 0;
		return res;
	}
	hashnum operator + () const { return hashnum(*this); }

	hashnum& operator /= (const hashnum& o) { return *this *= bin_pow(o, MOD - 2); }
	hashnum& operator *= (const hashnum& o) {
		v = modulo(__uint128_t(v) * o.v); // Can optimized this ???
		return *this;
	}

	friend hashnum operator + (const hashnum& x, const hashnum& y) { return hashnum(x) += y; }
	friend hashnum operator - (const hashnum& x, const hashnum& y) { return hashnum(x) -= y; }
	friend hashnum operator * (const hashnum& x, const hashnum& y) { return hashnum(x) *= y; }
	friend hashnum operator / (const hashnum& x, const hashnum& y) { return hashnum(x) /= y; }

	template <typename U> friend hashnum bin_pow(const hashnum& a, const U& b) {
		// assert(b >= 0);
		hashnum x = a, res = 1;
		U p = b;
		while (p > 0) {
			if (p & 1) res *= x;
			x *= x;
			p >>= 1;
		}
		return res;
	}
};

template <typename U, typename V = U> struct pairnum {
	U first;
	V second;

	pairnum() : first(), second() {}
	template <typename M> pairnum(const M& u) { first = U(u); second = V(u); }
	template <typename M, typename N>
	pairnum(const M& u, const N& v) { first = U(u); second = V(v); }
	template <typename M>
	explicit operator M() const { return M({first, second}); }

	template <typename S> friend S& operator << (S& out, const pairnum& p) { return out << p.first << ' ' << p.second; }
	template <typename S> friend S& operator >> (S& in, pairnum& p) { return in >> p.first >> p.second; }
	friend void __print(const pairnum& p) { std::cerr << '{' << p.first << ',' << p.second << '}'; }

	friend bool operator == (const pairnum& x, const pairnum& y) {
		return x.first == y.first && x.second == y.second;
	}
	template <typename S> friend bool operator == (const pairnum& x, const S& y) { return x == pairnum(y); }
	template <typename S> friend bool operator == (const S& x, const pairnum& y) { return pairnum(x) == y; }

	friend bool operator != (const pairnum& x, const pairnum& y) { return !(x == y); }
	template <typename S> friend bool operator != (const pairnum& x, const S& y) { return x != pairnum(y); }
	template <typename S> friend bool operator != (const S& x, const pairnum& y) { return pairnum(x) != y; }

	friend bool operator < (const pairnum& x, const pairnum& y) {
		return x.first != y.first ? x.first < y.first : x.second < y.second;
	}
	template <typename S> friend bool operator < (const pairnum& x, const S& y) { return x < pairnum(y); }
	template <typename S> friend bool operator < (const S& x, const pairnum& y) { return pairnum(x) < y; }

	friend bool operator <= (const pairnum& x, const pairnum& y) { return (x < y) || (x == y); }
	template <typename S> friend bool operator <= (const pairnum& x, const S& y) { return x <= pairnum(y); }
	template <typename S> friend bool operator <= (const S& x, const pairnum& y) { return pairnum(x) <= y; }

	friend bool operator > (const pairnum& x, const pairnum& y) {
		return x.first != y.first ? x.first > y.first : x.second > y.second;
	}
	template <typename S> friend bool operator > (const pairnum& x, const S& y) { return x > pairnum(y); }
	template <typename S> friend bool operator > (const S& x, const pairnum& y) { return pairnum(x) > y; }

	friend bool operator >= (const pairnum& x, const pairnum& y) { return (x > y) || (x == y); }
	template <typename S> friend bool operator >= (const pairnum& x, const S& y) { return x >= pairnum(y); }
	template <typename S> friend bool operator >= (const S& x, const pairnum& y) { return pairnum(x) >= y; }

	pairnum& operator += (const pairnum& o) {
		first += o.first; second += o.second; return *this;
	}
	pairnum& operator -= (const pairnum& o) {
		first -= o.first; second -= o.second; return *this;
	}
	template <typename S> pairnum& operator += (const S& other) { return *this += pairnum(other); }
	template <typename S> pairnum& operator -= (const S& other) { return *this -= pairnum(other); }
	pairnum& operator ++ () {
		++ first, ++ second; return *this;
	}
	pairnum& operator -- () {
		-- first, -- second; return *this;
	}
	pairnum operator ++ (int) { pairnum result(*this); ++ (*this); return result; }
	pairnum operator -- (int) { pairnum result(*this); -- (*this); return result; }
	pairnum operator - () const {
		return pairnum(- first, - second);
	}
	pairnum operator + () const { return pairnum(+ first, + second); }

	pairnum& operator /= (const pairnum& o) {
		first /= o.first; second /= o.second; return *this;
	}
	pairnum& operator *= (const pairnum& o) {
		first *= o.first; second *= o.second; return *this;
	}

	friend pairnum operator + (const pairnum& x, const pairnum& y) { return pairnum(x) += y; }
	template <typename S> friend pairnum operator + (const pairnum& x, const S& y) { return pairnum(x) += y; }
	template <typename S> friend pairnum operator + (const S& x, const pairnum& y) { return pairnum(x) += y; }

	friend pairnum operator - (const pairnum& x, const pairnum& y) { return pairnum(x) -= y; }
	template <typename S> friend pairnum operator - (const pairnum& x, const S& y) { return pairnum(x) -= y; }
	template <typename S> friend pairnum operator - (const S& x, const pairnum& y) { return pairnum(x) -= y; }

	friend pairnum operator * (const pairnum& x, const pairnum& y) { return pairnum(x) *= y; }
	template <typename S> friend pairnum operator * (const pairnum& x, const S& y) { return pairnum(x) *= y; }
	template <typename S> friend pairnum operator * (const S& x, const pairnum& y) { return pairnum(x) *= y; }

	friend pairnum operator / (const pairnum& x, const pairnum& y) { return pairnum(x) /= y; }
	template <typename S> friend pairnum operator / (const pairnum& x, const S& y) { return pairnum(x) /= y; }
	template <typename S> friend pairnum operator / (const S& x, const pairnum& y) { return pairnum(x) /= y; }
};

template <typename T> struct base_power {
	static T base;
	static std::vector<T> base_pow;

	static void ensure_power(const int& sz) {
		int cur = int(base_pow.size());
		if (cur - 1 < sz) {
			base_pow.resize(sz + 1);
			for (int i = cur; i <= sz; ++ i) {
				base_pow[i] = base_pow[i - 1] * base;
			}
		}
	}
	static void clear() {
		base_pow = std::vector<T>(1, T(1));
	}
	static const T& power(const int& x) { ensure_power(x); return base_pow[x]; }
	const T& operator () (const int& x) { ensure_power(x); return base_pow[x]; }
	const T& operator [] (const int& x) { ensure_power(x); return base_pow[x]; }
};
template <typename T> T base_power<T>::base;
template <typename T> std::vector<T> base_power<T>::base_pow = std::vector<T>(1, T(1));

template <typename T> struct prefix_power_hash {
	std::vector<T> data;

	template <typename String> static constexpr T hash(const String& S) {
		T res = T(1);
		for (auto pt = S.begin(); pt != S.end(); ++ pt) {
			res = res * (base_power<T>::base) + (*pt);
		}
		return res;
	}

	prefix_power_hash() {}
	template <typename String> prefix_power_hash(const String& S) { build(S); }

	template <typename String> inline void build(const String& S) {
		int N = int(S.size());
		base_power<T>::ensure_power(N + 1);
		data.resize(N + 1);
		data[0] = T(1);
		auto pt = S.begin();
		for (int i = 0; i < N; ++ i) {
			data[i + 1] = data[i] * (base_power<T>::base) + (*pt); ++ pt;
		}
	}

	inline void reserve(const int& N) { data.reserve(N); }
	constexpr size_t size() const { return data.size() - size_t(1); }
	inline T& operator [] (const int& x) { return data[x]; }
	inline const T& operator [] (const int& x) const { return data[x]; }
	inline T& back(const int& x) { return data.back(); }
	inline const T& back(const int& x) const { return data.back(); }

	template <typename Char> inline void push_back(const Char& c) {
		data.push_back(data.back() * (base_power<T>::base) + c);
	}
	template <typename Char> inline void emplace_back(const Char& c) {
		data.push_back(data.back() * (base_power<T>::base) + c);
	}

	// NOTE: [L, R)
	inline T range_hash(const int& frm, const int& to) const {
		return data[to] - (data[frm] * base_power<T>::power(to - frm));
	}
	inline T operator () (const int& frm, const int& to) const {
		return data[to] - (data[frm] * base_power<T>::power(to - frm));
	}
	inline T operator [] (const std::array<int, 2>& r) const {
		return data[r[1]] - (data[r[0]] * base_power<T>::power(r[1] - r[0]));
	}

	friend int find_lcp(const prefix_power_hash& A, const int& s, const prefix_power_hash& B, const int& t) {
		int mi = 0, ma = std::min<int>(A.size() - s, B.size() - t) + 1;
		while (ma - mi > 1) {
			int md = (ma + mi) >> 1;
			if (A.range_hash(s, s + md) == B.range_hash(t, t + md)) mi = md;
			else ma = md;
		}
		return mi;
	}
	friend int compare(const prefix_power_hash& A, const int& as, const int& ae, const prefix_power_hash& B, const int& bs, const int& be) {
		// Return - 1 (A < B), 0 (A == B), 1 (A > B)
		int mi = 0, ma = std::min<int>(ae - as, be - bs) + 1;
		while (ma - mi > 1) {
			int md = (ma + mi) >> 1;
			if (A.range_hash(as, as + md) == B.range_hash(bs, bs + md)) mi = md;
			else ma = md;
		}
		int pta = as + mi, ptb = bs + mi;
		if (pta == ae && ptb == be) return 0; // They are equal
		else if (pta == ae) return 1;
		else if (ptb == be) return - 1;
		return A.range_hash(pta, pta + 1) < B.range_hash(ptb, ptb + 1) ? - 1 : 1;
	}
};

// NOTE: Hash(A) = a, Hash(B) = b, size(b) = sz_b ---> Hash(A + B) = ???
template <typename hnum, typename sz_t>
static inline hnum merge(const hnum& a, const hnum& b, const sz_t& sz_b) {
	return a * sz_b + b;	
}
