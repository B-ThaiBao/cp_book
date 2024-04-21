// NOTE: This doesn't support montgometry multiplication (quite slow if you don't apply it lazily)
// If time limit is strict, copy from here (if you want):
//            https://github.com/NyaanNyaan/library/blob/master/modint/montgomery-modint.hpp

template <typename T, const T num> struct constant {
	using value_type = T;
	static constexpr value_type value = num;

	template <typename U> friend U& operator << (U& i, const constant& c) { return i << c.value; }
	friend void __print(const constant& c) { std::cerr << c.value; }
};

template <typename T> struct inconstant {
	using value_type = T;
	static T value;

	template <typename U> friend U& operator << (U& i, const inconstant& c) { return i << c.value; }
	template <typename U> friend U& operator >> (U& i, inconstant& c) { return i >> c.value; }
	friend void __print(const inconstant& c) { std::cerr << c.value; }

	static void set_value(const T& t) { value = t; }
};

template <typename T> T inconstant<T>::value = 0;

template <typename T> struct naive_multiplier {
	template <typename U, typename P, typename Q, typename K = T>
	constexpr static typename std::enable_if<sizeof(K) == 4, T>::type multiply(const P& a, const Q& b, U md) {
		auto x = int64_t(a) * b;
		if (x < md) return x;
		return x % md;
	}
	template <typename U, typename P, typename Q, typename K = T>
	constexpr static typename std::enable_if<sizeof(K) == 8, T>::type multiply(const P& a, const Q& b, U md) {
		T q = T(static_cast<long double>(a) * b / md);
		auto x = a * b - q * md;
		if (x < md) return x;
		return x % md;
	}
};

template <typename T> struct barrett_multiplier {
	using barrett_t = typename std::conditional<sizeof(T) == 4, uint64_t, __uint128_t>::type;
	static T MOD;
	static barrett_t BARRETT_M;
	
	template <typename U, typename P, typename Q, typename K = T>
	constexpr static typename std::enable_if<sizeof(K) == 4, T>::type multiply(const P& a, const Q& b, U md) {
		if (md != MOD) {
			MOD = md;
			BARRETT_M = uint64_t(- 1) / md;
		}
		auto p = uint64_t(a) * b;
		uint64_t q = uint64_t((__uint128_t(BARRETT_M) * p) >> 64);
		q = p - q * MOD;
		if (q >= MOD) q -= MOD;
		return q;
	}
	template <typename U, typename P, typename Q, typename K = T>
	constexpr static typename std::enable_if<sizeof(K) == 8, P>::type multiply(const P& A, const Q& B, U md) {
		if (md != MOD) {
			MOD = md;
			BARRETT_M = __uint128_t(- 1) / md;
		}
		__uint128_t p = __uint128_t(A) * B;
		auto q = [](const __uint128_t& x, const __uint128_t& y) -> __uint128_t {
			// Multiply them and get the 128-highest bit in __uint256_t
			// Copied from CP-algorithm: https://cp-algorithms.com/algebra/montgomery_multiplication.html
			uint64_t a = x >> 64, b = x;
			uint64_t c = y >> 64, d = y;

			__uint128_t ac = __uint128_t(a) * c;
			__uint128_t ad = __uint128_t(a) * d;
			__uint128_t bc = __uint128_t(b) * c;
			__uint128_t bd = __uint128_t(b) * d;
			__uint128_t carry = (__uint128_t)(uint64_t)ad + (__uint128_t)(uint64_t)bc + (bd >> 64u);
			__uint128_t high = ac + (ad >> 64u) + (bc >> 64u) + (carry >> 64u);
			return high;
		}(BARRETT_M, p);
		q = p - q * MOD;
		if (q >= MOD) q -= MOD;
		return q;
	}
};

template <typename T> T barrett_multiplier<T>::MOD = 0;
template <typename T> typename barrett_multiplier<T>::barrett_t barrett_multiplier<T>::BARRETT_M = 0;

template <typename T, typename multiplier> struct modnum {
	using value_type = typename T::value_type;

	template <typename U>
	static value_type normalize(const U& x) {
		value_type v;
		if (- mod() <= x && x < mod()) v = static_cast<value_type>(x);
		else v = static_cast<value_type>(x % mod());
		if (v < 0) v += mod();
		return v;
	}
	template <typename K> K inverse(K a, K m) {
		K u = 0, v = 1;
		while (a != 0) {
			K t = m / a;
			m -= t * a; std::swap(a, m);
			u -= t * v; std::swap(u, v);
		}
		// assert(m == 1);
		return u;
	}
	const value_type& operator () () const { return value; }
	template <typename U>
	explicit operator U() const { return static_cast<U>(value); }
	constexpr static value_type mod() { return T::value; }
	template <typename U> friend U& operator << (U& out, const modnum& o) { return out << o.value; }
	template <typename U> friend U& operator >> (U& in, modnum& o) { int64_t v; in >> v; o.value = normalize(v); return in; }
	friend std::string to_string(const modnum& o) { return std::to_string(o.value); }
	friend bool is_zero(const modnum& o) { return o.value == value_type(0); }
	friend const value_type& abs(const modnum& x) { return x.value; }
	friend void __print(const modnum& o) { std::cerr << o.value; }

	value_type value;

	constexpr modnum() : value() {}
	template <typename U>
	modnum(const U& x) : value(normalize(x)) {}

	// This is comparision (Rare but still happen)
	friend bool operator == (const modnum& x, const modnum& y) { return x.value == y.value; }
	template <typename U> friend bool operator == (const modnum& x, const U& y) { return x == modnum(y); }
	template <typename U> friend bool operator == (const U& x, const modnum& y) { return modnum(x) == y; }

	friend bool operator != (const modnum& x, const modnum& y) { return x.value != y.value; }
	template <typename U> friend bool operator != (const modnum& x, const U& y) { return x != modnum(y); }
	template <typename U> friend bool operator != (const U& x, const modnum& y) { return modnum(x) != y; }

	friend bool operator < (const modnum& x, const modnum& y) { return x.value < y.value; }
	template <typename U> friend bool operator < (const modnum& x, const U& y) { return x < modnum(y); }
	template <typename U> friend bool operator < (const U& x, const modnum& y) { return modnum(x) < y; }

	friend bool operator <= (const modnum& x, const modnum& y) { return x.value <= y.value; }
	template <typename U> friend bool operator <= (const modnum& x, const U& y) { return x <= modnum(y); }
	template <typename U> friend bool operator <= (const U& x, const modnum& y) { return modnum(x) <= y; }

	friend bool operator > (const modnum& x, const modnum& y) { return x.value > y.value; }
	template <typename U> friend bool operator > (const modnum& x, const U& y) { return x > modnum(y); }
	template <typename U> friend bool operator > (const U& x, const modnum& y) { return modnum(x) > y; }

	friend bool operator >= (const modnum& x, const modnum& y) { return x.value >= y.value; }
	template <typename U> friend bool operator >= (const modnum& x, const U& y) { return x >= modnum(y); }
	template <typename U> friend bool operator >= (const U& x, const modnum& y) { return modnum(x) >= y; }

	modnum& operator += (const modnum& other) {
		value -= mod() - other.value;
		value = (value < 0) ? value + mod() : value;
		return *this;
	}
	modnum& operator -= (const modnum& other) {
		value -= other.value;
		value = (value < 0) ? value + mod() : value;
		return *this;
	}
	template <typename U> modnum& operator += (const U& other) { return *this += modnum(other); }
	template <typename U> modnum& operator -= (const U& other) { return *this -= modnum(other); }
	modnum& operator ++ () {
		++ value;
		if (value == mod()) value = 0;
		return *this;
	}
	modnum& operator -- () {
		if (value == 0) value = mod();
		-- value;
		return *this;
	}
	modnum operator ++ (int) { modnum result(*this); ++ (*this); return result; }
	modnum operator -- (int) { modnum result(*this); -- (*this); return result; }
	modnum operator - () const {
		modnum res;
		res.value = value ? mod() - value : 0;
		return res;
	}
	modnum operator + () const { return modnum(*this); }

	modnum& operator /= (const modnum& other) { return *this *= modnum(inverse(other.value, mod())); }
	modnum& operator *= (const modnum& o) {
		value = multiplier::multiply(value, o.value, (*this).mod()); // Can optimized this ???
		return *this;
	}

	friend modnum operator + (const modnum& x, const modnum& y) { return modnum(x) += y; }
	template <typename U> friend modnum operator + (const modnum& x, const U& y) { return modnum(x) += y; }
	template <typename U> friend modnum operator + (const U& x, const modnum& y) { return modnum(x) += y; }

	friend modnum operator - (const modnum& x, const modnum& y) { return modnum(x) -= y; }
	template <typename U> friend modnum operator - (const modnum& x, const U& y) { return modnum(x) -= y; }
	template <typename U> friend modnum operator - (const U& x, const modnum& y) { return modnum(x) -= y; }

	friend modnum operator * (const modnum& x, const modnum& y) { return modnum(x) *= y; }
	template <typename U> friend modnum operator * (const modnum& x, const U& y) { return modnum(x) *= y; }
	template <typename U> friend modnum operator * (const U& x, const modnum& y) { return modnum(x) *= y; }

	friend modnum operator / (const modnum& x, const modnum& y) { return modnum(x) /= y; }
	template <typename U> friend modnum operator / (const modnum& x, const U& y) { return modnum(x) /= y; }
	template <typename U> friend modnum operator / (const U& x, const modnum& y) { return modnum(x) /= y; }

	template <typename U> friend modnum bin_pow(const modnum& a, const U& b) {
		// assert(b >= 0);
		modnum x = a, res = 1;
		U p = b;
		while (p > 0) {
			if (p & 1) res *= x;
			x *= x;
			p >>= 1;
		}
		return res;
	}
};
