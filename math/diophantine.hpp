// NOTE: Solve the linear equation: ax + by = gcd(a, b)
// Pass x, y are the parameters and return gcd(a, b)
template <typename T> T extended_gcd(const T& a, const T& b, T& ay, T& by) {
	T x = a, y = b;
	T ax = 1, bx = 0;
	ay = 0, by = 1;
	while (x != 0) {
		T k = y / x;
		y -= k * x;
		ay -= k * ax;
		by -= k * bx;
		std::swap(x, y);
		std::swap(ax, ay);
		std::swap(bx, by);
	}
	return y;
}

// NOTE: Solve more general linear equation: ax + by = c (avoid overflow)
// Pass x, y, g are the parameters and return if equaltion has feasible solution
// All roots: x + k * (b / g) and y - k * (a / g) with abitary k
template <typename T> bool diophantine(T a, T b, T c, T& x, T& y, T& g) {
	if (a == 0 && b == 0) {
		if (c == 0) {
			x = y = g = 0;
			return true;
		}
		return false;
	}
	auto num_abs = [](const auto& x) { return x < 0 ? (- x) : x; };
	if (a == 0) {
		if (c % b == 0) {
			x = 0, y = c / b, g = num_abs(b);
			return true;
		}
		return false;
	}
	if (b == 0) {
		if (c % a == 0) {
			x = c / a, y = 0, g = num_abs(a);
			return true;
		}
		return false;
	}
	g = extended_gcd(a, b, x, y);
	if (c % g != 0) return false;
	T dx = c / a;
	c -= dx * a;
	T dy = c / b;
	c -= dy * b;
	x = dx + (T) ((__int128_t) x * (c / g) % b);
	y = dy + (T) ((__int128_t) y * (c / g) % a);
	g = num_abs(g);
	return true;
	// NOTE: |x|, |y| <= max(|a|, |b|, |c|) [tested]
}

// NOTE: a mod mx = kx, a mod my = ky
// The result is that a mod m = k, return true if exist m and k like that
template <typename T, typename R>
bool crt_diophantine(const T& mx, R kx, const T& my, R ky, T& m, R& k) {
	if (kx > mx) kx %= mx;
	if (kx < 0) kx %= mx, kx += mx;
	if (ky > my) ky %= my;
	if (ky < 0) ky %= my, ky += my;
	T x, y, g;
	if (!diophantine(mx, - my, ky - kx, x, y, g)) return false;
	T dx = my / g;
	T delta = x / dx - (x % dx < 0);
	k = mx * (x - dx * delta) + kx;
	m = mx / g * my;
	// assert(0 <= k && k < m);
	return true;
}

// Idea: https://cp-algorithms.com/algebra/garners-algorithm.html
// NOTE: Just apply for distinct prime modulos
// Time: O(N ^ 2)
template <typename P, typename Values, typename T>
P crt_garner(const P& p, const Values& a, T& res) {
	// assert(p.size() == a.size());
	if (p.empty()) { res = 0; return {}; }
	using value_t = typename std::decay<decltype(p[0])>::type;
	auto garner_inv = [](auto q, const auto& m) {
		if (q > m) q %= m;
		if (q < 0) q %= m, q += m;
		value_t b = m, u = 0, v = 1;
		while (q != 0) {
			value_t t = b / q;
			b -= t * q; std::swap(q, b);
			u -= t * v; std::swap(u, v);
		}
		// assert(b == 1);
		if (u < 0) u += m;
		return u;
	};

	using product_value_t = typename std::conditional<
		sizeof(value_t) < sizeof(int64_t),
		int64_t,
		__int128_t
	>::type;

	P x(p.size());
	for (int i = 0; i < int(p.size()); ++ i) {
		// assert(0 <= a[i] && a[i] < p[i]);
		x[i] = a[i];
		for (int j = 0; j < i; ++ j) {
			x[i] = (value_t) (product_value_t(x[i] - x[j]) * garner_inv(p[j], p[i]) % p[i]);
			if (x[i] < 0) x[i] += p[i];
		}
	}
	res = 0;
	for (int i = int(p.size()) - 1; i >= 0; -- i) {
		res = res * p[i] + x[i];
	}
	return x;
}