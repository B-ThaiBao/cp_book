/**
 * FAST FOURIER TRANSFORM (FFT) !!!
 *
 * Explanation : http://neerc.ifmo.ru/trains/toulouse/2017/fft2.pdf
 * Code source : https://nyaannyaan.github.io/library
 * Applications: https://www.informatika.bg/resources/StringMatchingWithFFT.pdf
 *
 * Usage:
 *  Everything has been wrapped in several following classes:
 *   * naive_multiplier<T>
 *   * fft_complex_multiplier<T>
 *   * fft_complex_double_multiplier<T>
 *   * fft_mod_multiplier<T>
 *   * fft_numeric_multiplier<T> (NTT) !!!
 *   * multiply_inverser<T, multiplier>
 *   * fft_numeric_inverser<T>
 *
 *  Also, please call it though the following functions:
 *   * naive_multiply(a, b)
 *   * fft_complex_multiply(a, b)
 *   * fft_complex_double_multiply(a, b)
 *   * fft_complex_mod_multiply(a, b)
 *   * fft_numeric_multiply(a, b)
 *   * naive_inverse(a)
 *   * fft_complex_inverse(a)
 *   * fft_complex_double_inverse(a)
 *   * fft_complex_mod_inverse(a)
 *   * fft_numeric_inverse(a)
**/
namespace fft {

static inline int ceil_log_2(const int& N) { return N > 4 ? 32 - __builtin_clz(N - 1) : 2; }

template <typename num> struct naive_multiplier {
	template <typename U, typename V, typename P>
	static inline void multiply(const U& a, const int& sza, const V& b, const int& szb, P& r) {
		for (int i = 0; i < sza; i++) {
			for (int j = 0; j < szb; j++) {
				r[i + j] += num(a[i]) * num(b[j]);
			}
		}
	}
};

template <typename num> struct fft {
	// TODO: Wrap variables of convolution in static
	// fft for the purpose of memory optimization
	static std::vector<num> scratch_a;
	static std::vector<num> scratch_b;
};
template <typename num> std::vector<num> fft<num>::scratch_a;
template <typename num> std::vector<num> fft<num>::scratch_b;

template <typename num> struct complex {
	num x, y;
	constexpr complex() : x(0), y(0) {}
	constexpr complex(num _x, num _y) : x(_x), y(_y) {}
	constexpr inline complex operator+(const complex& c) const {
		return complex(x + c.x, y + c.y);
	}
	constexpr inline complex operator-(const complex& c) const {
		return complex(x - c.x, y - c.y);
	}
	constexpr inline complex operator*(const complex& c) const {
		return complex(x * c.x - y * c.y, x * c.y + y * c.x);
	}
	constexpr inline complex operator-() const { return complex(-x, -y); }
	constexpr inline complex conj() const { return complex(x, -y); }
	constexpr inline complex prep_ccw() const { return complex(-y, x); }
	template <typename U> friend U& operator<<(U& os, const complex& c) {
		os << "(" << c.x << ", " << c.y << ")";
		return os;
	}
};

template <typename num> struct fft_complex {
	// Alias for convienience
	using cnum = complex<double>;

	static std::vector<cnum> roots;

	static void init(int k) {
		static long double PI = acosl(-1);
		--k;
		if (int(roots.size()) >= (1 << k)) return;
		roots.resize(1 << k);
		std::vector<complex<long double>> base(k);
		const long double arg = PI / (1 << k);
		for (int i = 0, j = 1 << (k - 1); j; i++, j >>= 1) {
			const long double p = arg * j;
			base[i] = complex<long double>(cosl(p), sinl(p));
		}
		make_root(0, k - 1, complex<long double>{1, 0}, base);
	}
	static void make_root(int i, int b, complex<long double> z, const std::vector<complex<long double>>& base) {
		if (b == -1) {
			roots[i].x = double(z.x), roots[i].y = double(z.y);
		} else {
			make_root(i, b - 1, z, base);
			make_root(i | (1 << b), b - 1, z * base[b], base);
		}
	}

	template <typename Vector> static void fft(Vector& A, int k) {
		if (k <= 0) return;
		if (k == 1) {
			cnum v = A[1];
			A[1] = A[0] - A[1];
			A[0] = A[0] + v;
			return;
		}
		if (k & 1) {
			int v = 1 << (k - 1);
			for (int j = 0; j < v; ++j) {
				cnum a = A[j + v];
				A[j + v] = A[j] - a;
				A[j] = A[j] + a;
			}
		}
		int u = 1 << (k & 1), v = 1 << (k - 2 - (k & 1));
		while (v) {
			{
				int ja = 0;
				int jb = v;
				int jc = jb + v;
				int jd = jc + v;
				int je = v;
				for (; ja < je; ++ja, ++jb, ++jc, ++jd) {
					cnum ta = A[ja], tb = A[jb], tc = A[jc], td = A[jd];
					cnum x = ta + tc, y = tb + td;
					cnum z = ta - tc, t = (tb - td) * roots[1];
					A[ja] = x + y, A[jb] = x - y;
					A[jc] = z + t, A[jd] = z - t;
				}
			}
			// jh >= 1
			for (int jh = 1; jh < u; ++jh) {
				int ja = jh * v * 4;
				int jb = ja + v;
				int jc = jb + v;
				int jd = jc + v;
				int je = jb;
				cnum jw = roots[jh];
				cnum jx = roots[jh << 1];
				cnum wx = jw * jx;
				for (; ja < je; ++ja, ++jb, ++jc, ++jd) {
					cnum ta = A[ja], tb = A[jb] * jx, tc = A[jc] * jw, td = A[jd] * wx;
					cnum x = ta + tc, y = tb + td;
					cnum z = ta - tc, t = (tb - td) * roots[1];
					A[ja] = x + y, A[jb] = x - y;
					A[jc] = z + t, A[jd] = z - t;
				}
			}
			u <<= 2, v >>= 2;
		}
	}
	template <typename Vector> static void ifft(Vector& A, int k) {
		if (int(A.size()) <= 1) return;
		if (k == 1) {
			cnum v = A[1];
			A[1] = A[0] - A[1];
			A[0] = A[0] + v;
			return;
		}
		int u = 1 << (k - 2);
		int v = 1;
		while (u) {
			// jh = 0
			{
				int ja = 0;
				int jb = v;
				int jc = jb + v;
				int jd = jc + v;
				for (; ja < v; ++ja, ++jb, ++jc, ++jd) {
					cnum ta = A[ja], tb = A[jb], tc = A[jc], td = A[jd];
					cnum x = ta + tb, y = tc + td;
					cnum z = ta - tb, t = (tc - td) * roots[1].conj();
					A[ja] = x + y, A[jc] = x - y;
					A[jb] = z + t, A[jd] = z - t;
				}
			}
			// jh >= 1
			for (int jh = 1; jh < u; ++jh) {
				int ja = (jh * v) << 2;
				int jb = ja + v;
				int jc = jb + v;
				int jd = jc + v;
				int je = jb;
				cnum jw = roots[jh].conj();
				cnum jx = roots[jh << 1].conj();
				cnum wj = roots[(jh << 1) + 1].conj();
				for (; ja < je; ++ja, ++jb, ++jc, ++jd) {
					cnum ta = A[ja], tb = A[jb], tc = A[jc], td = A[jd];
					cnum x = ta + tb, y = tc + td;
					cnum z = (ta - tb) * jx, t = (tc - td) * wj;
					A[ja] = x + y, A[jc] = (x - y) * jw;
					A[jb] = z + t, A[jd] = (z - t) * jw;
				}
			}
			u >>= 2, v <<= 2;
		}
		if (k & 1) {
			u = 1 << (k - 1);
			for (int j = 0; j < u; j++) {
				cnum val = A[j] - A[j + u];
				A[j] = A[j] + A[j + u];
				A[j + u] = val;
			}
		}
	}

	template <typename U, typename V> static void fft_double(U& A, V& B, int k) {
		fft(A, k);
		B[0] = cnum{A[0].y * 2.0, 0};
		A[0] = cnum{A[0].x * 2.0, 0};
		B[1] = cnum{A[1].y * 2.0, 0};
		A[1] = cnum{A[1].x * 2.0, 0};
		for (int i = 2, y = 2; y < (1 << k); y <<= 1) {
			for (; i < 2 * y; i += 2) {
				int j = i ^ (y - 1);
				B[i] = (A[j].conj() - A[i]).prep_ccw();
				A[i] = (A[j].conj() + A[i]);
				B[j] = B[i].conj();
				A[j] = A[i].conj();
			}
		}
	}
};
template <typename num> std::vector<typename fft_complex<num>::cnum> fft_complex<num>::roots;

template <typename num> struct fft_complex_multiplier {
	template <typename U, typename V, typename P>
	static void multiply(const U& a, const int& sza, const V& b, const int& szb, P& res) {
		using cnum = complex<double>;
		std::vector<cnum>& fa = fft<cnum>::scratch_a;
		std::vector<cnum>& fb = fft<cnum>::scratch_b;

		int s = sza + szb - 1;
		int c = ceil_log_2(s);
		int n = 1 << c;
		fft_complex<num>::init(c);
		auto round = [](double x) -> num { return num(x + (x > 0 ? 0.5 : -0.5)); };

		if (int(fa.size()) < n) fa.resize(n);
		std::fill(fa.begin(), fa.begin() + n, cnum());
		for (int i = 0; i < sza; i++) fa[i].x = double(a[i]);
		for (int i = 0; i < szb; i++) fa[i].y = double(b[i]);
		fft_complex<num>::fft(fa, c);

		fa[0].y = 4.0 * fa[0].x * fa[0].y;
		fa[1].y = 4.0 * fa[1].x * fa[1].y;
		fa[0].x = fa[1].x = 0.0;
		for (int i = 2; i < n; i += 2) {
			int m = 1 << (31 - __builtin_clz(i));
			int j = i ^ (m - 1);
			fa[i] = (fa[i] + fa[j].conj()) * (fa[i] - fa[j].conj());
			fa[j] = -fa[i].conj();
		}

		if (int(fb.size()) < (n >> 1)) fb.resize(n >> 1);
		for (int j = 0; j < (n >> 1); j++) {
			int i = j << 1;
			cnum x = fa[i] + fa[i + 1];
			cnum y = (fa[i] - fa[i + 1]) * fft_complex<num>::roots[j].conj();
			fb[j] = x + y.prep_ccw();
		}
		fft_complex<num>::ifft(fb, c - 1);

		for (int i = 0; i < s; i++) {
			if (i & 1) {
				res[i] = round(-fb[i >> 1].x / (4.0 * n));
			} else {
				res[i] = round(fb[i >> 1].y / (4.0 * n));
			}
		}
	}
};

template <typename num> struct fft_complex_double_multiplier {
	template <typename U, typename V, typename P>
	static void multiply(const U& a, const int& sza, const V& b, const int& szb, P& res) {
		using cnum = complex<double>;
		std::vector<cnum>& fa = fft<cnum>::scratch_a;
		std::vector<cnum>& fb = fft<cnum>::scratch_b;

		int s = sza + szb - 1;
		int c = ceil_log_2(s);
		int n = 1 << c;
		fft_complex<num>::init(c);

		if (int(fa.size()) < n) fa.resize(n);
		std::fill(fa.begin(), fa.begin() + n, cnum());
		for (int i = 0; i < sza; i++) fa[i].x = double(a[i]);
		for (int i = 0; i < szb; i++) fa[i].y = double(b[i]);
		fft_complex<num>::fft(fa, c);

		fa[0].y = 4.0 * fa[0].x * fa[0].y;
		fa[1].y = 4.0 * fa[1].x * fa[1].y;
		fa[0].x = fa[1].x = 0.0;
		for (int i = 2; i < n; i += 2) {
			int m = 1 << (31 - __builtin_clz(i));
			int j = i ^ (m - 1);
			fa[i] = (fa[i] + fa[j].conj()) * (fa[i] - fa[j].conj());
			fa[j] = -fa[i].conj();
		}

		if (int(fb.size()) < (n >> 1)) fb.resize(n >> 1);
		for (int j = 0; j < (n >> 1); j++) {
			int i = j << 1;
			cnum x = fa[i] + fa[i + 1];
			cnum y = (fa[i] - fa[i + 1]) * fft_complex<num>::roots[j].conj();
			fb[j] = x + y.prep_ccw();
		}
		fft_complex<num>::ifft(fb, c - 1);

		for (int i = 0; i < s; i++) {
			if (i & 1) {
				res[i] = -fb[i >> 1].x / (4.0 * n);
			} else {
				res[i] = fb[i >> 1].y / (4.0 * n);
			}
		}
	}
};

template <typename num> struct fft_complex_mod_multiplier {
	template <typename U, typename V, typename P>
	static void multiply(const U& a, const int& sza, const V& b, const int& szb, P& res) {
		static constexpr uint64_t B = 32000;
		using cnum = complex<double>;

		std::vector<cnum>& fa = fft<cnum>::scratch_a;
		std::vector<cnum>& fb = fft<cnum>::scratch_b;
		static std::vector<cnum> fc;
		static std::vector<cnum> fd;

		int s = sza + szb - 1;
		int c = ceil_log_2(s);
		int n = 1 << c;
		fft_complex<num>::init(c);
		auto round = [](double x) { return uint64_t(x + 0.5); };

		if (int(fa.size()) < n) fa.resize(n);
		if (int(fb.size()) < n) fb.resize(n);
		if (int(fc.size()) < n) fc.resize(n);
		if (int(fd.size()) < n) fd.resize(n);
		using num_t = typename num::value_type;
		for (int i = 0; i < sza; i++) {
			auto v = num_t(a[i]);
			fa[i] = cnum{double(v % B), double(v / B)};
		}
		std::fill(fa.begin() + sza, fa.begin() + n, cnum());
		for (int i = 0; i < szb; i++) {
			auto v = num_t(b[i]);
			fc[i] = cnum{double(v % B), double(v / B)};
		}
		std::fill(fc.begin() + szb, fc.begin() + n, cnum());

		fft_complex<num>::fft_double(fa, fb, c);
		fft_complex<num>::fft_double(fc, fd, c);
		for (int i = 0; i < n; i++) {
			cnum v = fa[i] * fc[i] + (fb[i] * fd[i]).prep_ccw();
			fd[i] = fa[i] * fd[i] + (fb[i] * fc[i]).prep_ccw();
			fc[i] = v;
		}
		fft_complex<num>::ifft(fc, c);
		fft_complex<num>::ifft(fd, c);

		double v = 1.0 / (4.0 * n);
		for (int i = 0; i < s; i++) {
			fc[i].x *= v, fc[i].y *= v;
			fd[i].x *= v, fd[i].y *= v;
			auto sa = round(fc[i].x);
			auto sb = round(fd[i].x) + round(fd[i].y);
			auto sc = round(fc[i].y);

			res[i] += num(sa);
			res[i] += num(sb) * B;
			res[i] += num(sc) * B * B;
		}
	}
};

// TODO: Wrap everything inside numeric
template <typename num> struct numeric {
	typename num::value_type mod;
	num root;
	int base;
	std::vector<num> a;
	std::vector<num> b;
	std::vector<num> fa;
	std::vector<num> fb;

	constexpr void primitive_root() {
		mod = num::mod();
		base = __builtin_ctzll(mod - 1);
		a.resize(base); fa.resize(base);
		b.resize(base); fb.resize(base);
		// NOTE: Some special cases and not need to be calculated
		if (mod == 998244353) { root = 3; return; }
		if (mod == 18446744069414584321LLU) { root = 7; return; }
		// TODO: Calculate the primitive root of mod
		uint64_t idxs[32] = {};
		int idx = 0;
		uint64_t m = mod - 1;
		for (uint64_t x = 2; x * x <= m; ++x) {
			if (m % x == 0) {
				idxs[idx++] = x;
				while (m % x == 0) m /= x;
			}
		}
		if (m != 1) idxs[idx++] = m;
		root = 2;
		while (true) {
			bool v = true;
			for (int i = 0; i < idx; ++i) {
				if (bin_pow(root, (mod - 1) / idxs[i]) == num(1)) {
					v = false; break;
				}
			}
			if (v) return;
			++root;
		}
	}

	constexpr void init(const int& k) {
		if (mod == num::mod()) return;
		primitive_root();
		a[k - 1] = bin_pow(root, (mod - 1) >> k);
		b[k - 1] = static_cast<num>(1) / a[k - 1];
		for (int i = k - 2; i > 0; --i) a[i] = a[i + 1] * a[i + 1], b[i] = b[i + 1] * b[i + 1];
		fa[1] = a[1], fb[1] = b[1], fa[2] = a[2], fb[2] = b[2];
		for (int i = 3; i < k; ++i) {
			fa[i] = fa[i - 1] * b[i - 2] * a[i];
			fb[i] = fb[i - 1] * a[i - 2] * b[i];
		}
	}

	numeric() { init(base); }
};

template <typename num> struct fft_numeric {
	static numeric<num> root;

	static constexpr void init(const int& k) { root.init(k); }

	template <typename Vector> static constexpr void fft(Vector& A, int k) {
		if (k == 0) return;
		if (k == 1) {
			auto v = A[1];
			A[1] = A[0] - A[1];
			A[0] = A[0] + v;
			return;
		}
		if (k & 1) {
			int v = 1 << (k - 1);
			for (int j = 0; j < v; ++j) {
				num u = A[j + v];
				A[j + v] = A[j] - u;
				A[j] += u;
			}
		}
		int u = 1 << (2 + (k & 1));
		int v = 1 << (k - 2 - (k & 1));
		num one = num(1);
		auto s = root.fa[1];
		while (v) {
			{
				int ja = 0;
				int jb = v;
				int jc = jb + v;
				int jd = jc + v;
				for (; ja < v; ++ja, ++jb, ++jc, ++jd) {
					auto ta = A[ja], tb = A[jb], tc = A[jc], td = A[jd];
					auto x = ta + tc, y = tb + td;
					auto z = ta - tc, t = (tb - td) * s;
					A[ja] = x + y, A[jb] = x - y;
					A[jc] = z + t, A[jd] = z - t;
				}
			}
			auto ow = one, ox = one * root.fa[2], xw = one;
			for (int jh = 4; jh < u;) {
				ow = ox * ox, xw = ow * ox;
				int ja = jh * v;
				int jb = ja + v;
				int jc = jb + v;
				for (; ja < jb; ++ja, ++jc) {
					auto ta = A[ja], tb = A[ja + v] * ox, tc = A[jc] * ow, td = A[jc + v] * xw;
					auto x = ta + tc, y = tb + td;
					auto z = ta - tc, t = (tb - td) * s;
					A[ja] = x + y, A[ja + v] = x - y;
					A[jc] = z + t, A[jc + v] = z - t;
				}
				ox *= root.fa[__builtin_ctzll((jh += 4))];
			}
			u <<= 2;
			v >>= 2;
		}
	}
	template <typename Vector> static constexpr void ifft(Vector& A, int k) {
		if (k == 0) return;
		if (k == 1) {
			auto v = A[1];
			A[1] = A[0] - A[1];
			A[0] = A[0] + v;
			return;
		}
		int u = 1 << (k - 2);
		int v = 1;
		auto one = num(1);
		auto s = root.fb[1];
		while (u) {
			{
				int ja = 0;
				int jb = v;
				int jc = v + v;
				int jd = jc + v;
				for (; ja < v; ++ja, ++jb, ++jc, ++jd) {
					auto ta = A[ja], tb = A[jb], tc = A[jc], td = A[jd];
					auto x = ta + tb, y = tc + td;
					auto z = ta - tb, t = (tc - td) * s;
					A[ja] = x + y, A[jc] = x - y;
					A[jb] = z + t, A[jd] = z - t;
				}
			}
			auto ow = one, xw = one * root.fb[2], wy = one;
			u <<= 2;
			for (int jh = 4; jh < u;) {
				ow = xw * xw, wy = xw * s;
				int ja = jh * v;
				int jb = ja + v;
				int jc = jb + v;
				for (; ja < jb; ++ja, ++jc) {
					auto ta = A[ja], tb = A[ja + v], tc = A[jc], td = A[jc + v];
					auto x = ta + tb, y = tc + td;
					auto z = (ta - tb) * xw, t = (tc - td) * wy;
					A[ja] = x + y, A[jc] = (x - y) * ow;
					A[ja + v] = z + t, A[jc + v] = (z - t) * ow;
				}
				xw *= root.fb[__builtin_ctzll(jh += 4)];
			}
			u >>= 4;
			v <<= 2;
		}
		if (k & 1) {
			u = 1 << (k - 1);
			for (int j = 0; j < u; ++j) {
				auto ajv = A[j] - A[j + u];
				A[j] += A[j + u];
				A[j + u] = ajv;
			}
		}
	}
};
template <typename num> numeric<num> fft_numeric<num>::root;

template <typename num> struct fft_numeric_multiplier {
	template <typename U, typename V, typename P>
	static void multiply(const U& a, const int& sza, const V& b, const int& szb, P& r) {
		std::vector<num>& fa = fft<num>::scratch_a;
		std::vector<num>& fb = fft<num>::scratch_b;

		int s = sza + szb - 1;
		int c = ceil_log_2(s);
		int n = 1 << c;
		fft_numeric<num>::init(c);
		if (int(fa.size()) < n) fa.resize(n);
		std::copy(a.begin(), a.begin() + sza, fa.begin());
		std::fill(fa.begin() + sza, fa.begin() + n, num(0));
		fft_numeric<num>::fft(fa, c);
		if (sza == szb && std::equal(a.begin(), a.begin() + sza, b.begin(), b.begin() + szb)) {
			for (int i = 0; i < n; ++i) fa[i] *= fa[i];
		} else {
			if (int(fb.size()) < n) fb.resize(n);
			std::copy(b.begin(), b.begin() + szb, fb.begin());
			std::fill(fb.begin() + szb, fb.begin() + n, num(0));
			fft_numeric<num>::fft(fb, c);
			for (int i = 0; i < n; ++i) fa[i] *= fb[i];
		}
		fft_numeric<num>::ifft(fa, c);
		num inv_sz = num(1) / num(n);
		for (int i = 0; i < s; ++i) r[i] = fa[i] * inv_sz;
	}
};

template <typename num, typename multiplier> struct multiply_inverser {
	template <typename U, typename P>
	static void inverse(const U& a, const int& sza, P& b) {
		if (sza == 0) return;
		auto k = ceil_log_2(sza);
		int s = 1 << k;
		b[0] = num(1) / num(a[0]);
		std::vector<num> tmp(sza << 1);
		for (int n = 1, m = 1; n < sza; ++m) {
			multiplier::multiply(b, n, b, n, tmp);
			int nn = std::min(n << 1, sza);
			multiplier::multiply(tmp, nn, a, nn, tmp);
			for (int i = n; i < nn; ++i) b[i] = -tmp[i];
			n = nn;
		}
	}
};

template <typename num> struct fft_numeric_inverser {
	template <typename U, typename P>
	static void inverse(const U& a, const int& sza, P& b) {
		std::vector<num>& fa = fft<num>::scratch_a;
		std::vector<num>& fb = fft<num>::scratch_b;

		if (sza == 0) return;
		int c = ceil_log_2(sza) + 1;
		int s = 1 << c;
		fft_numeric<num>::init(c);
		if (int(fa.size()) < s) fa.resize(s);
		if (int(fb.size()) < s) fb.resize(s);
		fb[0] = num(1) / num(a[0]);
		for (int n = 1, m = 2; n < sza; ++m) {
			std::fill(fb.begin() + n, fb.begin() + (n << 2), num(0));
			n <<= 1;
			int nn = std::min(n, sza);
			std::copy(a.begin(), a.begin() + nn, fa.begin());
			std::fill(fa.begin() + nn, fa.begin() + (n << 1), num(0));
			fft_numeric<num>::fft(fb, m);
			fft_numeric<num>::fft(fa, m);
			num d = num(1) / num(n << 1);
			for (int i = 0; i < (n << 1); ++i) fb[i] = fb[i] * (2 - fa[i] * fb[i]) * d;
			fft_numeric<num>::ifft(fb, m);
		}
		std::copy(fb.begin(), fb.begin() + sza, b.begin());
	}
};

template <typename T, typename multiplier, typename U, typename V>
static inline std::vector<T> multiply(const U& a, const V& b) {
	int N = int(a.size()); int M = int(b.size());
	std::vector<T> res(N + M - 1);
	multiplier::multiply(a, N, b, M, res);
	return res;
}
template <typename T, typename U, typename V>
static inline std::vector<T> naive_multiply(const U& a, const V& b) {
	int N = int(a.size()); int M = int(b.size());
	std::vector<T> res(N + M - 1);
	naive_multiplier<T>::multiply(a, N, b, M, res);
	return res;
}
template <typename T, typename U, typename V>
static inline std::vector<T> fft_complex_multiply(const U& a, const V& b) {
	int N = int(a.size()); int M = int(b.size());
	std::vector<T> res(N + M - 1);
	fft_complex_multiplier<T>::multiply(a, N, b, M, res);
	return res;
}
template <typename T, typename U, typename V>
static inline std::vector<T> fft_complex_double_multiply(const U& a, const V& b) {
	int N = int(a.size()); int M = int(b.size());
	std::vector<T> res(N + M - 1);
	fft_complex_double_multiplier<T>::multiply(a, N, b, M, res);
	return res;
}
template <typename T, typename U, typename V>
static inline std::vector<T> fft_complex_mod_multiply(const U& a, const V& b) {
	int N = int(a.size()); int M = int(b.size());
	std::vector<T> res(N + M - 1);
	fft_complex_mod_multiplier<T>::multiply(a, N, b, M, res);
	return res;
}
template <typename T, typename U, typename V>
static inline std::vector<T> fft_numeric_multiply(const U& a, const V& b) {
	int N = int(a.size()); int M = int(b.size());
	std::vector<T> res(N + M - 1);
	fft_numeric_multiplier<T>::multiply(a, N, b, M, res);
	return res;
}

template <typename T, typename inverser, typename U> static inline std::vector<T> inverse(const U& a) {
	int N = int(a.size());
	std::vector<T> res(N, T(0));
	inverser::inverse(a, N, res);
	return res;
}
template <typename T, typename U> static inline std::vector<T> naive_inverse(const U& a) {
	int N = int(a.size());
	std::vector<T> res(N, T(0));
	multiply_inverser<T, naive_multiplier<T>>::inverse(a, N, res);
	return res;
}
template <typename T, typename U> static inline std::vector<T> fft_complex_inverse(const U& a) {
	int N = int(a.size());
	std::vector<T> res(N, T(0));
	multiply_inverser<T, fft_complex_multiplier<T>>::inverse(a, N, res);
	return res;
}
template <typename T, typename U> static inline std::vector<T> fft_complex_double_inverse(const U& a) {
	int N = int(a.size());
	std::vector<T> res(N, T(0));
	multiply_inverser<T, fft_complex_double_multiplier<T>>::inverse(a, N, res);
	return res;
}
template <typename T, typename U> static inline std::vector<T> fft_complex_mod_inverse(const U& a) {
	int N = int(a.size());
	std::vector<T> res(N, T(0));
	multiply_inverser<T, fft_complex_mod_multiplier<T>>::inverse(a, N, res);
	return res;
}
template <typename T, typename U> static inline std::vector<T> fft_numeric_inverse(const U& a) {
	int N = int(a.size());
	std::vector<T> res(N, T(0));
	fft_numeric_inverser<T>::inverse(a, N, res);
	return res;
}

} // namespace fft
