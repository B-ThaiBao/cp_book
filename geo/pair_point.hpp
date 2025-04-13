template <typename T> struct pair_point {
	T x, y;
	pair_point() : x(0), y(0) {}
	pair_point(const T& x_) : x(x_), y(x_) {}
	pair_point(const T& x_, const T& y_) : x(x_), y(y_) {}
	template <typename U> explicit pair_point(const pair_point<U>& p) : x(p.x), y(p.y) {}
	pair_point(const std::pair<T, T>& p) : x(p.first), y(p.second) {}
	pair_point(const std::complex<T>& p) : x(p.real()), y(p.imag()) {}
	pair_point(const std::array<T, 2>& p) : x(p[0]), y(p[1]) {}
	pair_point(const std::tuple<T, T>& p) : x(std::get<0>(p)), y(std::get<1>(p)) {}

	template <typename U, typename V>
	friend pair_point make_pair_point(const U& x, const V& y) { return pair_point(x, y); }

	// Change corresponding type
	template <typename U, typename Q> explicit operator std::pair<U, Q> () const { return {x, y}; }
	template <typename U> explicit operator std::complex<U> () const { return {x, y}; }
	template <typename U> explicit operator std::array<U, 2> () const { return {x, y}; }
	template <typename U, typename Q> explicit operator std::tuple<U, Q> () const { return {x, y}; }

	// Stream and print debug
	friend std::ostream& operator << (std::ostream &os, const pair_point &p) {
		return os << '(' << p.x << ',' << p.y << ')';
	}
	friend std::istream& operator >> (std::istream &is, pair_point &p) { return is >> p.x >> p.y; }
	friend void __print(const pair_point &p) { std::cerr << p; }

	friend bool operator == (const pair_point &a, const pair_point &b) { return a.x == b.x && a.y == b.y; }
	template <typename U>
	friend bool operator == (const pair_point &a, const U& b) { return a == pair_point(b); }
	template <typename U>
	friend bool operator == (const U& a, const pair_point &b) { return pair_point(a) == b; }
	friend bool operator != (const pair_point &a, const pair_point &b) { return !(a == b); }
	template <typename U>
	friend bool operator != (const pair_point &a, const U& b) { return !(a == pair_point(b)); }
	template <typename U>
	friend bool operator != (const U& a, const pair_point &b) { return !(pair_point(a) == b); }

	friend bool operator < (const pair_point& a, const pair_point& b) {
		return std::tie(a.x, a.y) < std::tie(b.x, b.y);
	}
	template <typename U> friend bool operator < (const pair_point& a, const U& b) { return a < pair_point(b); }
	template <typename U> friend bool operator < (const U& a, const pair_point& b) { return pair_point(a) < b; }
	friend bool operator <= (const pair_point& a, const pair_point& b) {
		return std::tie(a.x, a.y) <= std::tie(b.x, b.y);
	}
	template <typename U> friend bool operator <= (const pair_point& a, const U& b) { return a <= pair_point(b); }
	template <typename U> friend bool operator <= (const U& a, const pair_point& b) { return pair_point(a) <= b; }

	friend bool operator > (const pair_point& a, const pair_point& b) {
		return std::tie(a.x, a.y) > std::tie(b.x, b.y);
	}
	template <typename U> friend bool operator > (const pair_point& a, const U& b) { return a > pair_point(b); }
	template <typename U> friend bool operator > (const U& a, const pair_point& b) { return pair_point(a) > b; }
	friend bool operator >= (const pair_point& a, const pair_point& b) {
		return std::tie(a.x, a.y) >= std::tie(b.x, b.y);
	}
	template <typename U> friend bool operator >= (const pair_point& a, const U& b) { return a >= pair_point(b); }
	template <typename U> friend bool operator >= (const U& a, const pair_point& b) { return pair_point(a) >= b; }

	pair_point operator + () const { return pair_point(+x, +y); }
	pair_point operator - () const { return pair_point(-x, -y); }

	// Arithmetic
	template <typename U> pair_point& operator += (const pair_point<U> &p) { x += p.x, y += p.y; return *this; }
	template <typename U> pair_point& operator -= (const pair_point<U> &p) { x -= p.x, y -= p.y; return *this; }
	template <typename U> pair_point& operator *= (const U& t) { x *= t, y *= t; return *this; }
	template <typename U> pair_point& operator /= (const U& t) { x /= t, y /= t; return *this; }

	template <typename U>
	pair_point<typename std::common_type<T, U>::type> operator + (const pair_point<U> &p) const {
		return {x + p.x, y + p.y};
	}
	template <typename U>
	pair_point<typename std::common_type<T, U>::type> operator - (const pair_point<U> &p) const {
		return {x - p.x, y - p.y};
	}
	template <typename U>
	pair_point<typename std::common_type<T, U>::type> operator * (const U& t) const {
		return {x * t, y * t};
	}
	template <typename U>
	friend pair_point<typename std::common_type<T, U>::type> operator * (const U& t, const pair_point &p) {
		return {t * p.x, t * p.y};
	}
	template <typename U>
	pair_point<typename std::common_type<T, U>::type> operator / (const U& t) const {
		return {x / t, y / t};
	}

	// Choose bigger type for result of dist
	T norm() const { return x * x + y * y; }
	template <typename U> typename std::common_type<T, U>::type norm(const pair_point<U> &p) const {
		return (x - p.x) * (x - p.x) + (y - p.y) * (y - p.y);
	}
	double abs() const { return std::sqrt(this->norm()); }
	long double absl() const { return sqrtl(this->dist2()); }
	template <typename U> double abs(const pair_point<U> &p) const { return std::sqrt(this->norm(p)); }
	template <typename U> long double absl(const pair_point<U> &p) const { return sqrtl(this->norm(p)); }
	pair_point unit() const { return *this / T(this->abs()); }
	pair_point unitl() const { return *this / T(this->absl()); }
	double arg() const { return std::atan2(y, x); }
	long double argl() const { return atan2l(y, x); }
	template <typename U>
	double arg(const pair_point<U>& p) const { return atan2(this->cross(p), this->dot(p)); }
	template <typename U>
	long double argl(const pair_point<U>& p) const { return atan2l(this->cross(p), this->dot(p)); }

	T int_norm() const { return std::__gcd(x, y); }
	pair_point int_unit() const { if (!x && !y) return *this; return *this / this->int_norm(); }

	// Convenient free-functions, mostly for generic interop
	friend T norm(const pair_point &p) { return p.norm(); }
	template <typename U>
	friend typename std::common_type<T, U>::type norm(const pair_point &p, const pair_point<U> &q) { return p.norm(q); }
	friend double abs(const pair_point &p) { return p.abs(); }
	friend long double absl(const pair_point &p) { return p.absl(); }
	template <typename U> friend double abs(const pair_point &p, const pair_point<U> &q) { return p.abs(q); }
	template <typename U> friend long double absl(const pair_point &p, const pair_point<U> &q) { return p.absl(q); }
	friend pair_point unit(const pair_point &p) { return p.unit(); }
	friend pair_point unitl(const pair_point &p) { return p.unitl(); }
	friend double arg(const pair_point &p) { return p.arg(); }
	friend long double argl(const pair_point &p) { return p.argl(); }
	template <typename U>
	friend double arg(const pair_point &p, const pair_point<U> &q) { return p.arg(q); }
	template <typename U>
	friend long double argl(const pair_point &p, const pair_point<U> &q) { return p.argl(q); }
	friend T int_norm(const pair_point &p) { return p.int_norm(); }
	friend pair_point int_unit(const pair_point &p) { return p.int_unit(); }

	// Rotate pi / 2 and 3 * pi / 2
	pair_point perp_cw() const { return pair_point(y, -x); }
	pair_point perp_ccw() const { return pair_point(-y, x); }
	friend pair_point perp_cw(const pair_point &p) { return p.perp_cw(); }
	friend pair_point perp_ccw(const pair_point &p) { return p.perp_ccw(); }

	// Dot and cross
	template <typename U>
	typename std::common_type<T, U>::type dot(const pair_point<U> &p) const { return x * p.x + y * p.y; }
	template <typename U>
	typename std::common_type<T, U>::type cross(const pair_point<U> &p) const { return x * p.y - y * p.x; }
	template <typename U, typename V>
	typename std::common_type<T, U, V>::type cross3(const pair_point<U> &b, const pair_point<V> &c) const {
		return (b - *this).cross(c - *this);
	}
	template <typename U>
	friend typename std::common_type<T, U>::type dot(const pair_point<T> &a, const pair_point<U> &b) {
		return a.dot(b);
	}
	template <typename U>
	friend typename std::common_type<T, U>::type cross(const pair_point<T> &a, const pair_point<U> &b) {
		return a.cross(b);
	}
	template <typename U, typename V>
	friend typename std::common_type<T, U, V>::type cross3(const pair_point<T> &a, const pair_point<U> &b, const pair_point<V> &c) {
		return a.cross3(b, c);
	}

	// Complex numbers and rotation
	pair_point conj() const { return pair_point(x, -y); }
	friend pair_point conj(const pair_point &p) { return p.conj(); }

	// Returns conj(a) * b
	template <typename U>
	pair_point<typename std::common_type<T, U>::type> dot_cross(const pair_point<U> &b) const {
		return {this->dot(b), this->cross(b)};
	}
	template <typename U>
	pair_point<typename std::common_type<T, U>::type> cmul(const pair_point<U> &b) const {
		return this->conj().dot_cross(b);
	}
	template <typename U>
	pair_point<typename std::common_type<T, U>::type> cdiv(const pair_point<U> &b) const {
		return b.dot_cross(*this) / b.norm();
	}
	template <typename U>
	friend pair_point<typename std::common_type<T, U>::type> dot_cross(const pair_point<T> &a, const pair_point<U> &b) {
		return a.dot_cross(b);
	}
	template <typename U>
	friend pair_point<typename std::common_type<T, U>::type> cmul(const pair_point<T> &a, const pair_point<U> &b) {
		return a.cmul(b);
	}
	template <typename U>
	friend pair_point<typename std::common_type<T, U>::type> cdiv(const pair_point<T> &a, const pair_point<U> &b) {
		return a.cdiv(b);
	}

	// Must be a unit vector; otherwise multiplies the result by abs(u)
	// When rotate arg --> rotate by point{cos(arg), sin(arg)}
	template <typename U>
	pair_point<typename std::common_type<T, U>::type> rotate(const pair_point<U> &p) const {
		return p.conj().dot_cross(*this);
	}
	template <typename U>
	pair_point<typename std::common_type<T, U>::type> unrotate(const pair_point<U> &p) const {
		return p.dot_cross(*this);
	}

	template <typename U>
	bool same_dir(const pair_point<U> &p) const { return this->cross(p) == 0 && this->dot(p) > 0; }
	template <typename U>
	friend bool same_dir(const pair_point<T> &a, const pair_point<U> &b) { return a.same_dir(b); }

	// check if 180 <= s..t < 360
	template <typename U>
	bool is_reflex(const pair_point<U>& p) const {
		auto c = this->cross(p); return c < 0 || (c == 0 && this->dot(p) < 0);
	}
	template <typename U>
	friend bool is_reflex(const pair_point<T> &a, const pair_point<U> &b) { return a.is_reflex(b); }

	// operator < (s,t) for angles in [base, base + 2pi)
	template <typename U, typename V>
	bool less_angle(const pair_point<U>& s, const pair_point<V>& t) const {
		int r = this->is_reflex(s) - this->is_reflex(t);
		return r < 0 || (r == 0 && s.cross(t) > 0);
	}
	template <typename U, typename V>
	friend bool less_angle(const pair_point<T>& base, const pair_point<U>& s, const pair_point<V>& t) {
		return base.less_angle(s, t);
	}

	auto cmp_angle() const {
		return [this](const auto &s, const auto &t) { return this->less_angle(s, t); };
	}
	friend auto cmp_angle(const pair_point &base) { return base.cmp_angle(); }

	// Cmp with center and dir
	template <typename U> auto cmp_angle_center(const pair_point<U>& dir) const {
		return [this, dir](const auto &s, const auto &t) {
			return angle_less(dir, s - *this, t - *this);
		};
	}
	template <typename U> friend auto cmp_angle_center(const pair_point<T>& base, const pair_point<U>& dir) {
		return base.cmp_angle_center(dir);
	}
};
