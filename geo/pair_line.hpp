template <typename T> struct pair_line : public std::array<pair_point<T>, 2> {
	using std::array<pair_point<T>, 2>::array;
	template <typename P, typename Q, typename R, typename S>
	pair_line(const P& a, const Q& b, const R& c, const S& d)
		: std::array<pair_point<T>, 2>{pair_point<T>(a, b), pair_point<T>(c, d)} {}
	template <typename U> pair_line(const pair_point<U> &a)
		: std::array<pair_point<T>, 2>{pair_point<T>(a), pair_point<T>(a)} {}
	template <typename U, typename V> pair_line(const pair_point<U> &a, const pair_point<V> &b)
		: std::array<pair_point<T>, 2>{pair_point<T>(a), pair_point<T>(b)} {}

	friend std::istream& operator >> (std::istream& in, pair_line<T>& l) {
		in >> l[0] >> l[1]; return in;
	}

	template <typename Q, typename R, typename S>
	friend pair_line<typename std::common_type<T, Q, R, S>::type> make_pair_line(const T& a, const Q& b, const R& c, const S& d) {
		return {a, b, c, d};
	}
	template <typename U>
	friend pair_line<typename std::common_type<T, U>::type> make_pair_line(const pair_point<T> &a, const pair_point<U> &b) {
		return {a, b};
	}

	pair_point<T> dir() const { return this->at(1) - this->at(0); }
	friend pair_point<T> dir(const pair_line &l) { return l.dir(); }

	template <typename U> bool on_line(const pair_point<U>& p) {
		return cross(p - this->at(0), this->dir()) == 0;
	}
	template <typename U> bool on_segment(const pair_point<U>& p) {
		auto PA = p - this->at(0), PB = p - this->at(1);
		return cross(PA, PB) == 0 && dot(PA, PB) <= 0;
	}

	// Return square of distance from p to this line
	template <typename K, typename U> K norm(const pair_point<U> &p) const {
		auto d = this->dir();
		K c = K(cross(p - this->at(0), d));
		return c * c / K(d.norm());
	}
	template <typename U> double abs(const pair_point<U> &p) const {
		auto d = this->dir();
		return std::abs(cross(p - this->at(0), d)) / d.abs();
	}
	template <typename U> long double absl(const pair_point<U> &p) const {
		auto d = this->dir();
		return std::abs(cross(p - this->at(0), d)) / d.absl();
	}

	// NOTE: This dist functions purely supports for line (doesn't allow segment and ray)
	// For segment and ray, please check whether projection is satisfied by dot product
	template <typename K, typename U>
	friend K norm(const pair_point<U>& p, const pair_line &l) { return l.norm<K>(p); }
	template <typename U>
	friend double abs(const pair_point<U>& p, const pair_line &l) { return l.abs(p); }
	template <typename U>
	friend long double absl(const pair_point<U>& p, const pair_line &l) { return l.absl(p); }

	template <typename K, typename U> pair_point<K> proj(const pair_point<U> &p) const {
		pair_point<K> d = pair_point<K>(this->dir());
		K c = K(cross(p - this->at(0), d));
		return pair_point<K>(this->at(0)) + d * (c / K(d.norm()));
	}
	template <typename K, typename U> pair_point<K> reflect(const pair_point<U> &p) const {
		return this->proj<K>(p) * 2 - p;
	}

	template <typename K, typename U>
	friend pair_point<K> proj(const pair_point<U> &p, const pair_line &l) { return l.proj<K>(p); }
	template <typename K, typename U>
	friend pair_point<K> reflect(const pair_point<U> &p, const pair_line &l) { return l.reflect<K>(p); }

	std::array<T, 3> coef() const {
		auto d = this->dir();
		return {-d.y, d.x, perp_cw(d).dot(this->at(0))};
	}

	// NOTE: Return the fractions of intersection points (i = c[0] / c[2] and j = c[1] / c[2])
	// If c[2] == 0, two case: c[0] == 0 --> concide, c[0] != 0 --> parallel (no intersect)
	// Otherwise, lines strictly intersect: p[0] + i * dir(p) or l[0] + j * dir(l)
	// Check i, j such that: 0 <= i <= 1 (segment), 0 <= i (ray), nothing (line)
	// For easier check, just check between c[0..2] because the result ensures that c[2] > 0
	template <typename U>
	friend std::array<typename std::common_type<T, U>::type, 3> intersect(const pair_line<T>& p, const pair_line<U> &l) {
		std::array<typename std::common_type<T, U>::type, 3> res;
		auto p_dir = p.dir();
		auto l_dir = l.dir();
		auto d = l[0] - p[0];
		res[0] = cross(d, l_dir);
		res[1] = cross(d, p_dir);
		res[2] = cross(p_dir, l_dir);
		if (res[2] < 0) res[0] = -res[0], res[1] = -res[1], res[2] = -res[2];
		return res;
	}
};
