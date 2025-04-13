/**
 * LINE MULTISET !!!
 *
 * This container contains some lines for 'convex hull trick' in DP and work based on
 * line formula: ax + b (for max query: unchange, for min query: -a, -b instead):
 *   * Query: c = *lower_bound(x); max = c[0] * x + c[1], min = -c[0] * x - c[1]
**/
template <typename T> struct line : public std::array<T, 2> {
	mutable T x;
	friend inline bool operator < (const line &v, const T &o) { return v.x < o; }
	friend inline bool operator > (const line &v, const T &o) { return v.x > o; }
	friend inline bool operator <= (const line &v, const T &o) { return v.x <= o; }
	friend inline bool operator >= (const line &v, const T &o) { return v.x >= o; }
	friend inline bool operator == (const line &v, const T &o) { return v.x == o; }
	friend inline bool operator != (const line &v, const T &o) { return v.x != o; }

	friend inline bool operator < (const T &o, const line &v) { return o < v.x; }
	friend inline bool operator > (const T &o, const line &v) { return o > v.x; }
	friend inline bool operator <= (const T &o, const line &v) { return o <= v.x; }
	friend inline bool operator >= (const T &o, const line &v) { return o >= v.x; }
	friend inline bool operator == (const T &o, const line &v) { return o == v.x; }
	friend inline bool operator != (const T &o, const line &v) { return o != v.x; }
};

template <typename Line, typename Comp = std::less<>>
struct line_multiset : public std::multiset<Line, Comp> {
	using pt_t = typename Line::value_type;
	static constexpr pt_t INF = std::numeric_limits<pt_t>::max();

	template <typename T = pt_t>
	typename std::enable_if<std::is_integral<T>::value, T>::type divide(const T &a, const T &b) {
		return a / b - ((a ^ b) < 0 && a % b != 0);
	}
	template <typename T = pt_t>
	typename std::enable_if<std::is_floating_point<T>::value, T>::type divide(const T &a, const T &b) {
		return a / b;
	}

	template <typename Iterator> inline bool intersect(const Iterator &x, const Iterator &y) {
		if (y == this->end()) { x->x = INF; return false; }
		if (x->at(0) == y->at(0)) x->x = x->at(1) > y->at(1) ? INF : -INF;
		else x->x = divide(y->at(1) - x->at(1), x->at(0) - y->at(0));
		return x->x >= y->x;
	}
	template <typename T, typename U> inline void add_line(const T& a, const U &b) {
		auto z = this->insert({a, b, 0}), y = z++, x = y;
		while (intersect(y, z)) z = this->erase(z);
		if (x != this->begin() && intersect(--x, y)) intersect(x, y = this->erase(y));
		while ((y = x) != this->begin() && (--x)->x >= y->x) intersect(x, this->erase(y));
	}
};
