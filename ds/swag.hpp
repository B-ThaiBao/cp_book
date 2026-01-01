/**
 * SWAG (Sliding Window Aggregation) !!
 *
 * Sources:
 * - https://codeforces.com/blog/entry/97396 (Square Root Optimization)
 * - https://codeforces.com/blog/entry/122003 (Minimum deque)
 *
 * Overview:
 * - Maintains a sliding window of elements with support for:
 *   + push_back / push_front
 *   + pop_front / pop_back
 *   + O(1) (amortized) aggregate query over the entire window
 * - Requirement: the combine operation must be *associative*:
 *     f(f(a, b), c) == f(a, f(b, c))
 *   (commutativity is NOT required, and no identity element is needed)
 *
 * Template parameters:
 * - T : value type stored in the window
 * - P : aggregate (payload) type
 * - F : associative binary functor with signature:
 *       P f(const P&, const P&)
 *
 * Internal structure & aggregation direction:
 * - swag_stack:
 *    Stores vector<pair<T, P>> where `second` is the aggregate of the stack
 *    from bottom -> top (i.e., insertion order within that stack):
 *        agg[i] = f(agg[i-1], P(value_i))
 *
 * - swag_queue (two stacks: pref + suff):
 *    The logical queue order is:  front ... back  ==  (pref.top -> pref.bottom) + (suff.bottom -> suff.top)
 *    * pref uses reverse_args(f) so its stored aggregates represent:
 *       pref.back().second = aggregate over pref elements in queue order (front -> ... -> boundary)
 *    * suff uses f so its stored aggregates represent:
 *       suff.back().second = aggregate over suff elements in queue order (boundary -> ... -> back)
 *    Whole-window aggregate:
 *     if one side empty -> that side's .back().second
 *     else -> f(pref.back().second, suff.back().second)
 *
 * - swag_deque (two stacks: pref + suff):
 *    The logical deque order is: front ... back == (pref.top -> pref.bottom) + (suff.bottom -> suff.top)
 *    * push_front goes to pref (so front element is pref.back().first)
 *    * push_back  goes to suff (so back  element is suff.back().first when suff non-empty)
 *    Aggregation is identical to queue:
 *     pref.back().second aggregates pref in left->right (front->boundary)
 *     suff.back().second aggregates suff in left->right (boundary->back)
 *     whole = f(pref.back().second, suff.back().second) when both non-empty
 *
 * Complexity:
 * - push / pop operations: O(1) amortized
 * - aggregate queries (front / back / whole window): O(1)
 *
 * Usage:
 *   auto f = [](int a, int b) { return a + b; };
 *   auto q = make_swag_queue<int, int>(f);
 *   q.push_back(1);
 *   q.push_back(2);
 *   q.push_back(3);
 *   auto [v, tot] = q.front(); // v = 1, tot = 6
 */
template <typename T, typename P, typename F> struct swag_stack : public std::vector<std::pair<T, P>> {
	F f;

	swag_stack(const F& fun = F()): f(fun) {}
	inline void push_back(const T& v) {
		if (this->empty()) {
			this->emplace_back(v, P(v));
		} else {
			this->emplace_back(v, f(this->back().second, P(v)));
		}
	}
};

template <typename T, typename P, typename F>
static constexpr swag_stack<T, P, std::decay_t<F>> make_swag_stack(F&& f) {
	return swag_stack<T, P, std::decay_t<F>>(std::forward<F>(f));
}

template <typename T, typename P, typename F> struct swag_queue {
	swag_stack<T, P, std::reverse_args_t<F>> pref;
	swag_stack<T, P, F> suff;
	F f;

	swag_queue(const F& fun = F()): pref(std::reverse_args(fun)), suff(fun), f(fun) {}
	template <typename Q> inline void reserve(const Q& N) { pref.reserve(N), suff.reserve(N); }
	constexpr size_t size() const { return pref.size() + suff.size(); }
	constexpr bool empty() const { return pref.empty() && suff.empty(); }
	inline void clear() { pref.clear(), suff.clear(); }

	inline void push_back(const T& v) { suff.push_back(v); }
	inline void pop_front() {
		assert(!this->empty());
		if (!pref.empty()) {
			pref.pop_back();
			return;
		}
		while (!suff.empty()) {
			if (int(suff.size()) > 1) pref.push_back(suff.back().first);
			suff.pop_back();
		}
	}
	inline std::pair<T, P> front() const {
		assert(!this->empty());
		if (pref.empty()) return {suff.front().first, suff.back().second};
		if (suff.empty()) return pref.back();
		return {pref.back().first, f(pref.back().second, suff.back().second)};
	}
	inline std::pair<T, P> back() const {
		assert(!this->empty());
		if (pref.empty()) return suff.back();
		if (suff.empty()) return {pref.front().first, pref.back().second};
		return {suff.back().first, f(pref.back().second, suff.back().second)};
	}
};

template <typename T, typename P, typename F>
static constexpr swag_queue<T, P, std::decay_t<F>> make_swag_queue(F&& f) {
	return swag_queue<T, P, std::decay_t<F>>(std::forward<F>(f));
}

template <typename T, typename P, typename F> struct swag_deque {
	swag_stack<T, P, std::reverse_args_t<F>> pref;
	swag_stack<T, P, F> suff;
	std::vector<T> buff;
	F f;

	swag_deque(const F& fun = F()): pref(std::reverse_args(fun)), suff(fun), f(fun) {}
	template <typename Q> inline void reserve(const Q& N) { pref.reserve(N), suff.reserve(N); }
	inline size_t size() const { return pref.size() + suff.size(); }
	inline bool empty() const { return pref.empty() && suff.empty(); }
	inline void clear() { pref.clear(), suff.clear(); }

	inline void push_back(const T& v) { suff.push_back(v); }
	inline void push_front(const T& v) { pref.push_back(v); }
	inline void pop_back() {
		assert(!this->empty());
		if (!suff.empty()) {
			suff.pop_back();
			return;
		}
		int N = int(pref.size()) >> 1;
		while (N--) buff.push_back(pref.back().first), pref.pop_back();
		while (!pref.empty()) {
			if (int(pref.size()) > 1) suff.push_back(pref.back().first);
			pref.pop_back();
		}
		while (!buff.empty()) pref.push_back(buff.back()), buff.pop_back();
	}
	inline void pop_front() {
		assert(!this->empty());
		if (!pref.empty()) {
			pref.pop_back();
			return;
		}
		int N = int(suff.size()) >> 1;
		while (N--) buff.push_back(suff.back().first), suff.pop_back();
		while (!suff.empty()) {
			if (int(suff.size()) > 1) pref.push_back(suff.back().first);
			suff.pop_back();
		}
		while (!buff.empty()) suff.push_back(buff.back()), buff.pop_back();
	}

	inline std::pair<T, P> front() const {
		assert(!this->empty());
		if (pref.empty()) return {suff.front().first, suff.back().second};
		if (suff.empty()) return pref.back();
		return {pref.back().first, f(pref.back().second, suff.back().second)};
	}
	inline std::pair<T, P> back() const {
		assert(!this->empty());
		if (pref.empty()) return suff.back();
		if (suff.empty()) return {pref.front().first, pref.back().second};
		return {suff.back().first, f(pref.back().second, suff.back().second)};
	}
};

template <typename T, typename P, typename F>
static constexpr swag_deque<T, P, std::decay_t<F>> make_swag_deque(F&& f) {
	return swag_deque<T, P, std::decay_t<F>>(std::forward<F>(f));
}
