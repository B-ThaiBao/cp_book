namespace std {

template <typename T, const int NDIMS> struct tensor_view {
	static_assert(NDIMS >= 0, "NDIMS must be nonnegative");

	std::array<int, NDIMS> shape;
	std::array<int, NDIMS> strides;
	T* data;

	tensor_view(const std::array<int, NDIMS>& shape_, const std::array<int, NDIMS>& strides_, T* data_)
		: shape(shape_), strides(strides_), data(data_) {}
	tensor_view() : shape{0}, strides{0}, data(nullptr) {}

	inline int flatten_index(const std::array<int, NDIMS>& idx) const {
		int res = 0;
		for (int i = 0; i < NDIMS; i++) {
#ifdef _GLIBCXX_DEBUG
			assert(0 <= idx[i] && idx[i] < shape[i]);
#endif
			res += idx[i] * strides[i];
		}
		return res;
	}

	T& operator[] (const std::array<int, NDIMS>& idx) const { return data[flatten_index(idx)]; }
	T& at(const std::array<int, NDIMS>& idx) const { return data[flatten_index_checked(idx)]; }

	template <const int D = NDIMS>
	typename std::enable_if<(0 < D), tensor_view<T, NDIMS - 1>>::type operator[] (int idx) const {
		std::array<int, NDIMS - 1> nshape; std::copy(shape.begin() + 1, shape.end(), nshape.begin());
		std::array<int, NDIMS - 1> nstrides; std::copy(strides.begin() + 1, strides.end(), nstrides.begin());
		T* ndata = data + (strides[0] * idx);
		return tensor_view<T, NDIMS - 1>(nshape, nstrides, ndata);
	}
	template <const int D = NDIMS>
	typename std::enable_if<(0 < D), tensor_view<T, NDIMS - 1>>::type at(int idx) const {
#ifdef _GLIBCXX_DEBUG
		assert(0 <= idx && idx < shape[0]);
#endif
		return operator[](idx);
	}

	template <const int D = NDIMS>
	typename std::enable_if<(0 == D), T&>::type operator * () const {
		return *data;
	}
};

template <typename T, const int NDIMS> struct tensor {
	static_assert(NDIMS >= 0, "NDIMS must be nonnegative");

	std::array<int, NDIMS> shape;
	std::array<int, NDIMS> strides;
	int len;
	T* data;

	tensor() : shape{0}, strides{0}, len(0), data(nullptr) {}
	explicit tensor(const std::array<int, NDIMS>& shape_, const T& t = T()) {
		shape = shape_;
		len = 1;
		for (int i = NDIMS - 1; i >= 0; i--) {
			strides[i] = len;
			len *= shape[i];
		}
		data = new T[len];
		std::fill(data, data + len, t);
	}
	tensor(const tensor& o) : shape(o.shape), strides(o.strides), len(o.len), data(new T[len]) {
		for (int i = 0; i < len; i++) {
			data[i] = o.data[i];
		}
	}

	inline void assign(const std::array<int, NDIMS>& shape_, const T& t = T()) {
		shape = shape_;
		len = 1;
		for (int i = NDIMS - 1; i >= 0; i--) {
			strides[i] = len;
			len *= shape[i];
		}
		delete[] data;
		data = new T[len];
		std::fill(data, data + len, t);
	}

	tensor& operator=(tensor&& o) noexcept {
		using std::swap;
		swap(shape, o.shape);
		swap(strides, o.strides);
		swap(len, o.len);
		swap(data, o.data);
		return *this;
	}
	tensor(tensor&& o) : tensor() {
		*this = std::move(o);
	}
	tensor& operator=(const tensor& o) {
		return *this = tensor(o);
	}
	~tensor() { delete[] data; }

	// View the data without copying it (similar to std::string_view)
	using view_t = tensor_view<T, NDIMS>;
	view_t view() { return tensor_view<T, NDIMS>(shape, strides, data); }
	operator view_t() { return view(); }
	friend inline view_t view(const tensor<T, NDIMS>& t) { return tensor_view<T, NDIMS>(t.shape, t.strides, t.data); }

	using const_view_t = tensor_view<const T, NDIMS>;
	const_view_t view() const { return tensor_view<const T, NDIMS>(shape, strides, data); }
	operator const_view_t() const { return view(); }
	friend inline const_view_t const_view(const tensor<T, NDIMS>& t) { return tensor_view<const T, NDIMS>(t.shape, t.strides, t.data); }

	inline int flatten_index(const std::array<int, NDIMS>& idx) const {
		int res = 0;
		for (int i = 0; i < NDIMS; i++) {
#ifdef _GLIBCXX_DEBUG
			assert(0 <= idx[i] && idx[i] < shape[i]);
#endif
			res += idx[i] * strides[i];
		}
		return res;
	}
	T& operator[] (const std::array<int, NDIMS>& idx) { return data[flatten_index(idx)]; }
	T& at(const std::array<int, NDIMS>& idx) { return data[flatten_index(idx)]; }
	const T& operator[] (const std::array<int, NDIMS>& idx) const { return data[flatten_index(idx)]; }
	const T& at(const std::array<int, NDIMS>& idx) const { return data[flatten_index(idx)]; }

	template <const int D = NDIMS>
	typename std::enable_if<(0 < D), tensor_view<T, NDIMS - 1>>::type operator[] (int idx) {
		std::array<int, NDIMS - 1> nshape; std::copy(shape.begin() + 1, shape.end(), nshape.begin());
		std::array<int, NDIMS - 1> nstrides; std::copy(strides.begin() + 1, strides.end(), nstrides.begin());
		T* ndata = data + (strides[0] * idx);
		return tensor_view<T, NDIMS - 1>(nshape, nstrides, ndata);
	}
	template <const int D = NDIMS>
	typename std::enable_if<(0 < D), tensor_view<T, NDIMS - 1>>::type at(int idx) {
		std::array<int, NDIMS - 1> nshape; std::copy(shape.begin() + 1, shape.end(), nshape.begin());
		std::array<int, NDIMS - 1> nstrides; std::copy(strides.begin() + 1, strides.end(), nstrides.begin());
		T* ndata = data + (strides[0] * idx);
		return tensor_view<T, NDIMS - 1>(nshape, nstrides, ndata);
	}

	template <const int D = NDIMS>
	typename std::enable_if<(0 < D), tensor_view<const T, NDIMS - 1>>::type operator[] (int idx) const {
		std::array<int, NDIMS - 1> nshape; std::copy(shape.begin() + 1, shape.end(), nshape.begin());
		std::array<int, NDIMS - 1> nstrides; std::copy(strides.begin() + 1, strides.end(), nstrides.begin());
		const T* ndata = data + (strides[0] * idx);
		return tensor_view<const T, NDIMS - 1>(nshape, nstrides, ndata);
	}
	template <const int D = NDIMS>
	typename std::enable_if<(0 < D), tensor_view<const T, NDIMS - 1>>::type at(int idx) const {
		std::array<int, NDIMS - 1> nshape; std::copy(shape.begin() + 1, shape.end(), nshape.begin());
		std::array<int, NDIMS - 1> nstrides; std::copy(strides.begin() + 1, strides.end(), nstrides.begin());
		const T* ndata = data + (strides[0] * idx);
		return tensor_view<const T, NDIMS - 1>(nshape, nstrides, ndata);
	}

	template <const int D = NDIMS>
	typename std::enable_if<(0 == D), T&>::type operator * () {
		return *view();
	}
	template <const int D = NDIMS>
	typename std::enable_if<(0 == D), const T&>::type operator * () const {
		return *view();
	}
};

} // namespace std
