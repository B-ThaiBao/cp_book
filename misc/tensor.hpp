namespace std {

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
		for (int i = NDIMS - 1 ; i >= 0; -- i) {
			strides[i] = len;
			len *= shape[i];
		}
		data = new T[len];
		std::fill(data, data + len, t);
	}
	explicit tensor(std::array<int, NDIMS>&& shape_, const T& t = T()) {
		std::swap(shape, shape_);
		len = 1;
		for (int i = NDIMS - 1 ; i >= 0; -- i) {
			strides[i] = len;
			len *= shape[i];
		}
		data = new T[len];
		std::fill(data, data + len, t);
	}
	inline void assign(const std::array<int, NDIMS>& shape_, const T& t = T()) {
		shape = shape_;
		len = 1;
		for (int i = NDIMS - 1 ; i >= 0; -- i) {
			strides[i] = len;
			len *= shape[i];
		}
		data = new T[len];
		std::fill(data, data + len, t);
	}
	inline void assign(std::array<int, NDIMS>&& shape_, const T& t = T()) {
		std::swap(shape, shape_);
		len = 1;
		for (int i = NDIMS - 1 ; i >= 0; -- i) {
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

	tensor& operator = (tensor&& o) noexcept {
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
	tensor& operator = (const tensor& o) {
		return *this = tensor(o);
	}
	~tensor() { delete[] data; }

	inline const std::array<int, NDIMS>& size() const { return shape; }
	inline const int& size(const int& dim) const { return shape[dim]; }
	template <typename U> friend U& operator >> (U& in, tensor& t) {
		for (int i = 0; i < t.len; ++ i) in >> t.data[i];
		return in;
	}

	inline bool operator == (const tensor& o) const {
		return shape == o.shape && std::equal(data, data + len, o.data);
	}
	inline bool operator != (const tensor& o) const { return !(*this == o); }

	inline int flatten_index(const std::array<int, NDIMS>& idx) {
		int res = 0;
		for (int i = 0; i < NDIMS; i++) { res += idx[i] * strides[i]; }
		return res;
	}
	inline T& operator [] (const std::array<int, NDIMS>& idx) { return data[flatten_index(idx)]; }
	inline T& at(const std::array<int, NDIMS>& idx) { return data[flatten_index(idx)]; }
	inline const T& operator[] (const std::array<int, NDIMS>& idx) const { return data[flatten_index(idx)]; }
	inline const T& at(const std::array<int, NDIMS>& idx) const { return data[flatten_index(idx)]; }

	inline int flatten_index(const std::array<int, NDIMS - 1>& idx) {
		int res = 0;
		for (int i = 0; i < NDIMS - 1; i++) { res += idx[i] * strides[i]; }
		return res;
	}
	inline T* operator [] (const std::array<int, NDIMS - 1>& idx) { return data + flatten_index(idx); }
	inline T* at(const std::array<int, NDIMS - 1>& idx) { return data + flatten_index(idx); }
};

} // namespace std