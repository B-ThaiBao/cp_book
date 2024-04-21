template <typename T, const size_t row_t, const size_t col_t>
struct matrix : public std::array<T, row_t * col_t> {
	using std::array<T, row_t * col_t>::array;
	using std::array<T, row_t * col_t>::operator [];

	constexpr size_t row() const { return row_t; }
	constexpr size_t col() const { return col_t; }

	matrix(const T& t) { (*this).fill(t); }
	matrix(std::initializer_list<T> c) {
		// Copied from: https://stackoverflow.com/questions/4118025
		std::uninitialized_copy(c.begin(), c.end(), (*this).begin());
	}

	inline T& operator [] (const std::array<int, 2>& c) { return (*this)[c[0] * col() + c[1]]; }
	inline const T& operator [] (const std::array<int, 2>& c) const { return (*this)[c[0] * col() + c[1]]; }
	inline T& operator () (const int& x, const int& y) { return (*this)[x * col() + y]; }
	inline const T& operator () (const int& x, const int& y) const { return (*this)[x * col() + y]; }

	template <typename U>
	friend U& operator >> (U& in, matrix& mat) {
		for (size_t i = 0; i < mat.size(); ++ i) in >> mat[i];
		return in;
	}
	template <typename U>
	friend U& operator << (U& out, const matrix& mat) {
		for (int i = 0; i < int(mat.row()); ++ i) {
			for (int j = 0; j < int(mat.col()); ++ j) {
				out << mat(i, j) << ' ';
			}
			if (i + 1 != int(mat.row())) out << '\n';
		}
		return out;
	}

	matrix operator - () const {
		matrix res;
		for (size_t i = 0; i < (*this).size(); ++ i) res[i] = - (*this)[i];
		return res;
	}
	matrix operator + () const { return matrix(*this); }

	matrix& operator += (const matrix& o) {
		for (size_t i = 0; i < (*this).size(); ++ i) (*this)[i] += o[i];
		return *this;
	}
	matrix& operator -= (const matrix& o) {
		for (size_t i = 0; i < (*this).size(); ++ i) (*this)[i] -= o[i];
		return *this;
	}
	matrix& operator /= (const matrix& o) {
		for (size_t i = 0; i < (*this).size(); ++ i) (*this)[i] /= o[i];
		return *this;
	}
	matrix& operator *= (const matrix& o) { return *this = (*this) * o; }

	friend matrix operator + (const matrix& a, const matrix& b) { return matrix(a) += b; }
	friend matrix operator - (const matrix& a, const matrix& b) { return matrix(a) -= b; }
	friend matrix operator / (const matrix& a, const matrix& b) { return matrix(a) /= b; }

	template <typename U> friend matrix bin_pow(const matrix& a, const U& b) {
		// assert(b >= 0);
		matrix x = a, res(0);
		for (int i = 0; i < res.row(); ++ i) res(i, i) = 1;
		U p = b;
		while (p > 0) {
			if (p & 1) res *= x;
			x *= x;
			p >>= 1;
		}
		return res;
	}

	// NOTE: Below functions support for square matrix only
	matrix& operator += (const T& o) {
		for (size_t i = 0; i < (*this).row(); ++ i) (*this)(i, i) += o;
		return *this;
	}
	matrix& operator -= (const T& o) {
		for (size_t i = 0; i < (*this).row(); ++ i) (*this)(i, i) -= o;
		return *this;
	}
	matrix& operator /= (const T& o) {
		for (size_t i = 0; i < (*this).size(); ++ i) (*this)[i] /= o;
		return *this;
	}
	matrix& operator *= (const T& o) {
		for (size_t i = 0; i < (*this).size(); ++ i) (*this)[i] *= o;
		return *this;
	}

	friend matrix operator + (const matrix& a, const T& b) { return matrix(a) += b; }
	friend matrix operator - (const matrix& a, const T& b) { return matrix(a) -= b; }
	friend matrix operator / (const matrix& a, const T& b) { return matrix(a) /= b; }
	friend matrix operator * (const matrix& a, const T& b) { return matrix(a) *= b; }

	matrix& operator ++ () { return (*this) += 1; }
	matrix& operator -- () { return (*this) -= 1; }
	matrix operator ++ (int) { matrix res(*this); ++ (*this); return res; }
	matrix operator -- (int) { matrix res(*this); -- (*this); return res; }
};

template <typename T, typename P, size_t A, size_t B, size_t C, size_t D>
inline matrix<T, A, D> operator * (const matrix<T, A, B>& a, const matrix<P, C, D>& b) {
	static_assert(B == C, "The col size of first must equal row size of second");
	matrix<T, A, D> res;
	for (int i = 0; i < int(A); ++ i) {
		for (int j = 0; j < int(D); ++ j) {
			res(i, j) =  0;
			for (int k = 0; k < int(B); ++ k) {
				res(i, j) += a(i, k) * b(k, j);
			}
		}
	}
	return res;
}

template <typename T, size_t A, size_t B>
inline matrix<T, B, A> tranpose(const matrix<T, A, B>& a) {
	matrix<T, B, A> res;
	for (int i = 0; i < int(A); ++ i) {
		for (int j = 0; j < int(B); ++ j) {
			res(i, j) = a(j, i);
		}
	}
	return res;
}