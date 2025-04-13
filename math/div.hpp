template <typename T, typename U>
static inline typename std::common_type<T, U>::type floor_div(const T &a, const U &b) {
	return a / b - ((a ^ b) < 0 && a % b != 0);
}

template <typename T, typename U>
static inline typename std::common_type<T, U>::type ceil_div(const T &a, const U &b) {
	return a / b + ((a ^ b) > 0 && a % b != 0);
}
