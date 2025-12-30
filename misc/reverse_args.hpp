namespace std {

template <typename F> struct reverse_args_t {
	F f;
	template <typename T, typename P> constexpr decltype(auto) operator () (T&& a, P&& b) & {
#if __cpp_lib_invoke >= 201411L
		return std::invoke(f, std::forward<P>(b), std::forward<T>(a));
#else
		return f(std::forward<P>(b), std::forward<T>(a));
#endif
	}
	template <typename T, typename P> constexpr decltype(auto) operator () (T&& a, P&& b) const& {
#if __cpp_lib_invoke >= 201411L
		return std::invoke(f, std::forward<P>(b), std::forward<T>(a));
#else
		return f(std::forward<P>(b), std::forward<T>(a));
#endif
	}
	template <typename T, typename P> constexpr decltype(auto) operator () (T&& a, P&& b) && {
#if __cpp_lib_invoke >= 201411L
		return std::invoke(std::move(f), std::forward<P>(b), std::forward<T>(a));
#else
		return std::move(f)(std::forward<P>(b), std::forward<T>(a));
#endif
	}
	template <typename T, typename P> constexpr decltype(auto) operator () (T&& a, P&& b) const&& {
#if __cpp_lib_invoke >= 201411L
		return std::invoke(std::move(f), std::forward<P>(b), std::forward<T>(a));
#else
		return std::move(f)(std::forward<P>(b), std::forward<T>(a));
#endif
	}
};

template <typename F> constexpr reverse_args_t<std::decay_t<F>> reverse_args(F&& f) {
	return { std::forward<F>(f) };
}

} // namespace std
