// NOTE: This supports checking prime or not for single numbers
// If for multiple checking, please use another ways (i.e Eratosthenes sieve, ...)
namespace miller_rabin {

template <typename T> struct miller_rabin_mod {
	using value_type = T;
	static T value;
};
template <typename T> T miller_rabin_mod<T>::value;

template <typename T> bool is_prime(const T& N) {
	if (N < 2) return false;
	for (const auto& x : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29}) {
		if (N == x) return true;
		if (N % x == 0) return false;
	}
	if (N < 961) return true;

	uint32_t s = __builtin_ctzll(N - 1);
	uint64_t d = (N - 1) >> s;

	auto get_bases = [&]() -> std::vector<uint64_t> {
		// Copied from: https://miller-rabin.appspot.com/
		if (N < 341531) return {9345883071009581737LLU};
		if (N < 1050535501) return {336781006125, 9639812373923155};
		if (N < 350269456337) return {4230279247111683200, 14694767155120705706LLU,
						16641139526367750375LLU};
		if (N < 55245642489451) return {2, 141889084524735, 1199124725622454117,
						11096072698276303650LLU};
		if (N < 7999252175582851) return {2, 4130806001517, 149795463772692060,
						186635894390467037, 3967304179347715805};
		if (N < 585226005592931977) return {2, 123635709730000, 9233062284813009,
						43835965440333360, 761179012939631437, 1263739024124850375};
		return {2, 325, 9375, 28178, 450775, 9780504, 1795265022};
	};

	miller_rabin_mod<T>::value = N;
	for (const auto& a : get_bases()) {
		modnum<miller_rabin_mod<T>, barrett_multiplier<T>> cur = a;
		if (cur() == 0) continue;
		cur = bin_pow(cur, d);
		if (cur() == 1) continue;
		bool seen = true;
		for (int i = 0; i < s; ++ i) {
			if (cur() == N - 1) {
				seen = false;
				break;
			}
			cur *= cur;
		}
		if (seen) return false;
	}
	return true;
}

} // namespace miller_rabin