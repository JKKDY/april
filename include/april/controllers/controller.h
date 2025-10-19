#pragma once

#include <concepts>

namespace april::controller  {

	class Controller {
	public:

		explicit Controller(const size_t call_frequency) : call_frequency_m(call_frequency) {}

		[[nodiscard]] size_t call_frequency() const { return call_frequency_m; }

		void apply(this auto&& self, size_t step) {

		}

	private:
		size_t call_frequency_m{};
	};



	// define controller concept
	template <class BC>
	concept IsController = std::derived_from<BC, Controller>;

	// define controller Pack
	template<IsController...>
	struct ControllerPack {};

	// constrained variable template
	template<class... Cs>
	requires (IsController<Cs> && ...)
	inline constexpr ControllerPack<Cs...> controllers {};
}
