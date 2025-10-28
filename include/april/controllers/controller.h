#pragma once

#include <concepts>

#include "april/shared/trigger.h"
#include "april/core/context.h"

namespace april::controller  {

	class Controller {
	public:
		explicit Controller(shared::Trigger trig) : trigger(std::move(trig)) {}

		[[nodiscard]] bool should_trigger(const core::SimulationContext & sys) {
			return trigger(sys);
		}

		void dispatch_init(this auto && self, core::SimulationContext & sys) {
			if constexpr (requires { self.apply(sys); }) {
				self.init(sys);
			}
		}

		void dispatch_apply(this auto&& self, core::SimulationContext & sys) {
			static_assert(
				requires { self.apply(sys); },
				"Controller subclass must implement: void apply(core::SimulationContext &)"
			);
			self.apply(sys);
		}

	private:
		shared::Trigger trigger;
	};



	// define controller concept
	template <class C>
	concept IsController = std::derived_from<C, Controller>;

	// define controller Pack
	template<IsController...>
	struct ControllerPack {};

	// constrained variable template
	template<class... Cs>
	requires (IsController<Cs> && ...)
	inline constexpr ControllerPack<Cs...> controllers {};


	// Concept to check if a type T is a ControllerPack
	template<typename T>
	inline constexpr bool is_controller_pack_v = false; // Default

	template<IsController... Cs>
	inline constexpr bool is_controller_pack_v<ControllerPack<Cs...>> = true; // Specialization

	template<typename T>
	concept IsControllerPack = is_controller_pack_v<std::remove_cvref_t<T>>;
}
