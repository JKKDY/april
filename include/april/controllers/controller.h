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
				self.apply(sys);
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
