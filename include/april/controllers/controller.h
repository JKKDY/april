#pragma once

#include <utility>

#include "april/shared/trigger.h"
#include "april/system/context.h"

namespace april::controller  {

	class Controller {
	public:
		explicit Controller(shared::Trigger trig) : trigger(std::move(trig)) {}

		[[nodiscard]] bool should_trigger(const shared::TriggerContext & sys) const {
			return trigger(sys);
		}

		// TODO add static constexpr bool override flags so the user can explicitly state their intent to implement a CRTP method
		template<class S>
		void dispatch_init(this auto && self, core::SystemContext<S> & sys) {
			if constexpr (requires { self.init(sys); }) {
				self.template init<S>(sys);
			}
		}

		template<class S>
		void dispatch_apply(this auto&& self, core::SystemContext<S> & sys) {
			static_assert(
				requires { self.apply(sys); },
				"Controller subclass must implement: void apply(const core::SystemContext<S> & sys)"
			);
			self.template apply<S>(sys);
		}

	private:
		shared::Trigger trigger;
	};




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

} // namespace april::controller
