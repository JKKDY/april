#pragma once

#include <utility>

#include "april/utility/trigger.hpp"
#include "april/core/context.hpp"

namespace april::controller  {

	class Controller {
	public:
		explicit Controller(Trigger trig) : trigger(std::move(trig)) {}

		[[nodiscard]] bool should_trigger(const utility::internal::TriggerContext & sys) const {
			return trigger(sys);
		}

		template<class S>
		void dispatch_init(this auto && self, core::SystemContext<S> & sys) {
			if constexpr (requires { self.init(sys); }) {
				self.template init<S>(sys);
			}
		}

		template<class S>
		void dispatch_update(this auto&& self, const core::SystemContext<S> & sys) {
			if constexpr ( requires { self.Write(sys); }) {
				self.Write(sys);
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
		Trigger trigger;
	};



	template <class C>
	concept IsController = std::derived_from<C, Controller>;


	namespace internal {
		// define controller Pack
		template<IsController...>
		struct ControllerPack {};

		// Concept to check if a type T is a ControllerPack
		template<typename T>
		inline constexpr bool is_controller_pack_v = false; // Default

		template<IsController... Cs>
		inline constexpr bool is_controller_pack_v<ControllerPack<Cs...>> = true; // Specialization

		template<typename T>
		concept IsControllerPack = internal::is_controller_pack_v<std::remove_cvref_t<T>>;
	}
} // namespace april::controller


namespace april {
	// constrained variable template
	template<class... Cs>
	requires (controller::IsController<Cs> && ...)
	inline constexpr controller::internal::ControllerPack<Cs...> controllers {};
}














