#pragma once

#include "april/shared/trigger.h"
#include "april/core/context.h"

namespace april::controller  {

	template<trigger::IsTrigger Trig>
	class Controller {
	public:
		explicit Controller(const Trig & trig) : trigger(trig) {}

		template<class S>
		[[nodiscard]] bool should_trigger(const core::SystemContext<S> & sys) {
			return trigger(sys);
		}

		template<class S>
		void dispatch_init(this auto && self, const core::SystemContext<S> & sys) {
			if constexpr (requires { self.apply(sys); }) {
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
		Trig trigger;
	};



	template<typename C>
	struct has_controller_base {
	private:
		template<trigger::IsTrigger Trig>
		static std::true_type test(const Controller<Trig>*);
		static std::false_type test(...);
	public:
		static constexpr bool value = decltype(test(std::declval<C*>()))::value;
	};

	template <class C>
	concept IsController = has_controller_base<C>::value;

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
