#pragma once

#include <concepts>
#include "april/env/particle.h"

namespace april::field {

	class Field {
	public:
		void dispatch_init(this auto&& self, const core::SimulationContext & sys) {
			if constexpr ( requires { self.init(sys); }) {
				self.init(sys);
			}
		}

		void dispatch_update(this auto&& self, const core::SimulationContext & sys) {
			if constexpr ( requires { self.update(sys); }) {
				self.update(sys);
			}
		}

		void dispatch_apply(this const auto& self, env::RestrictedParticleRef particle) {
			static_assert(
				requires { self.apply(particle); },
				"Field must implement: void apply(env::RestrictedParticleRef particle) const"
			);
			self.apply(particle);
		}
	};

	// define field concept
	template <class FFs>
	concept IsField = std::derived_from<FFs, Field>;

	// define field Pack
	template<IsField...>
	struct FieldPack {};

	// constrained variable template
	template<class... FFs>
	requires (IsField<FFs> && ...)
	inline constexpr FieldPack<FFs...> fields {};

	// Concept to check if a type T is a ForcePack
	template<typename T>
	inline constexpr bool is_field_pack_v = false; // Default

	template<IsField... FFs>
	inline constexpr bool is_field_pack_v<FieldPack<FFs...>> = true; // Specialization

	template<typename T>
	concept IsFieldPack = is_field_pack_v<std::remove_cvref_t<T>>;
}