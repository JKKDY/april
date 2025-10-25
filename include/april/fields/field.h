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
			if constexpr ( requires { self.init(sys); }) {
				self.init(sys);
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

	// define controller concept
	template <class FFs>
	concept IsField = std::derived_from<FFs, Field>;

	// define controller Pack
	template<IsField...>
	struct FieldPack {};

	// constrained variable template
	template<class... FFs>
	requires (IsField<FFs> && ...)
	inline constexpr FieldPack<FFs...> fields {};
}