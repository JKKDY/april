#pragma once

#include <concepts>
#include "april/particle/fields.h"

namespace april::field {

	class Field {
	public:
		// TODO add filters for region, types and ids. Right now dispatch apply is applied to everything

		template<class S>
		void dispatch_init(this auto&& self, const core::SystemContext<S> & sys) {
			if constexpr ( requires { self.init(sys); }) {
				self.init(sys);
			}
		}

		template<class S>
		void dispatch_update(this auto&& self, const core::SystemContext<S> & sys) {
			if constexpr ( requires { self.update(sys); }) {
				self.update(sys);
			}
		}

		template<env::IsUserData U, env::HasFields Self>
		void dispatch_apply(this const Self& self, const env::RestrictedParticleRef<env::FieldOf<Self>, U> & particle) {
			static_assert(
				requires { self.apply(particle); },
				"Field must implement: void apply(env::RestrictedParticleRef<M, U> particle) const"
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