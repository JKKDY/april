#pragma once

#include <concepts>
#include "april/particle/scalar_access.hpp"

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

		template<particle::internal::HasFields Self>
		void dispatch_apply(this const Self& self, const auto & particle) {
			// static_assert(
			// 	requires { self.apply(particle); },
			// 	"Field must implement: void apply(auto particle) const"
			// );
			self.apply(particle);
		}
	};

	// define field concept
	template <class FFs>
	concept IsField = std::derived_from<FFs, Field>;

	namespace  internal {
		// define field Pack
		template<IsField...>
		struct FieldPack {};

	}
}

namespace april {
	// constrained variable template
	template<class... FFs>
	requires (field::IsField<FFs> && ...)
	inline constexpr field::internal::FieldPack<FFs...> fields {};

	namespace field::internal {
		// Concept to check if a type T is a ForcePack
		template<typename T>
		inline constexpr bool is_field_pack_v = false; // Default

		template<IsField... FFs>
		inline constexpr bool is_field_pack_v<FieldPack<FFs...>> = true; // Specialization

		template<typename T>
		concept IsFieldPack = is_field_pack_v<std::remove_cvref_t<T>>;
	}

}












