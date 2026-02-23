#pragma once

#include "april/base/types.hpp"
#include "april/fields/field.hpp"


namespace april::field {
	struct UniformField final : Field {
		static constexpr ParticleField fields = ParticleField::force;

		explicit UniformField(const vec3 & force_dir): force(force_dir) {}

		template<core::IsUserData U>
		void apply(const core::ScalarRestrictedParticleRef<fields, U> & particle) const {
			particle.force += force;
		}

	private:
		const vec3 force;
	};
}




