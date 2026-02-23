#pragma once

#include "april/base/types.hpp"
#include "april/fields/field.hpp"


namespace april {
	struct UniformField final : field::Field {
		static constexpr auto fields = ParticleField::force;

		explicit UniformField(const vec3 & force_dir): force(force_dir) {}

		template<particle::IsParticleAttributes U>
		void apply(const particle::internal::ScalarRestrictedParticleRef<fields, U> & particle) const {
			particle.force += force;
		}

	private:
		const vec3 force;
	};
}












