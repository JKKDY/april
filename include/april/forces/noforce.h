#pragma once

#include "april/common.h"
#include "april/particle/fields.h"
#include "april/forces/force.h"


namespace april::force {
	// No-op force: always returns zero vector and mixes to itself.
	struct NoForce : Force{
		static constexpr env::FieldMask fields = +env::Field::none;

		NoForce(): Force(0) {}


		template<env::IsUserData U>
		vec3 eval(const env::ParticleView<fields, U> &, const env::ParticleView<fields, U> &, const vec3&) const noexcept {
			return vec3{0.0, 0.0, 0.0};
		}

		[[nodiscard]] NoForce mix(NoForce const&) const noexcept {
			return {};
		}
	};
}