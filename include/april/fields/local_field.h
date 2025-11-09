#pragma once

#include "april/system/context.h"
#include "april/fields/field.h"
#include "april/particle/fields.h"

namespace april::field {
	struct LocalForceField  final : Field {
		static constexpr env::FieldMask fields = env::Field::position | env::Field::force;

		LocalForceField(const vec3 & force_dir, const env::Domain & domain, const double start_time, const double stop_time):
		force(force_dir), region(env::Box::from_domain(domain)), start(start_time), stop(stop_time), active(false) {}


		template<env::IsUserData U>
		void apply(const env::RestrictedParticleRef<fields, U> & particle) const {
			if (!active && region.contains(particle.position))
				particle.force += force;
		}

		template<class S>
		void update(const core::SystemContext<S> & ctx) {
			active = ctx.time() >= start && ctx.time() < stop;
		}

	private:
		const vec3 force;
		env::Box region;
		double start;
		double stop;
		bool active;
	};
}