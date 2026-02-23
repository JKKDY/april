#pragma once

#include "april/base/types.hpp"
#include "april/core/context.hpp"
#include "april/fields/field.hpp"


namespace april::field {
	struct LocalForceField  final : Field {
		static constexpr ParticleField fields = ParticleField::position | ParticleField::force;

		LocalForceField(const vec3 & force_dir, const Domain & domain, const double start_time, const double stop_time):
		force(force_dir), region(core::Box::from_domain(domain)), start(start_time), stop(stop_time), active(start_time == 0) {}


		template<particle::IsParticleAttributes U>
		void apply(const particle::internal::ScalarRestrictedParticleRef<fields, U> & particle) const {
			if (active && region.contains(particle.position))
				particle.force += force;
		}

		template<class S>
		void update(const core::SystemContext<S> & sys) {
			active = sys.time() >= start && sys.time() < stop;
		}

	private:
		const vec3 force;
		core::Box region;
		double start;
		double stop;
		bool active;
	};
}












