#pragma once

#include "april/core/context.h"
#include "april/fields/field.h"
#include "april/env/particle.h"

namespace april::field {
	struct LocalForceField  final : Field {

		LocalForceField(const vec3 & force_dir, const double start_time, const double stop_time):
		force(force_dir), start(start_time), stop(stop_time), active(false) {}

		void apply(env::RestrictedParticleRef particle) const {
			if (!active) return;
			particle.force += force;
		}

		void update(const core::SimulationContext & ctx) {
			active = ctx.time() >= start && ctx.time() < stop;
		}

	private:
		const vec3 force;
		double start;
		double stop;
		bool active;
	};
}