#pragma once
#include <random>
#include "april/base/types.hpp"

namespace april::math {

	// Helper to get a thread-safe, consistent engine
	inline std::default_random_engine& get_random_engine() {
		thread_local std::default_random_engine engine(42);
		return engine;
	}

	/**
	 * Generates a velocity vector according to the Maxwell-Boltzmann distribution.
	 * @param averageVelocity Scaling factor (related to sqrt(kT/m))
	 * @param dimensions Number of active dimensions (1, 2, or 3)
	 */
	inline vec3d maxwell_boltzmann_velocity(const double averageVelocity, const size_t dimensions = 3) {
		std::normal_distribution dist{0.0, 1.0};
		auto& engine = get_random_engine();

		vec3d velocity{0.0};
		for (size_t i = 0; i < dimensions; ++i) {
			velocity[static_cast<int>(i)] = averageVelocity * dist(engine);
		}
		return velocity;
	}

}