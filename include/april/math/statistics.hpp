#pragma once

#include <random>

#include "april/common.hpp"

// we use a constant seed for repeatability.
// random engine needs static lifetime otherwise it would be recreated for every call.
static std::default_random_engine randomEngine(42);

namespace april::shared {

	inline vec3 maxwell_boltzmann_velocity_distribution(const double averageVelocity, const size_t dimensions) {
		// when adding independent normally distributed values to all velocity components
		// the velocity change is maxwell boltzmann distributed
		std::normal_distribution<double> normalDistribution{0, 1};
		vec3 randomVelocity{};
		for (size_t  i = 0; i < dimensions; ++i) {
			randomVelocity[static_cast<int>(i)] = averageVelocity * normalDistribution(randomEngine);
		}
		return randomVelocity;
	}

}