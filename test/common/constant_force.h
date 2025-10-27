#pragma once

#include "april/april.h"

using namespace april;
// A tiny force that returns a constant vector and mixes by summing
struct ConstantForce final : force::Force {
	vec3 v;

	ConstantForce(double x, double y, double z, double cutoff = -1)
	: Force(cutoff), v{x,y,z} {}

	[[nodiscard]] vec3 eval(const env::internal::Particle&, const env::internal::Particle&, const vec3&) const noexcept {
		return v;
	}

	[[nodiscard]] ConstantForce mix(const ConstantForce& other) const noexcept {
		return {
			v.x + other.v.x,
			v.y + other.v.y,
			v.z + other.v.z,
			std::max(cutoff, other.cutoff)
		};
	}
};