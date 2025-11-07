#pragma once

#include "april/april.h"

using namespace april;
// A tiny force that returns a constant vector and mixes by summing
struct ConstantForce final : force::Force {
	vec3 v;

	ConstantForce(const double x, const double y, const double z, const double cutoff = -1)
	: Force(cutoff), v{x,y,z} {}

	template<env::IsConstFetcher F>
	vec3 operator()(const F&, const F&, const vec3&) const noexcept {
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