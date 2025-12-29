#pragma once

#include "april/april.h"

using namespace april;
// A tiny force that returns a constant vector and mixes by summing
struct ConstantForce final : force::Force {
	static constexpr env::FieldMask fields = +env::Field::none;

	vec3 v;

	ConstantForce(const double x, const double y, const double z, const double cutoff = force::no_cutoff)
	: Force(cutoff), v{x,y,z} {}

	template<env::FieldMask M, env::IsUserData U>
	vec3 operator()(const ParticleView<M, U> &, const ParticleView<M, U> &, const vec3&) const noexcept {
		return v;
	}

	bool operator==(const ConstantForce&) const = default;

	[[nodiscard]] ConstantForce mix(const ConstantForce& other) const noexcept {
		return {
			v.x + other.v.x,
			v.y + other.v.y,
			v.z + other.v.z,
			std::max(cutoff(), other.cutoff())
		};
	}
};