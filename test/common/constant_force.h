#pragma once

#include "april/forces/force.hpp"

using namespace april;
using namespace april::core;

// A tiny force that returns a constant vector and mixes by summing
struct ConstantForce final : force::Force {
	static constexpr ParticleField fields = ParticleField::none;

	vec3 v;

	ConstantForce(const vec3::type x, const vec3::type y, const vec3::type z, const double cutoff = force::no_cutoff)
	: Force(cutoff), v{x,y,z} {}

	auto operator()(auto, auto, auto) const noexcept {
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


