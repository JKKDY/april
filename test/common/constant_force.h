#pragma once

#include "april/interactions/interaction.hpp"

using namespace april;

// A tiny interaction that returns a constant force vector and mixes by summing.
struct ConstantForce final : interaction::Interaction {
	static constexpr auto symmetry = interaction::InteractionSymmetry::Symmetric;
	static constexpr auto fields = ParticleField::none;

	vec3 v;

	ConstantForce(vec3 u, const double cutoff = interaction::no_cutoff)
	: Interaction(cutoff), v{u.x,u.y,u.z} {}

	ConstantForce(const vec3::type x, const vec3::type y, const vec3::type z, const double cutoff = interaction::no_cutoff)
	: Interaction(cutoff), v{x,y,z} {}

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


