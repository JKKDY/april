#pragma once

#include "april/base/types.hpp"
#include "april/particle/fields.hpp"
#include "april/forces/force.hpp"


namespace april::force {
	// No-op force: always returns zero vector and mixes to itself.
	struct NoForce : Force{
		static constexpr env::FieldMask fields = +env::Field::none;

		NoForce(): Force(0) {}


		vec3 eval(auto, auto, const vec3&) const noexcept {
			return vec3{0.0, 0.0, 0.0};
		}

		[[nodiscard]] NoForce mix(NoForce const&) const noexcept {
			return {};
		}

		bool operator==(const NoForce&) const = default;
	};
}