#pragma once

#include "april/particle/fields.hpp"
#include "april/forces/force.hpp"


namespace april::force {
	// No-op force: always returns zero vector and mixes to itself.
	struct NoForce : Force{
		static constexpr env::FieldMask fields = +env::Field::none;

		NoForce(): Force(0) {}


		auto eval(auto, auto, const auto& r) const noexcept {
			using v = std::remove_cvref_t<decltype(r)>;
			return v{0, 0, 0};
		}

		[[nodiscard]] NoForce mix(NoForce const&) const noexcept {
			return {};
		}

		bool operator==(const NoForce&) const = default;
	};
}