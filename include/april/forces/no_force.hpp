#pragma once

#include "april/forces/force.hpp"


namespace april {
	// No-op force: always returns zero vector and mixes to itself.
	struct NoForce : force::Force{
		static constexpr auto fields = ParticleField::none;

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














