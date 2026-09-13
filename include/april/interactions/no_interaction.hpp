#pragma once

#include "april/interactions/interaction.hpp"


namespace april {
	// No-op interaction: always returns a zero vector and mixes to itself.
	struct NoInteraction : interaction::Interaction{
		static constexpr auto fields = ParticleField::none;

		NoInteraction(): Interaction(0) {}


		auto eval(auto, auto, const auto& r) const noexcept {
			using v = std::remove_cvref_t<decltype(r)>;
			return v{0, 0, 0};
		}

		[[nodiscard]] NoInteraction mix(NoInteraction const&) const noexcept {
			return {};
		}

		bool operator==(const NoInteraction&) const = default;
	};
}














