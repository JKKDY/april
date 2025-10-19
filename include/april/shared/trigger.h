#pragma once

#include <cstddef>
#include <functional>

#include "april/core/system.h"

namespace april {
	struct Trigger {
		using Fn = std::function<bool(size_t step, double t, const core::System&)>;

		explicit Trigger(Fn fn) : fn_(std::move(fn)) {}

		bool operator()(std::size_t step, double t, const core::System& sys) const {
			return fn_ ? fn_(step, t, sys) : false;
		}

		// Some convenience constructors:
		static Trigger every(std::size_t N) {
			return Trigger{[=](std::size_t s, double, const core::System&) { return s % N == 0; }};
		}
		static Trigger always() {
			return Trigger{[](auto..., const core::System<>&) { return true; }};
		}

	private:
		Fn fn_;
	};

}
