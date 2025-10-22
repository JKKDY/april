#pragma once

#include <cstddef>
#include <functional>

#include "april/core/context.h"

namespace april::shared {
	struct Trigger {
		using TriggerFn = std::move_only_function<bool(const core::SimulationContext&)>;

		bool operator()(const core::SimulationContext& sys) {
			return fn_(sys);
		}

		// ---- convenience constructors ----

		// step based triggers:
		// trigger every N steps
		static Trigger every(const size_t N, const size_t offset = 0) {
			return Trigger{[=](const core::SimulationContext&sys) { return (sys.step() + offset) % N == 0; }};
		}

		// start triggering after *step*
		static Trigger after(const size_t step) {
			return Trigger{[=](const core::SimulationContext& sys) { return sys.step() >= step; }};
		}

		// trigger only while current step in [start, end)
		static Trigger between(const size_t start, const size_t end) {
			return Trigger{[=](const core::SimulationContext& sys) {
				const auto s = sys.step();
				return s >= start && s < end;
			}};
		}

		/// Fires exactly at a single step
		static Trigger at_step(const std::size_t step) {
			return Trigger{[=](const core::SimulationContext& ctx) {
				return ctx.step() == step;
			}};
		}



		// time based triggers:
		// trigger after given time period
		static Trigger periodically(const double period, const double offset = 0.0) {
			return Trigger{[=, last = offset - period](const core::SimulationContext& sys) mutable {
				if (sys.time() - last >= period) {
					last = sys.time();
					return true;
				}
				return false;
			}};
		}

		// start triggering only after a given time t
		static Trigger after_time(const double t) {
			return Trigger{[=](const core::SimulationContext& sys) {
				return sys.time() >= t;
			}};
		}

		// trigger only while simulation time between t_start and t_end
		static Trigger between_time(const double t_start, const double t_end) {
			return Trigger{[=](const core::SimulationContext& sys) {
				const auto t = sys.time();
				return t >= t_start && t < t_end;
			}};
		}

		// generic triggers:
		// trigger every step
		static Trigger always() {
			return Trigger{[](const core::SimulationContext&) { return true; }};
		}
		// for custom triggers
		static Trigger when(std::move_only_function<bool(const core::SimulationContext&)> pred) {
			return Trigger{std::move(pred)};
		}



		// ---- Chaining operators ----
		// Logical AND
		friend Trigger operator&&(Trigger lhs, Trigger rhs) {
			return Trigger{[lhs = std::move(lhs), rhs = std::move(rhs)]
						   (const core::SimulationContext& sys) mutable {
				return lhs(sys) && rhs(sys);
			}};
		}

		// Logical OR
		friend Trigger operator||(Trigger lhs, Trigger rhs) {
			return Trigger{[lhs = std::move(lhs), rhs = std::move(rhs)]
						   (const core::SimulationContext& sys) mutable {
				return lhs(sys) || rhs(sys);
			}};
		}

		// Logical NOT
		friend Trigger operator!(Trigger t) {
			return Trigger{[t = std::move(t)](const core::SimulationContext& sys) mutable {
				return !t(sys);
			}};
		}

	private:
		explicit Trigger(TriggerFn fn) : fn_(std::move(fn)) {}

		 TriggerFn fn_;
	};
}
