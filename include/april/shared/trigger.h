#pragma once

#include <cstddef>

#include "april/system/context.h"

namespace april::shared {

	struct TriggerContext {
		virtual ~TriggerContext() = default;

		[[nodiscard]] virtual env::Box box() const noexcept = 0;
		[[nodiscard]] virtual double time() const noexcept = 0;
		[[nodiscard]] virtual size_t step() const noexcept = 0;
		[[nodiscard]] virtual size_t size() const noexcept = 0;
		[[nodiscard]] virtual size_t size(env::ParticleState state) const noexcept = 0;

		// Todo: add more accessors
	};

	template<class System>
	struct TriggerContextImpl final : TriggerContext {
		explicit TriggerContextImpl(const System & sys): system(sys) {}

		[[nodiscard]] env::Box box() const noexcept override { return system.box(); }
		[[nodiscard]] double time() const noexcept override { return system.time(); }
		[[nodiscard]] size_t step() const noexcept override { return system.step(); }
		[[nodiscard]] size_t size() const noexcept override {
			return system.size(env::ParticleState::ALL);
		}
		[[nodiscard]] size_t size(env::ParticleState state) const noexcept override {
			return system.size(state);
		}

	private:
		const System& system;
	};

	struct Trigger {
		using TriggerFn = std::function<bool(const TriggerContext&)>;

		bool operator()(const TriggerContext& sys) const {
			return fn_(sys);
		}

		explicit Trigger(TriggerFn fn) : fn_(std::move(fn)) {}

		// ---- convenience constructors ----

		// step based triggers:
		// trigger every N steps
		static Trigger every(const size_t N, const size_t offset = 0) {
			return Trigger{[=](const TriggerContext&sys) { return (sys.step() + offset) % N == 0; }};
		}

		// start triggering after *step*
		static Trigger after(const size_t step) {
			return Trigger{[=](const TriggerContext& sys) { return sys.step() >= step; }};
		}

		// trigger only while current step in [start, end)
		static Trigger between(const size_t start, const size_t end) {
			return Trigger{[=](const TriggerContext& sys) {
				const auto s = sys.step();
				return s >= start && s < end;
			}};
		}

		/// Fires exactly at a single step
		static Trigger at_step(const std::size_t step) {
			return Trigger{[=](const TriggerContext& ctx) {
				return ctx.step() == step;
			}};
		}



		// time based triggers:
		// trigger after given time period
		static Trigger periodically(const double period, const double offset = 0.0) {
			return Trigger{[=, last = offset - period](const TriggerContext& sys) mutable {
				if (sys.time() - last >= period) {
					last = sys.time();
					return true;
				}
				return false;
			}};
		}

		// start triggering only after a given time t
		static Trigger after_time(const double t) {
			return Trigger{[=](const TriggerContext& sys) {
				return sys.time() >= t;
			}};
		}

		// trigger only while simulation time between t_start and t_end
		static Trigger between_time(const double t_start, const double t_end) {
			return Trigger{[=](const TriggerContext& sys) {
				const auto t = sys.time();
				return t >= t_start && t < t_end;
			}};
		}

		// generic triggers:
		// trigger every step
		static Trigger always() {
			return Trigger{[](const TriggerContext&) { return true; }};
		}
		// for custom triggers
		static Trigger when(const std::function<bool(const TriggerContext&)>& pred) {
			return Trigger{pred};
		}



		// ---- Chaining operators ----
		// Logical AND
		friend Trigger operator&&(Trigger lhs, Trigger rhs) {
			return Trigger{[lhs = std::move(lhs), rhs = std::move(rhs)]
						   (const TriggerContext& sys) mutable {
				return lhs(sys) && rhs(sys);
			}};
		}

		// Logical OR
		friend Trigger operator||(Trigger lhs, Trigger rhs) {
			return Trigger{[lhs = std::move(lhs), rhs = std::move(rhs)]
						   (const TriggerContext& sys) mutable {
				return lhs(sys) || rhs(sys);
			}};
		}

		// Logical NOT
		friend Trigger operator!(Trigger t) {
			return Trigger{[t = std::move(t)](const TriggerContext& sys) mutable {
				return !t(sys);
			}};
		}

	private:
		TriggerFn fn_;
	};
}
