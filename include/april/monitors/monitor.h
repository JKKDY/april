#pragma once

#include "april/env/particle.h"

namespace april::monitor {

	template <typename M> concept IsMonitor = requires(M m,
		const double dt,
		const double start_t,
		const double end_t,
		const size_t step,
		const std::vector<env::impl::ParticleView>& particles) {
			{ m.dispatch_before_step(step, end_t, particles) } -> std::same_as<void>;
			{ m.dispatch_record(step, end_t, particles) } -> std::same_as<void>;
	        { m.call_frequency() } -> std::convertible_to<std::size_t>;
			{ m.init(dt, start_t, end_t, step) } -> std::same_as<void>;
		};


	class Monitor {
	public:
		using Particles = std::vector<env::impl::ParticleView>;
		explicit Monitor(const size_t call_frequency) : call_frequency_m(call_frequency) {}

		[[nodiscard]] size_t call_frequency() const { return call_frequency_m; }

		// Called once at the start
		void init(const double delta_t, const double start_t, const double end_t, const size_t steps) {
			dt = delta_t;
			start_time = start_t;
			end_time = end_t;
			num_steps = steps;
		}

		void dispatch_initialize(this auto&& self) {
			if constexpr (requires { self.initialize(); }) {
				self.initialize();
			}
		}

		// Optional: Called before a step
		void dispatch_before_step(this auto&& self, size_t step, double time, const Particles& particles) {
			if constexpr (requires { self.before_step(step, time, particles); }) {
				self.before_step(step, time, particles);
			}
		}

		// Required: Called after a step
		void dispatch_record(this auto&& self, size_t step, double time, const Particles& particles) {
			static_assert(
				requires { self.record(step, time, particles); },
				"Monitor subclass must implement: void dispatch_record(size_t, double, const Particles&)"
			);
			self.record(step, time, particles);
		}

		// Optional: Called once at the end
		void dispatch_finalize(this auto&& self) {
			if constexpr (requires { self.finalize(); }) {
				self.finalize();
			}
		}

	protected:
		double dt{};
		double start_time{};
		double end_time{};
		size_t num_steps{};
		size_t call_frequency_m;
	};

	template<IsMonitor... Ts> struct MonitorPack {};
	template<class... Ms> inline constexpr MonitorPack<Ms...> monitors{};

}