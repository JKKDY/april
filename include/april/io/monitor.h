#pragma once

#include "april/env/particle.h"

namespace april::io {

	template <typename M> concept IsMonitor = requires(M m,
		const double dt,
		const double start_t,
		const double end_t,
		const size_t step,
		const std::vector<env::impl::Particle>& particles) {
			{ m.record(step, end_t, particles) } -> std::same_as<void>;
	        { m.call_frequency() } -> std::convertible_to<std::size_t>;
			{ m.init(dt, start_t, end_t, step) } -> std::same_as<void>;
		};


	class Monitor {
	public:
		explicit Monitor(const size_t call_frequency) : call_frequency_m(call_frequency) {}

		[[nodiscard]] size_t call_frequency() const { return call_frequency_m; }

		void init(const double dt, const double start_t, const double end_t, const size_t num_steps) {
			this->dt = dt;
			this->start_time = start_t;
			this->end_time = end_t;
			this->num_steps = num_steps;
		}

		void record(this auto&& self, size_t step, double time, const std::vector<env::impl::Particle>& particles) {
			static_assert(
				requires { self.write_output(step, time, particles); },
				"OutputWriter requires a write_output(size_t, const std::vector<Particle>&) method"
			);
			self.write_output(step, time, particles);
		}

	protected:
		double dt{};
		double start_time{};
		double end_time{};
		size_t num_steps{};
		size_t call_frequency_m;
	};
}