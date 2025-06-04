#pragma once

#include "april/env/particle.h"

namespace april::io {

	template <typename M> concept IsMonitor = requires(M m, size_t step, double time,
		const std::vector<env::impl::Particle>& particles) {
		{ m.record(step, time, particles) } -> std::same_as<void>;
        { m.call_frequency() } -> std::convertible_to<std::size_t>;	};


	class Monitor {
	public:
		explicit Monitor(const size_t call_frequency) : call_frequency_m(call_frequency) {}

		[[nodiscard]] size_t call_frequency() const { return call_frequency_m; }

		void record(this auto&& self, size_t step, double time, const std::vector<env::impl::Particle>& particles) {
			static_assert(
				requires { self.write_output(step, time, particles); },
				"OutputWriter requires a write_output(size_t, const std::vector<Particle>&) method"
			);
			self.write_output(step, time, particles); // works if 'self' has a write(t, particles) method
		}

	private:
		size_t call_frequency_m;
	};
}