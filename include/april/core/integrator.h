#pragma once
#include <concepts>

#include "april/env/environment.h"
#include "april/io/output.h"
#include "april/io/monitor.h"

namespace april::env {
	class Environment;
}

namespace april::core::impl {

	template <typename T>
	concept IsIntegrator = requires(T t) {
			{ t.integration_step() } -> std::same_as<void>;
	};

	template <io::IsOutputWriter OutputW = io::NullOutput>
	class Integrator {
	public:
		explicit Integrator(env::Environment& env_ref)
			: env(env_ref) {
			env.build();
		}

		void set_writer(OutputW writer) {
			output_writer = writer;
		}

		// Call with total duration
		void run(this auto&& self, double dt, const double duration) {
			const auto steps = static_cast<std::size_t>(duration / dt);
			self.run_steps(dt, steps);
		}

		// Call with explicit number of steps
		void run_steps(this auto&& self, const double delta_t, const std::size_t num_steps) {
			self.dt = delta_t;
			self.t = 0;

			for (std::size_t i = 0, j = 0; i < num_steps; ++i, ++j) {
				self.integration_step();
				self.t += self.dt;

				if (++j >= self.output_writer.write_frequency) {
					j = 0;
					self.output_writer.write(i, self.env.export_particles());
				}
			}
		}

	protected:
		OutputW output_writer;
		env::Environment& env;
		double t = 0;
		double dt = 0;
	};

} // namespace april::core