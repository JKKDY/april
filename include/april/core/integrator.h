#pragma once
#include <concepts>

#include "april/env/environment.h"
#include "april/io/output.h"
#include "april/io/monitor.h"

namespace april::env {
	class Environment;
}

namespace april::core::impl {

	template <typename T> concept IsIntegrator = requires(T t) {
			{ t.integration_step() } -> std::same_as<void>;
	};

	template <typename T, typename... Ts> concept same_as_any = (... or std::same_as<T, Ts>);

	template <io::IsMonitor ... TMonitors>
	class Integrator {
	public:
		explicit Integrator(env::Environment& env_ref)
			: env(env_ref) {
			env.build();
		}

		template<typename T> requires same_as_any<T, TMonitors...>
		void add_monitor(T monitor) {
			std::get<std::vector<T>>(monitors).push_back(std::move(monitor));
		}

		// Call with total duration
		void run(this auto&& self, double dt, const double duration) {
			self.duration = duration;
			const auto steps = static_cast<std::size_t>(duration / dt);
			self.run_steps(dt, steps);
		}

		// Call with explicit number of steps
		void run_steps(this auto&& self, const double delta_t, const std::size_t num_steps) {
			self.dt = delta_t;
			self.time = 0;
			self.num_steps = num_steps;

			for (self.step = 0; self.step < num_steps; ++self.step) {
				self.integration_step();
				self.time += self.dt;

				self.dispatch_monitors();
				self.dispatch_controllers();
			}
		}

	protected:
		env::Environment& env;
		size_t num_steps;
		double duration = 0;
		double time = 0;
		double dt = 0;
		size_t step = 0;

	private:
		std::tuple<std::vector<TMonitors>...> monitors{};

		void init_monitors() {
			auto init  = [&] (auto &lst) {
				for (auto &monitor : lst) {
					if (step % monitor.call_frequency() == 0) {
						monitor.init(dt, 0, duration, num_steps);
					}
				}
			};

			std::apply(
			   [init](auto &... lst) {
				 (init(lst), ...);
				},
				monitors);
		}

		void dispatch_monitors() {
			auto dispatch  = [&] (auto &lst) {
				for (auto &monitor : lst) {
					if (step % monitor.call_frequency() == 0) {
						monitor.record(step, time, env.export_particles());
					}
				}
			};

			std::apply(
			   [dispatch](auto &... lst) {
				 (dispatch(lst), ...);
				},
				monitors);
		}

		void dispatch_controllers() {

		}
	};

} // namespace april::core