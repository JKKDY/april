#pragma once
#include <concepts>

#include "april/core/system.h"
#include "april/io/monitor.h"
#include "april/common.h"


namespace april::core::impl {

	template <typename T> concept IsIntegrator = requires(T t) {
			{ t.integration_step() } -> std::same_as<void>;
	};

	template<IsSystem Sys, class Pack> class Integrator;  // primary template

	template <IsSystem Sys, class... TMonitors>  // partial specialization
	class Integrator<Sys, io::MonitorPack<TMonitors...>> {
	public:
		explicit Integrator(Sys& sys_ref)
			: sys(sys_ref)
		{}

		explicit Integrator(Sys& s, io::MonitorPack<TMonitors...>) : sys(s) {}

		explicit Integrator(Sys& s, TMonitors... mons) : sys(s) {
			(this->add_monitor(std::move(mons)), ...);
		}

		template<typename T> requires same_as_any<T, TMonitors...>
		void add_monitor(T monitor) {
			std::get<std::vector<T>>(monitors).push_back(std::move(monitor));
		}

		// Call with total duration
		void run_for(this auto&& self, double dt, const double duration) {
			self.run_steps(dt, static_cast<std::size_t>(duration / dt));
		}

		// Call with explicit number of steps
		void run_steps(this auto&& self, const double delta_t, const std::size_t num_steps) {
			self.duration = static_cast<double>(num_steps) * delta_t;
			self.dt = delta_t;
			self.time = 0;
			self.num_steps = num_steps;

			self.init_monitors();
			self.dispatch_initialize_monitors();

			for (self.step = 0; self.step < num_steps; ++self.step, self.time+=self.dt) {
				self.dispatch_monitor_preparation();

				self.integration_step();

				self.dispatch_monitor_recording();
			}

			self.finalize_monitors();
		}

	protected:
		Sys & sys;
		size_t num_steps{};
		double duration = 0;
		double time = 0;
		double dt = 0;
		size_t step = 0;

	private:
		std::tuple<std::vector<TMonitors>...> monitors{};

		template<typename Func> void for_each_monitor(Func&& f) {
			std::apply(
				[&](auto&... lst) {
					(f(lst), ...);
				},
				monitors
			);
		}

		void init_monitors() {
			for_each_monitor([&](auto& list) {
				for (auto& monitor : list) {
					monitor.init(dt, 0, duration, num_steps);
				}
			});
		}

		void dispatch_initialize_monitors() {
			for_each_monitor([&](auto& list) {
				for (auto& monitor : list) {
					monitor.dispatch_initialize();
				}
			});
		}

		void dispatch_monitor_preparation() {
			for_each_monitor([&](auto& list) {
				for (auto& monitor : list) {
					if (step % monitor.call_frequency() == 0) {
						monitor.dispatch_before_step(step, time, sys.export_particles());
					}
				}
			});
		}

		void dispatch_monitor_recording() {
			for_each_monitor([&](auto& list) {
				for (auto& monitor : list) {
					if (step % monitor.call_frequency() == 0) {
						monitor.dispatch_record(step, time, sys.export_particles());
					}
				}
			});
		}

		void finalize_monitors() {
			for_each_monitor([&](auto& list) {
				for (auto& monitor : list) {
					monitor.dispatch_finalize();
				}
			});
		}
	};
} // namespace april::core