#pragma once
#include <concepts>

#include "april/system/system.hpp"
#include "april/monitors/monitor.hpp"
#include "april/common.hpp"

#include "april/shared/pack_storage.hpp"

namespace april::integrator {

	template<core::IsSystem Sys, class Pack> class Integrator;  // primary template

	template <core::IsSystem Sys, class... TMonitors>  // partial specialization
	class Integrator<Sys, monitor::MonitorPack<TMonitors...>> {
	public:

		explicit Integrator(Sys& sys_ref)
			: sys(sys_ref)
		{}

		explicit Integrator(Sys& s, monitor::MonitorPack<TMonitors...>) : sys(s) {}

		template<typename T> requires same_as_any<T, TMonitors...>
		void add_monitor(T monitor) {
			monitors.add(std::move(monitor));
		}

		template<typename... Ts>
		void add_monitors(Ts&&... ms) {
			(add_monitor(std::forward<Ts>(ms)), ...);
		}

		// DSL-style chaining
		template<typename T>
		requires same_as_any<T, TMonitors...>
		auto&& with_monitor(this auto&& self, T monitor) {
			self.add_monitor(std::move(monitor));
			return self;
		}

		template<typename... Ts>
		auto&& with_monitors(this auto&& self, Ts&&... ms) {
			(self.add_monitor(std::forward<Ts>(ms)), ...);
			return self;
		}


		void set_dt(const double delta_t) {
			dt = delta_t;
		}

		void set_duration(const double dur) {
			duration = dur;
			most_recently_set = DURATION;
		}

		void set_steps(const size_t steps) {
			num_steps = steps;
			most_recently_set = STEP;
		}

		auto&& with_dt(this auto&& self, const double delta_t) {
			self.set_dt(delta_t);
			return self;
		}

		auto&& for_duration(this auto&& self, const double duration) {
			self.set_duration(duration);
			return self;
		}

		auto&& for_steps(this auto&& self, size_t num_steps) {
			self.set_steps(num_steps);
			return self;
		}

		auto&& run(this auto&& self) {
			if (self.dt <= 0) {
				throw std::invalid_argument(
					"time step cannot be negative. Got delta_t=" + std::to_string(self.dt)
				);
			}

			if (self.most_recently_set == DURATION) {
				self.num_steps = static_cast<size_t>(self.duration / self.dt);
			} else if (self.most_recently_set == STEP) {
				self.duration = static_cast<double>(self.num_steps) * self.dt;
			} else {
				throw std::invalid_argument("neither duration nor number steps have been specified!");
			}

			self.init_monitors();
			self.dispatch_initialize_monitors();

			// simulation loop
			self.sys.update_forces(); // ensure valid force initialization
			for (self.step = 0; self.step < self.num_steps; ++self.step) {
				self.dispatch_monitor_preparation();
				self.integration_step();
				self.dispatch_monitor_recording();
				self.sys.update_time(self.dt);
				self.sys.increment_step();
			}

			self.finalize_monitors();

			return self;
		}

		auto&& run_for_duration(this auto&& self, double delta_t, const double duration) {
			return self.with_dt(delta_t).for_duration(duration).run();
		}

		// Integrate for explicit number of steps
		auto&& run_for_steps(this auto&& self, double delta_t, size_t num_steps) {
			return self.with_dt(delta_t).for_steps(num_steps).run();
		}


	protected:
		Sys & sys;
		size_t num_steps{};
		double duration = 0;
		double dt = 0;
		size_t step = 0;

	private:
		enum MostRecentlySet {
			DURATION = 1,
			NONE = 0,
			STEP = -1
		};
		MostRecentlySet most_recently_set = NONE;

		shared::internal::PackStorage<TMonitors...> monitors;

		void init_monitors() {
			monitors.for_each_item([&](auto& mon){mon.init(dt, 0, duration, num_steps); } );
		}

		void dispatch_initialize_monitors() {
			monitors.for_each_item([&](auto& mon){mon.dispatch_initialize(); } );
		}

		void dispatch_monitor_preparation() {
			monitors.for_each_item([&](auto & mon) {
				if (mon.should_trigger(sys.trigger_context())) {
					mon.dispatch_before_step(sys.context());
				}
			});
		}

		void dispatch_monitor_recording() {
			monitors.for_each_item([&](auto & mon) {
				if (mon.should_trigger(sys.trigger_context())) {
					mon.dispatch_record(sys.context());
				}
			});
		}

		void finalize_monitors() {
			monitors.for_each_item([&](auto& mon) {mon.dispatch_finalize(); } );
		}
	};


	template <typename T> concept IsIntegrator = requires(T t) {
		{ t.integration_step() } -> std::same_as<void>;
	};

} // namespace april::core