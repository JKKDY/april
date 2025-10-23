#pragma once
#include <concepts>

#include "april/core/system.h"
#include "april/monitors/monitor.h"
#include "april/common.h"

#include "april/shared/pack_storage.h"

namespace april::integrator {

	template<core::IsSystem Sys, class Pack> class Integrator;  // primary template

	template <core::IsSystem Sys, class... TMonitors>  // partial specialization
	class Integrator<Sys, monitor::MonitorPack<TMonitors...>> {
	public:

		explicit Integrator(Sys& sys_ref)
			: sys(sys_ref)
		{}

		explicit Integrator(Sys& s, monitor::MonitorPack<TMonitors...>) : sys(s) {}

		explicit Integrator(Sys& s, TMonitors... mons) : sys(s) {
			(this->add_monitor(std::move(mons)), ...);
		}

		template<typename T> requires same_as_any<T, TMonitors...>
		void add_monitor(T monitor) {
			// std::get<std::vector<T>>(monitors).push_back(std::move(monitor));
			m.add(std::move(monitor));
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
			return std::forward<decltype(self)>(self);
		}

		template<typename... Ts>
		auto&& with_monitors(this auto&& self, Ts&&... ms) {
			(self.add_monitor(std::forward<Ts>(ms)), ...);
			return std::forward<decltype(self)>(self);
		}

		auto&& run_for(this auto&& self, double dt, double duration) {
			self.run_steps(dt, static_cast<std::size_t>(duration / dt));
			return std::forward<decltype(self)>(self);
		}

		// Integrate for explicit number of steps
		auto&& run_steps(this auto&& self, double delta_t, std::size_t num_steps) {
			if (delta_t <= 0) {
				throw std::invalid_argument(
					"time step cannot be negative. Got delta_t=" + std::to_string(delta_t)
				);
			}

			self.duration = static_cast<double>(num_steps) * delta_t;
			self.dt = delta_t;
			self.num_steps = num_steps;

			self.init_monitors();
			self.dispatch_initialize_monitors();

			for (self.step = 0; self.step < num_steps; ++self.step) {
				self.dispatch_monitor_preparation();
				self.integration_step();
				self.dispatch_monitor_recording();
				self.sys.update_time(self.dt);
				self.sys.increment_step();
			}

			self.finalize_monitors();

			return std::forward<decltype(self)>(self);
		}


	protected:
		Sys & sys;
		size_t num_steps{};
		double duration = 0;
		double dt = 0;
		size_t step = 0;

	private:
		// std::tuple<std::vector<TMonitors>...> monitors{};
		using Monitors = shared::internal::PackStorage<TMonitors...>;
		Monitors m;

		void init_monitors() {
			m.for_each_item([&](auto& mon){mon.init(dt, 0, duration, num_steps); } );
		}

		void dispatch_initialize_monitors() {
			m.for_each_item([&](auto& mon){mon.dispatch_initialize(); } );
		}

		void dispatch_monitor_preparation() {
			m.for_each_item([&](auto & mon) {
				if (mon.should_trigger(sys.context())) {
					mon.dispatch_before_step(sys.context());
				}
			});
		}

		void dispatch_monitor_recording() {
			m.for_each_item([&](auto & mon) {
				if (mon.should_trigger(sys.context())) {
					mon.dispatch_record(sys.context());
				}
			});
		}

		void finalize_monitors() {
			m.for_each_item([&](auto& mon) {mon.dispatch_finalize(); } );
		}
	};


	template <typename T> concept IsIntegrator = requires(T t) {
		{ t.integration_step() } -> std::same_as<void>;
	};

} // namespace april::core