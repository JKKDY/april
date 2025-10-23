#pragma once

#include "april/core/context.h"
#include "april/shared/trigger.h"

namespace april::monitor {


	class Monitor {
	public:
		explicit Monitor(shared::Trigger trig) : trigger(std::move(trig)) {}

		[[nodiscard]] bool should_trigger(const core::SimulationContext & sys) {
			return trigger(sys);
		}

		// Called once at the start to set integration parameters
		void init(const double delta_t, const double start_t, const double end_t, const size_t steps) {
			dt = delta_t;
			start_time = start_t;
			end_time = end_t;
			num_steps = steps;
		}

		// Optional: overridable initialization for custom setup. Is called right after init()
		void dispatch_initialize(this auto&& self) {
			if constexpr (requires { self.initialize(); }) {
				self.initialize();
			}
		}

		// Optional: Called before a step
		void dispatch_before_step(this auto&& self, const core::SimulationContext & context) {
			if constexpr (requires { self.before_step(context); }) {
				self.before_step(context);
			}
		}

		// Required: Called after a step
		void dispatch_record(this auto&& self, const core::SimulationContext & context) {
			static_assert(
				requires { self.record(context); },
				"Monitor subclass must implement: void dispatch_record(const core::SimulationContext &)"
			);
			self.record(context);
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
		shared::Trigger trigger;
	};


	template <class M> concept IsMonitor = std::derived_from<M, Monitor>;


	template<IsMonitor... Ms> struct MonitorPack {};

	template<class... Ms> inline constexpr MonitorPack<Ms...> monitors{};

}