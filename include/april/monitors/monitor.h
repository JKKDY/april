#pragma once

#include "april/core/context.h"
#include "april/shared/trigger.h"

namespace april::monitor {


	template<trigger::IsTrigger Trig>
	class Monitor {
	public:
		explicit Monitor(const Trig & trig) : trigger(trig) {}

		template<class S>
		[[nodiscard]] bool should_trigger(const core::SystemContext<S> & sys) {
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
		template<class S>
		void dispatch_before_step(this auto&& self, const core::SystemContext<S> & sys) {
			if constexpr (requires { self.before_step(sys); }) {
				self.before_step(sys);
			}
		}

		// Required: Called after a step
		template<class S>
		void dispatch_record(this auto&& self, const core::SystemContext<S> & sys) {
			static_assert(
				requires { self.record(sys); },
				"Monitor subclass must implement: void dispatch_record(const core::SimulationContext &)"
			);
			self.record(sys);
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
		Trig trigger;
	};


	// define monitor concept
	template <typename T>
	struct has_monitor_base {
	private:
		template <trigger::IsTrigger Trig>
		static std::true_type test(const Monitor<Trig>*); // use ptr conversion from T* to Monitor*
		static std::false_type test(...);
	public:
		static constexpr bool value = decltype(test(std::declval<T*>()))::value;
	};

	template <typename T>
	concept IsMonitor = has_monitor_base<T>::value;

	template<IsMonitor... Ms> struct MonitorPack {};

	template<class... Ms> inline constexpr MonitorPack<Ms...> monitors{};

}