#pragma once

#include <utility>

#include "april/core/context.hpp"
#include "april/utility/trigger.hpp"

namespace april::monitor {

	/**
	 * @brief Base class for output and diagnostic monitors.
	 *
	 * Custom monitors derive from Monitor and implement:
	 *   void record(const auto& system_context);
	 *
	 * They may optionally implement:
	 *   void initialize();
	 *   void before_step(const auto& system_context);
	 *   void finalize();
	 */
	class Monitor {
	public:
		explicit Monitor(Trigger  trig) : trigger(std::move(trig)) {}

		[[nodiscard]] bool should_trigger(const utility::internal::TriggerContext & sys) const {
			return trigger(sys);
		}

		// Called once at the start to set integration parameters
		void init(
			const double delta_t,
			const double start_t,
			const double end_t,
			const size_t initial_step,
			const size_t steps
		) {
			dt = delta_t;
			start_time = start_t;
			end_time = end_t;
			start_step = initial_step;
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
				"APRIL ERROR: Monitor subclasses must implement record(const core::SystemContext<S>&)."
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
		size_t start_step{};
		size_t num_steps{};
		Trigger trigger;
	};


	template <class M> concept IsMonitor = std::derived_from<std::remove_cvref_t<M>, Monitor>;

	namespace internal {
		template<IsMonitor... Ms> struct MonitorPack {};

		// check if type is monitor pack
		template<typename T>
		inline constexpr bool is_monitor_pack_v = false;

		template<IsMonitor... Ms>
		inline constexpr bool is_monitor_pack_v<MonitorPack<Ms...>> = true;

		template<typename T>
		concept IsMonitorPack = is_monitor_pack_v<std::remove_cvref_t<T>>;

	}
}

namespace april {
	template<class... Ms> inline constexpr monitor::internal::MonitorPack<Ms...> monitors{};
}
















