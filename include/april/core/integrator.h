#pragma once
#include <concepts>

#include "april/env/environment.h"

namespace april::env {
	class Environment;
}

namespace april::core::impl {

	template <typename T>
	concept IsIntegrator = requires(T t) {
			{ t.integration_step() } -> std::same_as<void>;
	};

	// === CRTP Base Integrator ===
	template <typename Derived> class Integrator {
	public:
		explicit Integrator(env::Environment& env_ref)
			: env(env_ref) {
			static_assert(IsIntegrator<Derived>, "Derived must implement: void integration_step()");
		}

		void run(const double dt, const double duration) {
			const std::size_t steps = static_cast<std::size_t>(duration / dt);
			run_steps(dt, steps);
		}

		void run_steps(double dt, const size_t num_steps) {
			this->dt = dt;
			t = 0;

			for (size_t i = 0; i < num_steps; i++) {
				derived()->integration_step();
				t += dt;
			}
		}

	protected:
		Derived* derived() { return static_cast<Derived*>(this); }
		env::Environment& env;
		double t = 0;
		double dt = 0;
	};
} // namespace april::core