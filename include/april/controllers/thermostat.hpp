#pragma once

#include "april/controllers/controller.hpp"
#include "april/math/distributions.hpp"
namespace april::controller {

	static constexpr double temperature_not_set = -1.0;

	class VelocityScalingThermostat : public Controller {
		static constexpr env::FieldMask mass_vel = env::Field::velocity | env::Field::mass;
		static constexpr env::FieldMask vel = to_field_mask(env::Field::velocity);

	public:
		VelocityScalingThermostat(const double init_T, const double target_T, const double max_dT, const shared::Trigger & trig)
		: Controller(trig),
		 init_temp(init_T),
		 target_temp(target_T),
		 max_temp_change(max_dT) {}


		template<class S>
		void init(core::SystemContext<S> & sys) const {
			AP_ASSERT(sys.size() > 1, "For the thermostat to work correctly, there should be at least two particles");

			if (init_temp == temperature_not_set) return;
			for (size_t i = 0; i < sys.size(); ++i) {
				auto p = sys.template at<mass_vel>(i);
				const double sigma = std::sqrt(init_temp / p.mass);
				p.velocity = shared::maxwell_boltzmann_velocity_distribution(sigma, dimensions(sys));
			}
		}
		template<class S>
		void apply(core::SystemContext<S> & sys) const {
			AP_ASSERT(sys.size() > 1, "For the thermostat to work correctly, there should be at least two particles");

			if (target_temp == temperature_not_set) return;
			const vec3 avg_v = average_velocity<S>(sys);
			const double current_T = temperature(sys, avg_v);
			const double diff = target_temp - current_T;
			const double new_T = current_T + std::clamp(diff, -max_temp_change, max_temp_change);
			if (std::abs(new_T - current_T) < 1e-12) return;

			if (current_T < 1e-12) {
				// Cannot scale from T=0. We must "ignite" the system
				// by setting new thermal velocities, just like init(),
				// but preserving the average velocity.
				for (size_t i = 0; i < sys.size(); ++i) {
					auto p = sys.template at<mass_vel>(i);
					const double sigma = std::sqrt(new_T / p.mass);
					p.velocity = avg_v + shared::maxwell_boltzmann_velocity_distribution(sigma, dimensions(sys));
				}
			} else {
				const double factor = std::sqrt(new_T / current_T);
				scale_thermal_velocities<S>(sys, factor, avg_v);
			}
		}

	private:
		double init_temp;
		double target_temp;
		double max_temp_change;

		template<class S>
		static vec3 average_velocity(const core::SystemContext<S> & sys) {
			if (sys.size() == 0) return {0,0,0};

			vec3 sum{0, 0, 0};
			for (size_t i = 0; i < sys.size(); ++i) {
				auto p = sys.template view<vel>(i);
				sum += p.velocity;
			}
			return sum / static_cast<double>(sys.size());
		}

		template<class S>
		static double temperature(const core::SystemContext<S> & sys, const vec3& avg_v) {
			double kinetic = 0.0;
			for (size_t i = 0; i < sys.size(); ++i) {
				const auto p = sys.template view<mass_vel>(i);
				const vec3 dv = p.velocity - avg_v;
				kinetic += p.mass * dv.norm_squared();
			}

			const size_t D = dimensions(sys);
			if (D == 0) return 0.0;

			const size_t dof = D * sys.size();
			return kinetic / static_cast<double>(dof);
		}

		template<class S>
		static void scale_thermal_velocities(core::SystemContext<S> & sys, const double factor, const vec3& avg_v) {
			for (size_t i = 0; i < sys.size(); ++i) {
				auto p = sys.template at<vel>(i);
				p.velocity = avg_v + factor * (p.velocity - avg_v);
			}
		}

		template<class S>
		static uint8_t dimensions(const core::SystemContext<S> & sys) {
			return 3 -
				(sys.box().extent.x == 0) -
				(sys.box().extent.y == 0) -
				(sys.box().extent.z == 0);
		}
	};

} // namespace april::controller



