#pragma once

#include "april/controllers/controller.h"
#include "april/shared/distributions.h"
namespace april::controller {

	static constexpr double TemperatureNotSet = -1.0;

	class VelocityScalingThermostat : public Controller {
	public:
		VelocityScalingThermostat(const double init_T, const double target_T, const double max_dT, shared::Trigger trig)
	   : Controller(std::move(trig)),
		 init_temp(init_T),
		 target_temp(target_T),
		 max_temp_change(max_dT) {}

		void init(core::SimulationContext & ctx) const {
			if (init_temp == TemperatureNotSet) return;

			for (size_t i = ctx.index_start(); i < ctx.index_end(); ++i) {
				auto p = ctx.get_particle_by_index(i);
				const double sigma = std::sqrt(init_temp / p.mass);
				p.velocity = shared::maxwell_boltzmann_velocity_distribution(sigma, dimensions(ctx));
			}
		}

		void apply(core::SimulationContext & ctx) const {
			if (target_temp == TemperatureNotSet) return;

			const vec3 avg_v = average_velocity(ctx);
			const double current_T = temperature(ctx, avg_v);
			const double diff = target_temp - current_T;
			const double new_T = current_T + std::clamp(diff, -max_temp_change, max_temp_change);

			if (std::abs(new_T - current_T) < 1e-12) return;

			const double factor = std::sqrt(new_T / current_T);
			scale_thermal_velocities(ctx, factor, avg_v);
		}

	private:
		double init_temp;
		double target_temp;
		double max_temp_change;

		static vec3 average_velocity(const core::SimulationContext& ctx) {
			vec3 sum{0, 0, 0};
			for (size_t i = ctx.index_start(); i < ctx.index_end(); ++i)
				sum += ctx.get_particle_by_index(i).velocity;
			return sum / static_cast<double>(ctx.size());
		}

		static double temperature(const core::SimulationContext& ctx, const vec3& avg_v) {
			double kinetic = 0.0;
			for (size_t i = ctx.index_start(); i < ctx.index_end(); ++i) {
				const auto p = ctx.get_particle_by_index(i);
				const vec3 dv = p.velocity - avg_v;
				kinetic += p.mass * dv.norm_squared();
			}
			const size_t dof = dimensions(ctx) * ctx.size();
			return kinetic / static_cast<double>(dof);
		}

		static void scale_thermal_velocities(core::SimulationContext& ctx, const double factor, const vec3& avg_v) {
			for (size_t i = ctx.index_start(); i < ctx.index_end(); ++i) {
				auto p = ctx.get_particle_by_index(i);
				p.velocity = avg_v + factor * (p.velocity - avg_v);
			}
		}

		static uint8_t dimensions(const core::SimulationContext& ctx) {
			return 3 -
				(ctx.box().extent.x == 0) -
				(ctx.box().extent.y == 0) -
				(ctx.box().extent.z == 0);
		}
	};

} // namespace md::env



