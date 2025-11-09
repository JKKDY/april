#pragma once

#include <utility>
#include "april/common.h"
#include "april/env/domain.h"
#include "april/particle/descriptors.h"


namespace april::core::internal {
	// ---- Domain Validation & Setting ----

	// calculate the minimal bounding box that contains all particles
	inline env::Box particle_bounding_box(const std::vector<env::Particle>& particles) {
		if (particles.empty()) return {};

		vec3 min = particles[0].position;
		vec3 max = particles[0].position;

		for (const auto& p : particles) {
			min.x = std::min(min.x, p.position.x);
			min.y = std::min(min.y, p.position.y);
			min.z = std::min(min.z, p.position.z);

			max.x = std::max(max.x, p.position.x);
			max.y = std::max(max.y, p.position.y);
			max.z = std::max(max.z, p.position.z);
		}

		return {min, max};
	}


	// given particle bounding box and user set parameters, calculate the simulation box
	inline env::Box calculate_simulation_box(
		const env::Domain& desired_domain,
		const env::Box& required_box,
		const env::Box& particle_bbox
	) {
		// Case 1: fully manual: both user origin & extent are specified. overrides any margin set
		if (desired_domain.origin.has_value() && desired_domain.extent.has_value()) {
			return env::Box::from_domain(desired_domain);
		}
		// Case 2: fully automatic: both user origin & extent not set
		if (!desired_domain.origin.has_value() && !desired_domain.extent.has_value()) {
			return required_box;
		}
		// Case 3: user origin set, user extent not set
		if (desired_domain.origin.has_value() && !desired_domain.extent.has_value()) {
			// we know from validation that origin < bbox.min on all axis so we use that as min corner
			// max corner is chosen such that it satisfies margin requirements
			return {desired_domain.origin.value(), required_box.max};
		}
		// Case 4: user origin not set, user extent set
		if (!desired_domain.origin.has_value() && desired_domain.extent.has_value()) {
			// center the particle bounding box inside the simulation domain
			const vec3 bbox_center = (particle_bbox.min + particle_bbox.max) * 0.5;
			const vec3 origin = bbox_center - desired_domain.extent.value() / 2;
			return {origin, origin + desired_domain.extent.value()};
		}
		std::unreachable();
	}


	inline void verify_domain_consistency(const env::Box & simulation_box, const env::Box & particle_bbox) {
		if (simulation_box.extent.x < 0 || simulation_box.extent.y < 0 || simulation_box.extent.z < 0)
		{
			throw std::logic_error("Simulation domain has negative extent. Got extent " + simulation_box.extent.to_string());
		}

		if (simulation_box.extent.x == 0 && simulation_box.extent.y == 0 && simulation_box.extent.z == 0)
		{
			throw std::logic_error("Simulation domain size is zero. Got extent " + simulation_box.extent.to_string() +
				"\n If you have no particles or they all have the same position, you must specify a domain manually.");
		}

		// check that min corner is outside of particle bbox on all axis
		if (simulation_box.min.x > particle_bbox.min.x ||
			simulation_box.min.y > particle_bbox.min.y ||
			simulation_box.min.z > particle_bbox.min.z
		) {
			throw std::invalid_argument(
				"Specified Environment domain does not contain all particles: \n"
				"\tDomain box min corner: " + simulation_box.min.to_string() + "\n"
				"\tParticle bounding min corner: " + particle_bbox.min.to_string()
			);
		}

		// check that max corner is outside of particle bbox on all axis
		if (simulation_box.max.x < particle_bbox.max.x ||
			simulation_box.max.y < particle_bbox.max.y ||
			simulation_box.max.z < particle_bbox.max.z
		) {
			throw std::invalid_argument(
				"Specified Environment domain does not contain all particles: \n"
				"\tDomain box max corner: " + simulation_box.max.to_string() + "\n"
				"\tParticle bounding max corner: " + particle_bbox.max.to_string()
			);
		}
	}

	inline env::Box determine_simulation_box(
		const env::Domain& desired_domain,
		const env::Box& particle_bbox,
		const vec3 & margin_abs,
		const vec3 & margin_fac
	) {
		if (margin_abs.x < 0 || margin_abs.y < 0 || margin_abs.z < 0) {
			throw std::logic_error("Absolute margin was set to negative on at least one axis. Got: " + margin_abs.to_string());
		}

		if (margin_fac.x < 0 || margin_fac.y < 0 || margin_fac.z < 0) {
			throw std::logic_error("Margin factor was set to negative on at least one axis. Got: " + margin_fac.to_string());
		}

		const vec3 effective_margin = {
			std::max( particle_bbox.extent.x * margin_fac.x, margin_abs.x),
			std::max( particle_bbox.extent.y * margin_fac.y, margin_abs.y),
			std::max( particle_bbox.extent.z * margin_fac.z, margin_abs.z)
		};

		const env::Box required_box (particle_bbox.min-effective_margin, particle_bbox.max+effective_margin);

		const env::Box simulation_box = calculate_simulation_box(desired_domain, required_box, particle_bbox);
		verify_domain_consistency(simulation_box, particle_bbox);

		return {simulation_box.min, simulation_box.max};
	}

}