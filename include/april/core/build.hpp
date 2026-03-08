#pragma once

#include <vector>

#include "april/core/system.hpp"
#include "april/core/environment.hpp"
#include "april/core/domain.hpp"
#include "april/core/internal/build_helpers_particle.hpp"
#include "april/core/internal/build_helpers_domain.hpp"
#include "april/core/internal/build_helpers_boundary.hpp"
#include "april/containers/container.hpp"


namespace april {

	struct BuildInfo {
		std::unordered_map<ParticleType, ParticleType> type_map;
		std::unordered_map<ParticleID, ParticleID> id_map;
		Domain particle_box;
		Domain simulation_domain;
	};


	template <class Container, core::IsEnvironment EnvT>
	requires container::IsContainerDecl<Container, typename EnvT::traits>
	System<Container, typename EnvT::traits> build_system(
		const EnvT & environment,
		const Container& container_config,
		BuildInfo * build_info
	) {
		using BoundaryTable = EnvT::traits::boundary_table_t;
		using ForceTable = EnvT::traits::force_table_t;
		using ParticleRecord = EnvT::traits::particle_record_t;

		using EnvData = core::internal::EnvironmentData< // explicit type so the IDE can perform code completion
			typename EnvT::traits::force_variant_t,
			typename EnvT::traits::boundary_variant_t,
			typename EnvT::traits::controller_storage_t,
			typename EnvT::traits::field_storage_t>;

		// get a copy of the environment data
		EnvData env = core::internal::get_env_data(environment);

		// validate & set simulation domain
		const core::Box particle_bbox = core::internal::particle_bounding_box(env.particles);
		const core::Box simulation_box = core::internal::determine_simulation_box(
			env.domain, particle_bbox, env.margin_abs, env.margin_fac);

		// validate & create Particles
		auto [type_pairs, id_pairs] = core::internal::extract_interaction_parameters(
			env.type_interactions, env.id_interactions );

		core::internal::assign_missing_particle_ids(env.particles, env.user_particle_ids);

		auto [type_map, id_map] = core::internal::create_particle_mappings(
			env.particles,
			env.user_particle_types,
			env.user_particle_ids,
			type_pairs,
			id_pairs
		);

		const std::vector<ParticleRecord> particles =
			core::internal::build_particles<typename EnvT::traits::particle_attributes_t>(env.particles, type_map, id_map);

		// create boundary table
		core::internal::set_default_boundaries(env.boundaries);
		BoundaryTable boundaries(env.boundaries, simulation_box);
		auto topologies = core::internal::extract_topologies(boundaries);
		core::internal::validate_topologies(topologies);

		//  create force table
		ForceTable forces (env.type_interactions, env.id_interactions, type_map, id_map);

		// fill build info if given
		if (build_info) {
			build_info->type_map = type_map;
			build_info->id_map = id_map;
			build_info->particle_box = Domain(particle_bbox.min, particle_bbox.extent);
			build_info->simulation_domain = Domain(simulation_box.min, simulation_box.extent);
		}

		container::internal::ContainerCreateInfo container_info {
			.flags = core::internal::set_container_flags(topologies),
			.hints = container::internal::ContainerHints(),
			.force_schema = forces.generate_schema(),
			.domain = simulation_box
		};

		return System<Container, typename EnvT::traits> (
			container_config,
			container_info,
			simulation_box,
			particles,
			boundaries,
			forces,
			env.controllers,
			env.fields
		);
	}
}















