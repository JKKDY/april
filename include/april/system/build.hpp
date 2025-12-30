#pragma once

#include <vector>
#include "april/system/system.hpp"
#include "april/env/environment.hpp"
#include "april/containers/container.hpp"
#include "april/env/domain.hpp"
#include "april/system/build_helpers_particle.hpp"
#include "april/system/build_helpers_domain.hpp"
#include "april/system/build_helpers_boundary.hpp"


namespace april::core {

	struct BuildInfo {
		std::unordered_map<env::ParticleType, env::ParticleType> type_map;
		std::unordered_map<env::ParticleID, env::ParticleID> id_map;
		env::Domain particle_box;
		env::Domain simulation_domain;
	};


	template <class Container, env::IsEnvironment EnvT>
	requires container::IsContainerDecl<Container, typename EnvT::traits>
	System<Container, typename EnvT::traits> build_system(
		const EnvT & environment,
		const Container& container_config,
		BuildInfo * build_info
	) {
		using BoundaryTable = typename EnvT::traits::boundary_table_t;
		using ForceTable = typename EnvT::traits::force_table_t;
		using ParticleRecord = typename EnvT::traits::particle_record_t;

		using EnvData = env::internal::EnvironmentData< // explicit type so the IDE can perform code completion
			typename EnvT::traits::force_variant_t,
			typename EnvT::traits::boundary_variant_t,
			typename EnvT::traits::controller_storage_t,
			typename EnvT::traits::field_storage_t>;

		// get a copy of the environment data
		EnvData env = env::internal::get_env_data(environment);

		// validate & set simulation domain
		const env::Box particle_bbox = internal::particle_bounding_box(env.particles);
		const env::Box simulation_box = internal::determine_simulation_box(
			env.domain, particle_bbox, env.margin_abs, env.margin_fac);

		// validate & create Particles
		auto [type_pairs, id_pairs] = internal::extract_interaction_parameters(
			env.type_interactions, env.id_interactions );

		internal::assign_missing_particle_ids(env.particles, env.user_particle_ids);

		auto [type_map, id_map] = internal::create_particle_mappings(
			env.particles,
			env.user_particle_types,
			env.user_particle_ids,
			type_pairs,
			id_pairs
		);

		const std::vector<ParticleRecord> particles =
			internal::build_particles<typename EnvT::traits::user_data_t>(env.particles, type_map, id_map);

		// create boundary table
		internal::set_default_boundaries(env.boundaries);
		BoundaryTable boundaries(env.boundaries, simulation_box);
		auto topologies = internal::extract_topologies(boundaries);
		internal::validate_topologies(topologies);

		//  create force table
		ForceTable forces (env.type_interactions, env.id_interactions, type_map, id_map);

		// fill build info if given
		if (build_info) {
			build_info->type_map = type_map;
			build_info->id_map = id_map;
			build_info->particle_box = env::Domain(particle_bbox.min, particle_bbox.extent);
			build_info->simulation_domain = env::Domain(simulation_box.min, simulation_box.extent);
		}

		container::internal::ContainerCreateInfo container_info {
			.flags = internal::set_container_flags(topologies),
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
