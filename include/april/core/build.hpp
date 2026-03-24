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


	// constructs a system from an environment and a container.
	template <class ContainerCfg, core::IsEnvironment Env, exec::IsExecutionConfig ExecCfg>
	requires container::IsContainerDecl<ContainerCfg, typename Env::traits>
	auto build_system(
		const Env & environment,
		const ContainerCfg & container_config,
		const ExecCfg & execution_config,
		BuildInfo * build_info
	) {
		using BoundaryTable = Env::traits::boundary_table_t;
		using ForceTable = Env::traits::force_table_t;
		using ParticleRecord = Env::traits::particle_record_t;
		using ParticleAttributes = Env::traits::particle_attributes_t;

		using EnvData = core::internal::EnvironmentData< // explicit type so the IDE can perform code completion
			typename Env::traits::force_variant_t,
			typename Env::traits::boundary_variant_t,
			typename Env::traits::controller_storage_t,
			typename Env::traits::field_storage_t>;

		// get a copy of the environment data
		EnvData env = core::internal::get_env_data(environment);

		// validate & set simulation domain
		const core::Box particle_bbox = core::internal::particle_bounding_box(env.particles);
		const core::Box simulation_box = core::internal::determine_simulation_box(
			env.domain, particle_bbox, env.margin_abs, env.margin_fac);

		// validate & create Particles
		auto [type_pairs, id_pairs] =
			core::internal::extract_interaction_parameters(env.type_interactions, env.id_interactions );

		core::internal::assign_missing_particle_ids(env.particles, env.user_particle_ids);

		auto [type_map, id_map] = core::internal::create_particle_mappings(
			env.particles,
			env.user_particle_types,
			env.user_particle_ids,
			type_pairs,
			id_pairs
		);

		const std::vector<ParticleRecord> particles =
			core::internal::build_particles<typename Env::traits::particle_attributes_t>(env.particles, type_map, id_map);

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

		container::ContainerBuildConfig<ContainerCfg, ExecCfg, ParticleAttributes> container_info {
			.execution_config = execution_config,
			.config = container_config,
			.flags = core::internal::set_container_flags(topologies),
			.hints = container::ContainerHints(),
			.interaction_map = forces.generate_interaction_map(),
			.domain = simulation_box
		};

		return System<decltype(container_info), typename Env::traits, ExecCfg> (
			execution_config,
			container_info,
			particles,
			boundaries,
			forces,
			env.controllers,
			env.fields
		);
	}

	// convenience overloads
	template <class Container, core::IsEnvironment EnvT>
	auto build_system(
		const EnvT & environment,
		const Container & container_config
	) {
		return build_system(environment, container_config, ExecutionConfig(), nullptr);
	}

	template <class Container, core::IsEnvironment EnvT>
	auto build_system(
		const EnvT & environment,
		const Container & container_config,
		BuildInfo * build_info
	) {
		return build_system(environment, container_config, ExecutionConfig(), build_info);
	}

	template <class Container, core::IsEnvironment EnvT>
	auto build_system(
		const EnvT & environment,
		const Container & container_config,
		const exec::IsExecutionConfig auto & execution_config
	) {
		return build_system(environment, container_config, execution_config, nullptr);
	}
}


