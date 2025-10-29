#pragma once


#include <vector>
#include "april/core/system.h"
#include "april/env/particle.h"
#include "april/env/environment.h"
#include "april/containers/container.h"
#include "april/env/domain.h"


namespace april::core {



	namespace internal {

		// calculate the minimal bounding box that contains all particles
		env::Box particle_bounding_box(const std::vector<env::Particle>& particles);

		// given particle bounding box and user set parameters, calculate the simulation box
		env::Box determine_simulation_box(
			const env::Domain& desired_domain,
			const env::Box& particle_bbox,
			const vec3 & margin_abs,
			const vec3 & margin_fac
		);

		// give particles without and id a new id
		void assign_missing_particle_ids(
			std::vector<env::Particle>& particles,
			std::unordered_set<env::ParticleID>& user_ids
		);

		// map user set particle ids & types to dense internal ids & types and return mappings
		std::pair<std::unordered_map<env::ParticleType, env::internal::ParticleType>,
			std::unordered_map<env::ParticleID, env::internal::ParticleID>>
		create_particle_mappings(
			const std::vector<env::Particle>& particles,
			const std::unordered_set<env::ParticleType>& user_types,
			const std::unordered_set<env::ParticleID>& user_ids,
			const std::vector<std::pair<env::ParticleType, env::ParticleType>>& type_pairs,
			const std::vector<std::pair<env::ParticleID, env::ParticleID>>& id_pairs
		);

		// build internal particle representation from user data
		std::vector<env::internal::Particle> build_particles(
			const std::vector<env::Particle>& particle_infos,
			const std::unordered_map<env::ParticleType, env::internal::ParticleType> & type_map,
			const std::unordered_map<env::ParticleID, env::internal::ParticleID> & id_map
		);

		// get container flags from boundary topologies
		container::internal::ContainerFlags set_container_flags(
			const std::vector<boundary::Topology> & topologies
		);

		template<force::internal::IsForceVariant FV>
		auto extract_interaction_parameters(
			const std::vector<force::internal::TypeInteraction<FV>> & type_interactions,
			const std::vector<force::internal::IdInteraction<FV>> & id_interaction)
		{
			std::vector<std::pair<env::ParticleType, env::ParticleType>> type_pairs(type_interactions.size());
			std::vector<std::pair<env::ParticleID, env::ParticleID>> id_pairs(id_interaction.size());

			for (size_t i = 0; i < type_interactions.size(); i++)
				type_pairs[i] = {type_interactions[i].type1, type_interactions[i].type2};

			for (size_t i = 0; i < id_interaction.size(); i++)
				id_pairs[i] = {id_interaction[i].id1, id_interaction[i].id2};

			return std::pair {type_pairs, id_pairs};
		}

		template<boundary::internal::IsBoundaryVariant BV>
		auto set_default_boundaries(std::array<BV, 6>  & boundaries) {
			for (auto & v : boundaries)
				if (std::holds_alternative<boundary::internal::BoundarySentinel>(v))
					v.template emplace<boundary::Open>(); // default-construct Open boundary
		}

		template<class BoundaryTable>
		auto extract_topologies(const BoundaryTable  & boundaries) {
			std::vector<boundary::Topology> topologies;
			for (boundary::Face face : boundary::all_faces) {
				topologies.push_back(boundaries.get_boundary(face).topology);
			}
			return topologies;
		}

		void validate_topologies(const std::vector<boundary::Topology> & topologies);

	}



	template <container::IsContDecl Container, env::IsEnvironment EnvT>
	auto build_system(
		const EnvT & environment,
		const Container& container,
		BuildInfo * build_info
	) {
		using BoundaryTable = typename EnvT::traits::boundary_table_t;
		using ForceTable = typename EnvT::traits::force_table_t;

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
		const std::vector<env::internal::Particle> particles = internal::build_particles(env.particles, type_map, id_map);

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

		return System<Container, typename EnvT::traits> (
			container,
			internal::set_container_flags(topologies),
			simulation_box,
			particles,
			boundaries,
			forces,
			env.controllers,
			env.fields
		);
	}
}