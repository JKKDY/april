#pragma once


#include <vector>
#include "april/core/system.h"
#include "april/env/particle.h"
#include "april/env/environment.h"
#include "april/containers/container.h"
#include "april/env/domain.h"


namespace april::core {

	namespace internal {
		// swap origin and extent components such that for all axis i: origin[i] < extent[i]
		// env::Box normalize_domain( const env::Domain & domain);

		// calculate the minimal bounding box that contains all particles
		env::Box particle_bounding_box(const std::vector<env::Particle>& particles);

		// check if user set domain parameters are ok (e.g. if it contains all particles)
		// check that all particles are contained in the box specified by extent and origin

		// void verify_domain_consistency(const env::Box& domain_box, const env::Box& particle_bbox);

		// given particle bounding box and user set parameters, calculate the simulation box
		env::Box determine_simulation_box(
			const env::Domain& desired_domain,
			const env::Box& particle_bbox,
			const vec3 & margin_abs,
			const vec3 & margin_fac
		);


		// check if particle parameters are ok (e.g. no duplicate ids, every particle type has a force specified)
		// void validate_particle_data(
		// 	const std::vector<env::Particle>& particles,
		// 	const std::unordered_set<env::ParticleType>& user_types,
		// 	const std::unordered_set<env::ParticleID>& user_ids,
		// 	const std::vector<std::pair<env::ParticleType, env::ParticleType>>& type_pairs,
		// 	const std::vector<std::pair<env::ParticleID, env::ParticleID>>& id_pairs
		// );

		// map user set particle ids & types to dense internal ids & types and return mappings
		UserToInternalMappings create_particle_mappings(
			std::vector<env::Particle>& particles,
			const std::vector<InteractionParameters>& interactions,
			std::unordered_set<env::ParticleID>& usr_particle_ids,
			std::unordered_set<env::ParticleType>& usr_particle_types
		);


		// build internal particle representation from user data
		std::vector<env::internal::Particle> build_particles(
			const std::vector<env::Particle>& particle_infos,
			const UserToInternalMappings& mapping
		);

		// set container flags
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

			for (size_t i = 0; i < type_interactions.size(); i++) {
				type_pairs[i] = {type_interactions[i].type1, type_interactions[i].type2};
			}

			for (size_t i = 0; i < id_interaction.size(); i++) {
				id_pairs[i] = {id_interaction[i].id1, id_interaction[i].id2};
			}

			return std::pair {type_pairs, id_pairs};
		}

	}

	template <container::IsContDecl Container, env::IsEnvironment EnvT>
	auto build_system(
		const EnvT & environment,
		const Container& container,
		UserToInternalMappings* particle_mappings
	) {
		using BoundaryTable = typename EnvT::traits::boundary_table_t;
		using EnvData = env::internal::EnvironmentData< // explicit type so the IDE can perform code completion
			typename EnvT::traits::force_variant_t,
			typename EnvT::traits::boundary_variant_t,
			typename EnvT::traits::controller_storage_t,
			typename EnvT::traits::field_storage_t>;
		using namespace internal;

		EnvData env = env::internal::get_env_data(environment);

		// validate & set simulation domain
		const env::Box particle_bbox = particle_bounding_box(env.particles);
		const env::Box simulation_box = determine_simulation_box(env.domain, particle_bbox, env.margin_abs, env.margin_fac);




		// --- validate & create Particles ---
		auto [type_pairs, id_pairs] = internal::extract_interaction_parameters(
			env.type_interactions, env.id_interactions
		);


		validate_particle_params(
			env.particles,
			interactions,
			env.user_particle_ids,
			env.user_particle_types);

		const UserToInternalMappings mapping = create_particle_mappings(
			env.particles,
			interactions,
			env.user_particle_ids,
			env.user_particle_types
		);

		const std::vector<env::internal::Particle> particles = build_particles(env.particles, mapping);


		// --- create boundary table ---
		for (auto & v : env.boundaries)
			if (std::holds_alternative<boundary::internal::BoundarySentinel>(v))
				v.template emplace<boundary::Open>(); // default-construct Open boundary

		BoundaryTable boundaries(env.boundaries, domain);

		std::vector<boundary::Topology> topologies;
		for (boundary::Face face : boundary::all_faces) {
			topologies.push_back(boundaries.get_boundary(face).topology);
		}

		// TODO validate topologies e.g face linkage


		// return mappings
		if (particle_mappings) {
			particle_mappings->user_ids_to_impl_ids = mapping.user_ids_to_impl_ids;
			particle_mappings->user_types_to_impl_types = mapping.user_types_to_impl_types;
		}

		return System<Container, typename EnvT::traits> (
			container,
			set_container_flags(topologies),
			domain,
			particles,
			boundaries,
			env.controllers,
			env.fields,
			mapping.user_types_to_impl_types,
			mapping.user_ids_to_impl_ids,
			env.interactions
		);
	}
}