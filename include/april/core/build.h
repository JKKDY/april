#pragma once


#include <vector>
#include "april/core/system.h"
#include "april/env/particle.h"
#include "april/env/environment.h"
#include "april/containers/container.h"
#include "april/env/domain.h"


namespace april::core {

	namespace internal {
		struct InteractionParameters {
			bool pair_contains_types;
			std::pair<int,int> key_pair;
		};

		// swap origin and extent components such that for all axis i: origin[i] < extent[i]
		env::Domain clean_domain(
			const env::Domain & domain
		);

		// calculate the minimal bounding box that contains all particles
		env::Box calculate_bounding_box(
			const std::vector<env::Particle>& particles
		);

		// check if user set domain parameters are ok (e.g. if it contains all particles)
		void validate_domain_params(
			const env::Domain& domain,
			const env::Box& bbox
		);

		// check if particle parameters are ok (e.g. no duplicate ids, every particle type has a force specified)
		void validate_particle_params(
			const std::vector<env::Particle> & particles,
			std::vector<InteractionParameters> interactions,
			const std::unordered_set<env::ParticleID>& usr_particle_ids,
			const std::unordered_set<env::ParticleType>& usr_particle_types
		);

		// map user set particle ids & types to dense internal ids & types and return mappings
		UserToInternalMappings map_ids_and_types_to_internal(
			std::vector<env::Particle>& particles,
			const std::vector<InteractionParameters>& interactions,
			std::unordered_set<env::ParticleID>& usr_particle_ids,
			std::unordered_set<env::ParticleType>& usr_particle_types
		);

		// given particle bounding box and user set parameters, calculate the final simulation domain
		env::Box finalize_environment_domain(
			const env::Box& bbox,
			const env::Domain& usr_domain,
			const vec3 & margin_abs,
			const vec3 & margin_fac
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
	}

	template <container::IsContDecl Container, env::IsEnvironment Environment>
	auto build_system(
		const Environment & environment,
		const Container& container,
		UserToInternalMappings* particle_mappings
	) {
		using BoundaryTable = typename Environment::traits::boundary_table_t;
		using namespace internal;

		auto env = env::internal::get_env_data(environment);

		// --- validate & set simulation domain ---
		const env::Box bbox = calculate_bounding_box(env.particles);
		const env::Domain domain_cleaned = clean_domain(env.domain);

		validate_domain_params(domain_cleaned, bbox);
		const env::Box box = finalize_environment_domain(bbox, domain_cleaned, env.margin_abs, env.margin_fac);
		const env::Domain domain = {box.min, box.extent};


		// --- validate & create Particles ---
		std::vector<InteractionParameters> interactions(env.interactions.size());
		for (size_t i = 0; i < interactions.size(); i++) {
			interactions[i].key_pair = env.interactions[i].key_pair;
			interactions[i].pair_contains_types = env.interactions[i].pair_contains_types;
		}

		validate_particle_params(
			env.particles,
			interactions,
			env.user_particle_ids,
			env.user_particle_types);

		const UserToInternalMappings mapping = map_ids_and_types_to_internal(
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

		return System<Container, typename Environment::traits> (
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