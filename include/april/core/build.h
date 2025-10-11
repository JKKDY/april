#pragma once


#include <vector>
#include "april/core/system.h"
#include "april/env/particle.h"
#include "april/env/environment.h"
#include "april/containers/container.h"
#include "april/env/domain.h"


namespace april::core {
	namespace internal {
		env::Domain calculate_bounding_box(const std::vector<env::Particle>& particles);

		struct InteractionParams {
			bool pair_contains_types;
			std::pair<int,int> key_pair;
		};

		void validate_domain_params(
			const env::Domain& domain,
			const env::Domain& bbox
		);

		void validate_particle_params(
			const std::vector<env::Particle> & particles,
			std::vector<InteractionParams> interactions,
			const std::unordered_set<env::ParticleID>& usr_particle_ids,
			const std::unordered_set<env::ParticleType>& usr_particle_types
		);

		UserToInternalMappings map_ids_and_types_to_internal(
			std::vector<env::Particle>& particles,
			const std::vector<InteractionParams>& interactions,
			std::unordered_set<env::ParticleID>& usr_particle_ids,
			std::unordered_set<env::ParticleType>& usr_particle_types
		);

		env::Domain finalize_environment_domain(
			const env::Domain& bbox,
			const env::Domain& usr_domain
		);

		std::vector<env::internal::Particle> build_particles(
			const std::vector<env::Particle>& particle_infos,
			const UserToInternalMappings& mapping
		);
	}

	template <container::IsContDecl C, class FPack, class BPack>
	System<C, env::Environment<FPack, BPack>> build_system(
		env::Environment<FPack, BPack> & environment,
		const C& container,
		UserToInternalMappings* particle_mappings
	) {
		using EnvT = env::Environment<FPack, BPack>;
		using BoundaryTable = boundary::internal::BoundaryTable<typename EnvT::boundary_variant_t>;
		using namespace internal;

		auto & env = env::internal::get_env_data(environment);
		const env::Domain bbox = calculate_bounding_box(env.particles);

		std::vector<InteractionParams> interactions(env.interactions.size());
		for (size_t i = 0; i < interactions.size(); i++) {
			interactions[i].key_pair = env.interactions[i].key_pair;
			interactions[i].pair_contains_types = env.interactions[i].pair_contains_types;
		}

		validate_domain_params(env.domain, bbox);
		validate_particle_params(
			env.particles,
			interactions,
			env.usr_particle_ids,
			env.usr_particle_types);

		const UserToInternalMappings mapping = map_ids_and_types_to_internal(
			env.particles,
			interactions,
			env.usr_particle_ids,
			env.usr_particle_types
		);

		const env::Domain domain = finalize_environment_domain(bbox, env.domain);
		const std::vector<env::internal::Particle> particles = build_particles(env.particles, mapping);

		if (particle_mappings) {
			particle_mappings->usr_ids_to_impl_ids = mapping.usr_ids_to_impl_ids;
			particle_mappings->usr_types_to_impl_types = mapping.usr_types_to_impl_types;
		}

		BoundaryTable boundaries (env.boundaries, domain);

		return System<C, env::Environment<FPack, BPack>> (
			container,
			domain,
			particles,
			boundaries,
			mapping.usr_types_to_impl_types,
			mapping.usr_ids_to_impl_ids,
			env.interactions
		);
	}
}