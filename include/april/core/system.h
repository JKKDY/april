#pragma once

#include "april/env/particle.h"
#include "april/env/environment.h"
#include "april/env/interaction.h"
#include "april/algo/algorithm.h"

namespace april::core {

	template<algo::IsAlgoDecl Algo> class System {
		using Algorithm = typename Algo::impl;
		using Interactions = env::impl::InteractionManager;
	public:
		void update_forces() {
			algorithm.calculate_forces();
		}
		std::vector<const env::Particle*> particles(env::ParticleState state) {

		}
		[[nodiscard]] std::vector<const env::Particle*const> export_particles() const {

		}

	private:
		System(
			Algo algo,
			Interactions interactions,
			const env::Domain & domain,
			const std::vector<env::impl::Particle> & particles
		):
			interaction_manager(std::move(interactions)),
			algorithm(algo, interaction_manager, domain)
		{
			algorithm.build(particles);
		}


		env::impl::InteractionManager interaction_manager;
		Algorithm algorithm;

		double time{};


		friend System compile(env::Environment);
	};


	struct UserToInternalMappings  {
		std::unordered_map<env::ParticleType, env::impl::ParticleType> usr_types_to_impl_types;
		std::unordered_map<env::ParticleID, env::impl::ParticleID> usr_ids_to_impl_ids;
	};


	namespace impl {
		env::Domain calculate_bounding_box(const std::vector<env::Particle> & particles);

		void validate_domain_params(
			const env::Domain & domain,
			const env::Domain & bbox,
			const std::vector<env::Particle> & particles
		);

		void validate_particle_params(
			std::vector<env::impl::InteractionInfo>& interactions,
			const std::unordered_set<env::ParticleID> & usr_particle_ids,
			const std::unordered_set<env::ParticleType> & usr_particle_types
		);

		UserToInternalMappings map_ids_and_types_to_internal(
			 std::vector<env::Particle> & particles,
			 std::vector<env::impl::InteractionInfo> & interactions,
			 std::unordered_set<env::ParticleID> & usr_particle_ids,
			 std::unordered_set<env::ParticleType> & usr_particle_types
		 );

		env::Domain finalize_environment_domain(
			const env::Domain & bbox,
			const env::Domain & usr_domain
		);

		std::vector<env::impl::Particle> build_particles(
			const std::vector<env::Particle> & particle_infos,
			const UserToInternalMappings& mapping
		);
	}


	template<algo::IsAlgoDecl Algo>
	auto compile(
		const env::Environment & environment,
		const Algo & algorithm,
		UserToInternalMappings * particle_mappings = nullptr
		) -> System<typename Algo::impl>
	{
		using namespace impl;
		auto env = env::impl::get_env_data(environment);
		const env::Domain bbox = calculate_bounding_box(env.particles);

		validate_domain_params(env.domain, bbox, env.particles);
		validate_particle_params(env.interactions, env.usr_particle_ids, env.usr_particle_types);

		UserToInternalMappings mapping = map_ids_and_types_to_internal(env.particles,env.interactions, env.usr_particle_ids, env.usr_particle_types);
		env::Domain domain = finalize_environment_domain(bbox, env.domain);
		std::vector<env::impl::Particle> particles = build_particles(env.particles, mapping);

		env::impl::InteractionManager interaction_manager;
		interaction_manager.build(env.interactions, mapping.usr_types_to_impl_types, mapping.usr_ids_to_impl_ids);

		if (particle_mappings) {
			particle_mappings->usr_ids_to_impl_ids = mapping.usr_ids_to_impl_ids;
			particle_mappings->usr_types_to_impl_types = mapping.usr_types_to_impl_types;
		}

		System<Algo> system (algorithm, interaction_manager, domain);
		return system;
	}

}