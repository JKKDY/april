#pragma once

#include "april/env/particle.h"
#include "april/env/environment.h"
#include "april/env/interaction.h"
#include "april/algo/algorithm.h"

namespace april::core {


	class System {
		// using Algorithm = typename Algo::impl;
		using Interactions = env::impl::InteractionManager;
		using Particle = env::impl::Particle;
		using ParticleRef = env::impl::ParticleRef;
		using ParticleView = env::impl::ParticleView;
		using ParticleID = env::impl::ParticleID;

	public:
		void update_forces() const {
			algorithm->calculate_forces();
		}

		[[nodiscard]] Particle & get_particle_by_id(const ParticleID id) const noexcept {
			return algorithm->get_particle_by_id(id);
		}

		[[nodiscard]] ParticleID id_start() const noexcept{
			return algorithm->id_start();
		}

		[[nodiscard]] ParticleID id_end() const noexcept {
			return algorithm->id_end();
		}

		[[nodiscard]] Particle & get_particle_by_index(const size_t index) const noexcept {
			return algorithm->get_particle_by_index(index);
		}

		[[nodiscard]] size_t index_start() const noexcept {
			return algorithm->index_start();
		}

		[[nodiscard]] size_t index_end() const noexcept {
			return algorithm->index_end();
		}

		[[nodiscard]] double t() const noexcept {
			return time;
		}

		[[nodiscard]] std::vector<ParticleView> export_particles() const {
			std::vector<ParticleView> particles;
			particles.reserve(index_end());
			for (auto i = index_start(); i < index_end(); i++) {
				particles.emplace_back(get_particle_by_index(i));
			}
			return {};
		}

	private:
		System(
			std::unique_ptr<algo::impl::IAlgorithm> algorithm, Interactions interactions,
			const env::Domain& domain, const std::vector<env::impl::Particle>& particles)
			: algorithm(std::move(algorithm)),
			  interaction_manager(std::move(interactions)) {
			algorithm->init(interactions, domain);
			algorithm->build(particles);
		}

		std::unique_ptr<algo::impl::IAlgorithm> algorithm;
		env::impl::InteractionManager interaction_manager;

		double time{};


		friend System compile(env::Environment);
	};


	struct UserToInternalMappings {
		std::unordered_map<env::ParticleType, env::impl::ParticleType> usr_types_to_impl_types;
		std::unordered_map<env::ParticleID, env::impl::ParticleID> usr_ids_to_impl_ids;
	};


	namespace impl {
		env::Domain calculate_bounding_box(const std::vector<env::Particle>& particles);

		void validate_domain_params(
			const env::Domain& domain,
			const env::Domain& bbox
		);

		void validate_particle_params(
			std::vector<env::impl::InteractionInfo>& interactions,
			const std::unordered_set<env::ParticleID>& usr_particle_ids,
			const std::unordered_set<env::ParticleType>& usr_particle_types
		);

		UserToInternalMappings map_ids_and_types_to_internal(
			std::vector<env::Particle>& particles,
			std::vector<env::impl::InteractionInfo>& interactions,
			std::unordered_set<env::ParticleID>& usr_particle_ids,
			std::unordered_set<env::ParticleType>& usr_particle_types
		);

		env::Domain finalize_environment_domain(
			const env::Domain& bbox,
			const env::Domain& usr_domain
		);

		std::vector<env::impl::Particle> build_particles(
			const std::vector<env::Particle>& particle_infos,
			const UserToInternalMappings& mapping
		);
	}


	template <algo::impl::IsAlgoDecl Algo>
	auto compile(
		env::Environment& environment,
		const Algo& algorithm,
		UserToInternalMappings* particle_mappings = nullptr
	) -> System {
		using namespace impl;
		auto & env = env::impl::get_env_data(environment);
		const env::Domain bbox = calculate_bounding_box(env.particles);

		validate_domain_params(env.domain, bbox);
		validate_particle_params(env.interactions, env.usr_particle_ids, env.usr_particle_types);

		UserToInternalMappings mapping = map_ids_and_types_to_internal(
			env.particles,
			env.interactions,
			env.usr_particle_ids,
			env.usr_particle_types
		);

		env::Domain domain = finalize_environment_domain(bbox, env.domain);
		std::vector<env::impl::Particle> particles = build_particles(env.particles, mapping);

		env::impl::InteractionManager interaction_manager;
		interaction_manager.build(env.interactions, mapping.usr_types_to_impl_types, mapping.usr_ids_to_impl_ids);

		auto algo = std::make_unique<Algo::impl>(algorithm);

		if (particle_mappings) {
			particle_mappings->usr_ids_to_impl_ids = mapping.usr_ids_to_impl_ids;
			particle_mappings->usr_types_to_impl_types = mapping.usr_types_to_impl_types;
		}

		System system(std::move(algo), algorithm, interaction_manager, domain);
		return system;
	}
}
