#pragma once

#include "april/env/particle.h"
#include "april/env/environment.h"
#include "april/env/interaction.h"
#include "april/algo/algorithm.h"

namespace april::core {

	struct UserToInternalMappings {
		std::unordered_map<env::ParticleType, env::impl::ParticleType> usr_types_to_impl_types;
		std::unordered_map<env::ParticleID, env::impl::ParticleID> usr_ids_to_impl_ids;
	};

	template <algo::impl::IsAlgoDecl Algo> class System;

	template <algo::impl::IsAlgoDecl Algo> System<Algo> compile(
		env::Environment& environment,
		const Algo& algorithm,
		UserToInternalMappings* particle_mappings = nullptr
	);


	template<class S>
	concept IsSystem = requires(S s, size_t i, typename S::ParticleID pid) {
		// force update
		{ s.update_forces() } -> std::same_as<void>;

		// particle access
		{ s.get_particle_by_id(pid) } -> std::same_as<typename S::Particle&>;
		{ s.get_particle_by_index(i) } -> std::same_as<typename S::Particle&>;

		{ s.id_start() } -> std::same_as<typename S::ParticleID>;
		{ s.id_end()   } -> std::same_as<typename S::ParticleID>;

		{ s.index_start() } -> std::same_as<size_t>;
		{ s.index_end()   } -> std::same_as<size_t>;

		// time query
		{ s.t() } -> std::convertible_to<double>;

		// export
		{ s.export_particles() } -> std::same_as<std::vector<typename S::ParticleView>>;
	};

	template <algo::impl::IsAlgoDecl Algo> class System {
	public:
		using Algorithm = typename Algo::impl;
		using Interactions = env::impl::InteractionManager;
		using Particle = env::impl::Particle;
		using ParticleRef = env::impl::ParticleRef;
		using ParticleView = env::impl::ParticleView;
		using ParticleID = env::impl::ParticleID;

		void update_forces() {
			algorithm.dispatch_calculate_forces();
		}

		[[nodiscard]] Particle & get_particle_by_id(const ParticleID id) noexcept {
			return algorithm.dispatch_get_particle_by_id(id);
		}

		[[nodiscard]] ParticleID id_start() const noexcept{
			return algorithm.dispatch_id_start();
		}

		// inclusive bound
		[[nodiscard]] ParticleID id_end() const noexcept {
			return algorithm.dispatch_id_end();
		}

		[[nodiscard]] Particle & get_particle_by_index(const size_t index) noexcept {
			return algorithm.dispatch_get_particle_by_index(index);
		}

		[[nodiscard]] size_t index_start() const noexcept {
			return algorithm.dispatch_index_start();
		}

		// inclusive bound
		[[nodiscard]] size_t index_end() const noexcept {
			return algorithm.dispatch_index_end();
		}

		[[nodiscard]] double t() const noexcept {
			return time;
		}

		[[nodiscard]] std::vector<ParticleView> export_particles(const env::ParticleState state = env::ParticleState::ALL) {
			std::vector<ParticleView> particles;
			if (algorithm.dispatch_particle_count() == 0) {
				return {};
			}
			particles.reserve(index_end() - index_start() + 1);
			for (auto i = index_start(); i <= index_end(); ++i) {
				auto & p = get_particle_by_index(i);
				if (static_cast<int>(p.state & state))
					particles.emplace_back(p);
			}
			return particles;
		}

		const env::Domain domain;

	private:
		System(
			const Algo & algo_cfg,
			const env::Domain& domain, const std::vector<env::impl::Particle>& particles,
			std::vector<env::impl::InteractionInfo> & interaction_infos,
			const std::unordered_map<env::ParticleType, env::impl::ParticleType> & usr_types_to_impl_types,
			const std::unordered_map<env::ParticleID, env::impl::ParticleID> & usr_ids_to_impl_ids)
			: domain(domain), algorithm(algo_cfg)
		{
			interaction_manager.build(interaction_infos, usr_types_to_impl_types, usr_ids_to_impl_ids);
			algorithm.init(interaction_manager, domain);
			algorithm.dispatch_build(particles);
		}

		Algorithm algorithm;
		env::impl::InteractionManager interaction_manager;

		double time{};

		template <algo::impl::IsAlgoDecl AlgoCfg> friend System<AlgoCfg> compile(
			env::Environment& environment,
			const AlgoCfg& algorithm,
			UserToInternalMappings* particle_mappings
		);
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


	template <algo::impl::IsAlgoDecl AlgoCfg> System<AlgoCfg> compile(
		env::Environment& environment,
		const AlgoCfg& algorithm,
		UserToInternalMappings* particle_mappings
	) {
		using namespace impl;
		auto & env = env::impl::get_env_data(environment);
		const env::Domain bbox = calculate_bounding_box(env.particles);

		validate_domain_params(env.domain, bbox);
		validate_particle_params(env.interactions, env.usr_particle_ids, env.usr_particle_types);

		const UserToInternalMappings mapping = map_ids_and_types_to_internal(
			env.particles,
			env.interactions,
			env.usr_particle_ids,
			env.usr_particle_types
		);

		const env::Domain domain = finalize_environment_domain(bbox, env.domain);
		const std::vector<env::impl::Particle> particles = build_particles(env.particles, mapping);

		if (particle_mappings) {
			particle_mappings->usr_ids_to_impl_ids = mapping.usr_ids_to_impl_ids;
			particle_mappings->usr_types_to_impl_types = mapping.usr_types_to_impl_types;
		}

		return System (algorithm, domain, particles,
			env.interactions, mapping.usr_types_to_impl_types, mapping.usr_ids_to_impl_ids);
	}
}
