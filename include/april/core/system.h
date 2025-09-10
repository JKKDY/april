#pragma once

#include "april/env/particle.h"
#include "april/env/force.h"
#include "april/env/environment.h"
#include "april/env/interaction.h"
#include "april/containers/container.h"

namespace april::core {

	struct UserToInternalMappings {
		std::unordered_map<env::ParticleType, env::impl::ParticleType> usr_types_to_impl_types;
		std::unordered_map<env::ParticleID, env::impl::ParticleID> usr_ids_to_impl_ids;
	};

	template <cont::impl::IsContDecl C, class Env>
	class System;

	template <cont::impl::IsContDecl C, class FPack>
	System<C, env::Environment<FPack>> build_system(
		env::Environment<FPack>& environment,
		const C& container,
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
		{ s.time() } -> std::convertible_to<double>;

		// export
		{ s.export_particles() } -> std::same_as<std::vector<typename S::ParticleView>>;
	};

	template <cont::impl::IsContDecl C, env::IsForce ... Fs>
	class System <C, env::Environment<env::ForcePack<Fs...>>> {
	public:
		using EnvT         = env::Environment<env::ForcePack<Fs...>>;
		using Container    = typename C::template impl<env::Environment<env::ForcePack<Fs...>>>;

		using Particle     = env::impl::Particle;
		using ParticleRef  = env::impl::ParticleRef;
		using ParticleView = env::impl::ParticleView;
		using ParticleID   = env::impl::ParticleID;


		void update_forces() {
			container.dispatch_calculate_forces();
		}

		[[nodiscard]] Particle & get_particle_by_id(const ParticleID id) noexcept {
			return container.dispatch_get_particle_by_id(id);
		}

		[[nodiscard]] ParticleID id_start() const noexcept{
			return container.dispatch_id_start();
		}

		// inclusive bound
		[[nodiscard]] ParticleID id_end() const noexcept {
			return container.dispatch_id_end();
		}

		[[nodiscard]] Particle & get_particle_by_index(const size_t index) noexcept {
			return container.dispatch_get_particle_by_index(index);
		}

		[[nodiscard]] size_t index_start() const noexcept {
			return container.dispatch_index_start();
		}

		// inclusive bound
		[[nodiscard]] size_t index_end() const noexcept {
			return container.dispatch_index_end();
		}

		[[nodiscard]] double time() const noexcept {
			return time_;
		}

		void update_time(const double dt) noexcept {
			time_ += dt;
		}

		void reset_time() noexcept {
			time_ = 0;
		}

		[[nodiscard]] std::vector<ParticleView> export_particles(const env::ParticleState state = env::ParticleState::ALL) {
			std::vector<ParticleView> particles;
			if (container.dispatch_particle_count() == 0) {
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
			const C & algo_cfg,
			const env::Domain& domain, const std::vector<env::impl::Particle>& particles,
			std::vector<env::impl::InteractionInfo<typename EnvT::force_variant_t>> & interaction_infos,
			const std::unordered_map<env::ParticleType, env::impl::ParticleType> & usr_types_to_impl_types,
			const std::unordered_map<env::ParticleID, env::impl::ParticleID> & usr_ids_to_impl_ids)
			: domain(domain), container(algo_cfg), time_(0)
		{
			interaction_manager.build(interaction_infos, usr_types_to_impl_types, usr_ids_to_impl_ids);
			container.init(interaction_manager, domain);
			container.dispatch_build(particles);
		}

		Container container;
		env::impl::InteractionManager<EnvT> interaction_manager;

		// TODO make System time aware
		double time_;

		template <cont::impl::IsContDecl Cont, class FPack>
		friend System<Cont, env::Environment<FPack>>
		build_system(
			env::Environment<FPack>& environment,
			 const Cont& container,
			 UserToInternalMappings* particle_mappings
			 );
	};


	namespace impl {
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
			std::vector<InteractionParams> interactions,
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

	template <cont::impl::IsContDecl C, class FPack>
	System<C, env::Environment<FPack>> build_system(
		env::Environment<FPack>& environment,
		const C& container,
		UserToInternalMappings* particle_mappings
	) {
		using namespace impl;
		auto & env = env::impl::get_env_data(environment);
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
		const std::vector<env::impl::Particle> particles = build_particles(env.particles, mapping);

		if (particle_mappings) {
			particle_mappings->usr_ids_to_impl_ids = mapping.usr_ids_to_impl_ids;
			particle_mappings->usr_types_to_impl_types = mapping.usr_types_to_impl_types;
		}

		return System<C, env::Environment<FPack>> (container, domain, particles,
			env.interactions, mapping.usr_types_to_impl_types, mapping.usr_ids_to_impl_ids);
	}
}
