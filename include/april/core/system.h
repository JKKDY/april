#pragma once

#include "april/env/particle.h"
#include "april/forces/force.h"
#include "april/env/environment.h"
#include "april/forces/force_table.h"
#include "april/containers/container.h"
#include "april/env/domain.h"
#include "april/boundaries/boundary_table.h"

namespace april::core {

	struct UserToInternalMappings {
		using TypeMap = std::unordered_map<env::ParticleType, env::internal::ParticleType>;
		using IdMap = std::unordered_map<env::ParticleID, env::internal::ParticleID>;
		TypeMap usr_types_to_impl_types;
		IdMap usr_ids_to_impl_ids;
	};

	template <container::IsContDecl C, class Env>
	class System;

	template <container::IsContDecl Cont, class FPack, class BPack>
	auto build_system(
		env::Environment<FPack, BPack>& environment,
		const Cont& container,
		UserToInternalMappings* particle_mappings = nullptr
	) -> System<Cont, env::Environment<FPack, BPack>>;


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


	template <container::IsContDecl C, force::IsForce ... Fs, boundary::IsBoundary ... BCs>
	class System <C, env::Environment<force::ForcePack<Fs...>, boundary::BoundaryPack<BCs...>>> {
	public:
		using EnvT          = env::Environment<force::ForcePack<Fs...>, boundary::BoundaryPack<BCs...>>;
		using Container     = typename C::template impl<EnvT>;
		using BoundaryTable = boundary::internal::BoundaryTable<typename EnvT::boundary_variant_t>;
		using ForceTable    = force::internal::ForceTable<EnvT>;
		using Interaction   = force::internal::InteractionInfo<typename EnvT::force_variant_t>;
		using Particle      = env::internal::Particle;
		using ParticleRef   = env::internal::ParticleRef;
		using ParticleView  = env::ParticleView;
		using ParticleID    = env::internal::ParticleID;

		const env::Domain domain;

	private:
		System(
			const C & container_cfg,
			const env::Domain& domain,
			const std::vector<Particle> & particles,
			const BoundaryTable boundaries,
			const UserToInternalMappings::TypeMap & usr_types_to_impl_types,
			const UserToInternalMappings::IdMap & usr_ids_to_impl_ids,
			std::vector<Interaction> & interaction_infos)
			: domain(domain), container(container_cfg), boundary_table(boundaries), time_(0)
		{
			force_table.build(interaction_infos, usr_types_to_impl_types, usr_ids_to_impl_ids);
			container.init(force_table, domain);
			container.dispatch_build(particles);
		}

		Container container;
		BoundaryTable boundary_table;
		ForceTable force_table;

		double time_;

		template <container::IsContDecl Cont, class FPack, class BPack>
		friend System<Cont, env::Environment<FPack, BPack>>
		build_system(
			env::Environment<FPack, BPack>& environment,
			 const Cont& container,
			 UserToInternalMappings* particle_mappings
		);

	public:
		void update_forces() {
			container.dispatch_calculate_forces();
		}

		void apply_boundary_conditions() {
			for (boundary::Face _ : boundary::faces) {

				// env::impl::CompiledBoundary<typename EnvT::boundary_variant_t> & boundary = boundary_table.get_boundary(face);
				// std::vector<size_t> particle_ids = container.dispatch_collect_indices_in_region(boundary.region);

			}
		}

		void apply_controllers() {

		}

		void apply_force_fields() {

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
			for (auto i = index_start(); i < index_end(); ++i) {
				auto & p = get_particle_by_index(i);
				if (static_cast<int>(p.state & state))
					particles.emplace_back(p);
			}
			return particles;
		}
	};
}
