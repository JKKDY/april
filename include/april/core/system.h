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
		// private constructor since System should only be creatable through build_system(...)
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

		// System factory. Only valid way to create a System.
		template <container::IsContDecl Cont, class FPack, class BPack>
		friend System<Cont, env::Environment<FPack, BPack>>
		build_system(
			env::Environment<FPack, BPack>& environment,
			 const Cont& container,
			 UserToInternalMappings* particle_mappings
		);

	public:
		// call to update all pairwise forces between particles
		void update_forces() {
			container.dispatch_calculate_forces();
		}

		// call to register particle movements. This may cause container internals to change/be rebuilt
		void register_all_particle_movements() {
			container.dispatch_register_all_particle_movements();
		}

		void register_particle_movement(ParticleID id) {
			size_t idx = container.id_to_index(id);
			container.dispatch_register_particle_movement(idx);
		}

		// call to apply boundary conditions to all particles. should not be called before register_particle_movements
		void apply_boundary_conditions() {
			using Boundary = boundary::internal::CompiledBoundary<typename EnvT::boundary_variant_t>;

			auto box = env::Box(domain);

			for (boundary::Face face : boundary::all_faces) {

				Boundary & boundary = boundary_table.get_boundary(face);
				std::vector<size_t> particle_ids = container.dispatch_collect_indices_in_region(boundary.region);

				if (boundary.topology.boundary_thickness >= 0) {
					for (auto p_idx : particle_ids) {
						env::internal::Particle & p = container.dispatch_get_particle_by_index(p_idx);
						boundary.apply(p, box, face);

						if (boundary.topology.may_change_particle_position) {
							container.register_particle_movement(p_idx);
						}
					}
				} else {
					for (auto p_idx : particle_ids) {
						env::internal::Particle & p = container.dispatch_get_particle_by_index(p_idx);

						// make sure the particle exited through the current boundary face
						// solve for intersection of the particles path with the boundary face
						// with the equation y = t * diff + p where:
						// diff is the path traveled, p is the particles starting position and y is the face
						const int ax = axis_of_face(face);
						const vec3 diff = p.position - p.old_position;
						const double y = diff[ax] < 0 ? box.min[ax] : box.max[ax];
						const double t = (y - p.old_position[ax]) / diff[ax];

						const vec3 intersection = t * diff + p.old_position;

						// and check if that point is on the domains surface
						auto [ax1, ax2] = non_face_axis(face);
						if (box.max[ax1] >= intersection[ax1] && box.min[ax1] <= intersection[ax1] &&
							box.max[ax2] >= intersection[ax2] && box.min[ax2] <= intersection[ax2]) {
							boundary.apply(p, box, face);

							if (boundary.topology.may_change_particle_position) {
								container.register_particle_movement(p_idx);
							}
						}
					}
				}
			}
		}

		void apply_controllers() {

		}

		void apply_force_fields() {

		}

		// get a particle reference by its id. Usually slower than getting it by its index.
		// Useful for stable iterations and accessing a specific particle
		[[nodiscard]] Particle & get_particle_by_id(const ParticleID id) noexcept {
			return container.dispatch_get_particle_by_id(id);
		}

		// get the first particle id (usually 0)
		[[nodiscard]] ParticleID id_start() const noexcept{
			return container.dispatch_id_start();
		}

		// get the last particle id (usually n-1 with n = #particles)
		[[nodiscard]] ParticleID id_end() const noexcept {
			return container.dispatch_id_end();
		}

		// get a particle by its container specific id. useful for non-stable (but fast) iteration over particles
		[[nodiscard]] Particle & get_particle_by_index(const size_t index) noexcept {
			return container.dispatch_get_particle_by_index(index);
		}

		// get the first particle index (usually 0)
		[[nodiscard]] size_t index_start() const noexcept {
			return container.dispatch_index_start();
		}

		// get the last particle index (usually n-1 with n = #particles)
		[[nodiscard]] size_t index_end() const noexcept {
			return container.dispatch_index_end();
		}

		// returns the systems time
		[[nodiscard]] double time() const noexcept {
			return time_;
		}

		// propagate the systems time by a time step dt
		void update_time(const double dt) noexcept {
			time_ += dt;
		}

		// reset the systems time to 0
		void reset_time() noexcept {
			time_ = 0;
		}

		// get read access to all internal particles based on their state. Useful for snapshots and analysis.
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


