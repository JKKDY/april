#pragma once

#include <memory>
#include "april/env/particle.h"
#include "april/forces/force.h"
#include "april/env/environment.h"
#include "april/forces/force_table.h"
#include "april/containers/container.h"
#include "april/env/domain.h"
#include "april/boundaries/boundary_table.h"
#include "april/core/context.h"
#include "april/shared/pack_storage.h"

namespace april::core {

	struct UserToInternalMappings {
		using TypeMap = std::unordered_map<env::ParticleType, env::internal::ParticleType>;
		using IdMap = std::unordered_map<env::ParticleID, env::internal::ParticleID>;
		TypeMap usr_types_to_impl_types;
		IdMap usr_ids_to_impl_ids;
	};

	template <container::IsContDecl C, class Env>
	class System;

	template <container::IsContDecl Cont, class FPack, class BPack, class CPack, class FFPack>
	auto build_system(
		env::Environment<FPack, BPack, CPack, FFPack>& environment,
		const Cont& container,
		UserToInternalMappings* particle_mappings = nullptr
	) -> System<Cont, env::Environment<FPack, BPack, CPack, FFPack>>;


	// Default: assume any type is not a System
	template<typename>
	inline constexpr bool is_system_v = false;

	// Specialization: mark all System<C, Env> instantiations as true
	template<container::IsContDecl C, class Env>
	inline constexpr bool is_system_v<System<C, Env>> = true;

	// Concept: true if T (after removing cv/ref) is a System specialization
	template<typename T>
	concept IsSystem = is_system_v<std::remove_cvref_t<T>>;


	template <
		container::IsContDecl C,
		force::IsForce ... Fs,
		boundary::IsBoundary ... BCs,
		controller::IsController ... Cs,
		field::IsField ... FFs
	>
	class System <C,
		env::Environment<
			force::ForcePack<Fs...>,
			boundary::BoundaryPack<BCs...>,
			controller::ControllerPack<Cs...>,
			field::FieldPack<FFs...>
		>
	> {
		using Controllers    = shared::internal::PackStorage<Cs...>;
		using Fields	     = shared::internal::PackStorage<FFs...>;
	public:
		using EnvT           = env::Environment<force::ForcePack<Fs...>, boundary::BoundaryPack<BCs...>, controller::ControllerPack<Cs...>, field::FieldPack<FFs...>>;
		using Container      = typename C::template impl<EnvT>;
		using ContainerFlags = container::internal::ContainerFlags;
		using BoundaryTable  = boundary::internal::BoundaryTable<typename EnvT::boundary_variant_t>;
		using ForceTable     = force::internal::ForceTable<EnvT>;
		using Interaction    = force::internal::InteractionInfo<typename EnvT::force_variant_t>;
		using Particle       = env::internal::Particle;
		using ParticleRef    = env::ParticleRef;
		using ParticleView   = env::ParticleView;
		using ParticleID     = env::internal::ParticleID;


		const env::Domain domain;

	private:
		// private constructor since System should only be creatable through build_system(...)
		System(
			const C & container_cfg,
			const ContainerFlags & container_flags,
			const env::Domain& domainIn,
			const std::vector<Particle> & particles,
			const BoundaryTable & boundaries,
			const Controllers & controllersIn,
			const Fields & fieldsIn,
			const UserToInternalMappings::TypeMap & usr_types_to_impl_types,
			const UserToInternalMappings::IdMap & usr_ids_to_impl_ids,
			std::vector<Interaction> & interaction_infos):
		domain(domainIn),
		container(container_cfg, container_flags),
		boundary_table(boundaries),
		controllers(controllersIn),
		fields(fieldsIn)
		{
			force_table.build(interaction_infos, usr_types_to_impl_types, usr_ids_to_impl_ids);
			container.init(force_table, domain);
			container.dispatch_build(particles);

			using SystemType = std::remove_cvref_t<decltype(*this)>;
			simulation_context = std::make_unique<internal::SimulationContextImpl<SystemType>>(*this);

			controllers.for_each_item([&](auto & controller) {controller.dispatch_init(context()); });
			fields.for_each_item([&](auto & field) {field.dispatch_init(context()); });
		}

		Container container;
		BoundaryTable boundary_table;
		ForceTable force_table;
		Controllers controllers;
		Fields fields;

		double time_ = 0;
		size_t step_ = 0;

		std::unique_ptr<SimulationContext> simulation_context;

		// System factory. Only valid way to create a System.
		template <container::IsContDecl Cont, class FPack, class BPack, class CPack, class FFPack>
		friend System<Cont, env::Environment<FPack, BPack, CPack, FFPack>>
		build_system(
			env::Environment<FPack, BPack, CPack, FFPack>& environment,
			 const Cont& container,
			 UserToInternalMappings* particle_mappings
		);

	public:
		[[nodiscard]] SimulationContext & context() const {
			return *simulation_context;
		}

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

				const Boundary & boundary = boundary_table.get_boundary(face);
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
			controllers.for_each_item([this](auto & controller) {
				if (controller.should_trigger(context())) {
					controller.dispatch_apply(context());
				}
			});
		}

		void apply_force_fields() {
			fields.for_each_item([this](auto & field) {
				for (size_t i = index_start(); i < index_end(); ++i) {
					field.dispatch_apply(env::RestrictedParticleRef(get_particle_by_index(i)));
				}
			});

			fields.for_each_item([this](auto & field) {
				field.dispatch_update(context());
			});
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

		[[nodiscard]] std::vector<size_t> collect_indices_in_region(const env::Box & region) {
			return container.dispatch_collect_indices_in_region(region);
		}

		[[nodiscard]] std::vector<size_t> collect_indices_in_region(const env::Domain & region) {
			return collect_indices_in_region(env::Box(region));
		}

		// returns the systems time
		[[nodiscard]] double time() const noexcept {
			return time_;
		}

		void update_time(const double dt) noexcept {
			time_ += dt;
		}

		[[nodiscard]] size_t step() const noexcept {
			return step_;
		}

		void increment_step() noexcept {
			++step_;
		}

		void reset_time() noexcept {
			time_ = 0;
			step_ = 0;
		}

		[[nodiscard]] size_t size(const env::ParticleState = env::ParticleState::ALL) const noexcept {
			// return container.size(state);
			// TODO implement this method properly
			return index_end() - index_start();
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


