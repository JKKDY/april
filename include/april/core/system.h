#pragma once

#include <memory>
#include "april/env/particle.h"
#include "april/forces/force.h"
#include "april/env/environment.h"
#include "april/containers/container.h"
#include "april/env/domain.h"
#include "april/boundaries/boundary_table.h"
#include "april/core/context.h"

namespace april::core {

	struct BuildInfo {
		std::unordered_map<env::ParticleType, env::internal::ParticleType> type_map;
		std::unordered_map<env::ParticleID, env::internal::ParticleID> id_map;
		env::Domain particle_box;
		env::Domain simulation_domain;
	};


	template <container::IsContDecl C, env::IsEnvironment Env>
	auto build_system(
		const Env & environment,
		const C& container,
		BuildInfo* build_info = nullptr
	);


	template <container::IsContDecl C, env::internal::IsEnvironmentTraits Traits>
	class System final {
		using Controllers    = typename Traits::controller_storage_t;
		using Fields	     = typename Traits::field_storage_t;
		using ForceTable     = typename Traits::force_table_t;
		using BoundaryTable  = typename Traits::boundary_table_t;
		using Container      = typename C::template impl<ForceTable>;
		using ContainerFlags = container::internal::ContainerFlags;
		using TypeInteraction= force::internal::TypeInteraction<typename Traits::force_variant_t>;
		using IdInteraction	 = force::internal::IdInteraction<typename Traits::force_variant_t>;

		using Particle       = env::internal::Particle;
		using ParticleRef    = env::ParticleRef;
		using ParticleView   = env::ParticleView;
		using ParticleID     = env::internal::ParticleID;

		// private constructor since System should only be creatable through build_system(...)
		System(
			const C & container_cfg,
			const ContainerFlags & container_flags,
			const env::Box & domainIn,
			const std::vector<Particle> & particles,
			const BoundaryTable & boundariesIn,
			const ForceTable & forcesIn,
			const Controllers & controllersIn,
			const Fields & fieldsIn
			):
		simulation_box(domainIn),
		boundary_table(boundariesIn),
		force_table(forcesIn),
		controllers(controllersIn),
		fields(fieldsIn),
		container(container_cfg, container_flags, domainIn, &force_table),
		simulation_context(std::make_unique<internal::SimulationContextImpl<decltype(*this)>>(*this))
		{
			container.dispatch_build(particles);
			controllers.for_each_item([&](auto & controller) {controller.dispatch_init(context()); });
			fields.for_each_item([&](auto & field) {field.dispatch_init(context()); });
		}

		env::Box simulation_box;
		BoundaryTable boundary_table;
		ForceTable force_table;
		Controllers controllers;
		Fields fields;
		Container container;

		double time_ = 0;
		size_t step_ = 0;

		std::unique_ptr<SimulationContext> simulation_context;

		// System factory. Only valid way to create a System.
		template <container::IsContDecl Cont, env::IsEnvironment Env>
		friend auto
		build_system(
			 const Env & environment,
			 const Cont& container,
			 BuildInfo * build_info
		);

	public:

		[[nodiscard]] env::Domain domain() const {
			return {simulation_box.min, simulation_box.extent};
		}

		[[nodiscard]] env::Box box() const {
			return simulation_box;
		}

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
			using Boundary = boundary::internal::CompiledBoundary<typename Traits::boundary_variant_t>;

			auto sim_box = box();

			for (boundary::Face face : boundary::all_faces) {

				const Boundary & boundary = boundary_table.get_boundary(face);
				std::vector<size_t> particle_ids = container.dispatch_collect_indices_in_region(boundary.region);

				if (boundary.topology.boundary_thickness >= 0) {
					for (auto p_idx : particle_ids) {
						env::internal::Particle & p = container.dispatch_get_particle_by_index(p_idx);
						boundary.apply(p, sim_box, face);

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
						const double y = diff[ax] < 0 ? sim_box.min[ax] : sim_box.max[ax];
						const double t = (y - p.old_position[ax]) / diff[ax];

						const vec3 intersection = t * diff + p.old_position;

						// and check if that point is on the domains surface
						auto [ax1, ax2] = non_face_axis(face);
						if (sim_box.max[ax1] >= intersection[ax1] && sim_box.min[ax1] <= intersection[ax1] &&
							sim_box.max[ax2] >= intersection[ax2] && sim_box.min[ax2] <= intersection[ax2]) {
							boundary.apply(p, sim_box, face);

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
			return collect_indices_in_region(env::Box(region.min_corner().value(), region.max_corner().value()));
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


	// Default: assume any type is not a System
	template<typename>
	inline constexpr bool is_system_v = false;

	// Specialization: mark all System<C, Env> instantiations as true
	template<container::IsContDecl C, env::internal::IsEnvironmentTraits Traits>
	inline constexpr bool is_system_v<System<C, Traits>> = true;

	// Concept: true if T (after removing cv/ref) is a System specialization
	template<typename T>
	concept IsSystem = is_system_v<std::remove_cvref_t<T>>;
}


