#pragma once

#include "april/particle/fields.hpp"
#include "april/forces/force.hpp"
#include "april/env/environment.hpp"
#include "april/containers/container.hpp"
#include "april/env/domain.hpp"
#include "april/system/context.hpp"
#include "april/containers/batching.hpp"
#include "april/boundaries/boundary.hpp"



namespace april::core {
	struct BuildInfo;

	template <class C, env::internal::IsEnvironmentTraits Traits>
	requires container::IsContainerDecl<C, Traits>
	class System;


	template <class Container, env::IsEnvironment EnvT>
	requires container::IsContainerDecl<Container, typename EnvT::traits>
	System<Container, typename EnvT::traits> build_system(
		const EnvT & environment,
		const Container & container_config,
		BuildInfo * build_info = nullptr
	);


	template <class ContainerDecl, env::internal::IsEnvironmentTraits Traits>
	requires container::IsContainerDecl<ContainerDecl, Traits>
	class System final {
		// --------------
		// INTERNAL TYPES
		// --------------
		using ControllerStorage = typename Traits::controller_storage_t;
		using FieldStorage	    = typename Traits::field_storage_t;
		using ForceTable     	= typename Traits::force_table_t;
		using BoundaryTable  	= typename Traits::boundary_table_t;

	public:

		// ----------------
		// PUBLIC API TYPES
		// ----------------
		using SysContext = SystemContext<System>;
		using TrigContext = shared::TriggerContextImpl<System>;
		using Container = typename ContainerDecl::template impl<typename Traits::user_data_t>;
		using ParticleRec = typename Traits::particle_record_t;

		template<env::FieldMask M>
		using ParticleRef = typename Traits::template particle_ref_t<M>;

		template<env::FieldMask M>
		using ParticleView = typename Traits::template particle_view_t<M>;

		template<env::FieldMask M>
		using RestrictedParticleRef = typename Traits::template restricted_particle_ref_t<M>;


		// -----------------
		// LIFECYCLE & STATE
		// -----------------
		[[nodiscard]] double time() const noexcept { return time_; }
		[[nodiscard]] size_t step() const noexcept { return step_; }
		[[nodiscard]] env::Domain domain() const { return {simulation_box.min, simulation_box.extent}; }
		[[nodiscard]] env::Box box() const { return simulation_box; }

		void update_time(const double dt) noexcept { time_ += dt; }
		void increment_step() noexcept { ++step_; }
		void reset_time() noexcept { time_ = 0; step_ = 0; }


		// ------------------
		// PARTICLE ACCESSORS
		// ------------------
		// "at" implies mutable access, "view" implies read-only, "restricted_at" implies only force mutable

		// INDEX ACCESSORS (fast)
		template<env::FieldMask M>
		[[nodiscard]] auto at(const size_t index) {
			return particle_container.template at<M>(index);
		}

		template<env::FieldMask M>
		[[nodiscard]] auto view(const size_t index) const {
			return particle_container.template view<M>(index);
		}

		template<env::FieldMask M>
		[[nodiscard]] auto restricted_at(const size_t index) {
			return particle_container.template restricted_at<M>(index);
		}

		// ID ACCESSORS (stable)
		template<env::FieldMask M>
		[[nodiscard]] auto at_id(const env::ParticleID id) {
			return particle_container.template at_id<M>(id);
		}

		template<env::FieldMask M>
		[[nodiscard]] auto view_id(const env::ParticleID id) const {
			return particle_container.template view_id<M>(id);
		}

		template<env::FieldMask M>
		[[nodiscard]] auto restricted_at_id(const env::ParticleID id) {
			return particle_container.template restricted_at_id<M>(id);
		}


		// -----------
		// ID INDEXING
		// -----------
		// get the lowest particle id
		[[nodiscard]] env::ParticleID min_id() const noexcept{
			return particle_container.min_id();
		}

		// get the largest particle id
		[[nodiscard]] env::ParticleID max_id() const noexcept {
			return particle_container.max_id();
		}

		// check if particle id is valid
		[[nodiscard]] bool contains(const env::ParticleID id) const noexcept {
			return particle_container.invoke_contains(id);
		}

		// convert id to index
		[[nodiscard]] size_t id_to_index(const env::ParticleID id) const noexcept {
			return particle_container.invoke_id_to_index(id);
		}


		// -------
		// QUERIES
		// -------
		[[nodiscard]] size_t size(const env::ParticleState = env::ParticleState::ALL) const noexcept {
			// TODO implement this method (system::size) properly
			return particle_container.invoke_particle_count();
		}

		[[nodiscard]] std::vector<size_t> query_region(const env::Box & region) const {
			return particle_container.invoke_collect_indices_in_region(region);
		}

		[[nodiscard]] std::vector<size_t> query_region(const env::Domain & region) const {
			return query_region(env::Box(region.min_corner().value(), region.max_corner().value()));
		}


		// --------------
		// FUNCTIONAL OPS
		// --------------
		template<env::FieldMask M, typename Func, bool parallelize=false>
		void for_each_particle(Func && func, env::ParticleState state = env::ParticleState::ALL) {
			particle_container.template invoke_for_each_particle<M, Func, parallelize>(std::forward<Func>(func), state);
		}

		template<typename Func>
		void for_each_interaction_batch(Func && func) {
			particle_container.invoke_for_each_interaction_batch(std::forward<Func>(func));
		}

		template<env::FieldMask M, typename T, typename Mapper, typename Reducer = std::plus<T>>
		[[nodiscard]] T invoke_reduce(
			T initial_value,
			Mapper&& map_func,
			Reducer&& reduce_func = {},
			env::ParticleState state = env::ParticleState::ALIVE
		) const {
			return particle_container.template invoke_reduce<M>(initial_value, map_func, reduce_func, state);
		}


		// -----------------
		// STRUCTURE UPDATES
		// -----------------
		// call to register particle movements. This may cause container internals to change/be rebuilt
		void rebuild_structure() {
			particle_container.invoke_rebuild_structure();
		}

		void notify_moved(const std::vector<size_t> & indices) {
			particle_container.invoke_notify_moved(indices);
		}

		void notify_moved_id(const std::vector<env::ParticleID> & ids) {
			std::vector<size_t> indices;
			indices.reserve(ids.size());

			for (auto id : ids) {
				indices.push_back(particle_container.invoke_id_to_index(id));
			}

			particle_container.invoke_notify_moved(indices);
		}


		// -------
		// PHYSICS
		// -------
		void update_forces();
		void apply_boundary_conditions();
		void apply_controllers();
		void apply_force_fields();
		void update_all_components();


		// --------
		// CONTEXTS
		// --------
		[[nodiscard]] SysContext & context() { return system_context; }
		[[nodiscard]] TrigContext & trigger_context() { return trig_context; }

		[[nodiscard]] const SysContext & context() const { return system_context; }
		[[nodiscard]] const TrigContext & trigger_context() const { return trig_context; }


	private:
		env::Box simulation_box; // TODO rename to domain
		BoundaryTable boundary_table;
		ForceTable force_table;
		ControllerStorage controllers;
		FieldStorage fields;
		Container particle_container;

		std::vector<size_t> particles_to_update_buffer;

		double time_ = 0;
		size_t step_ = 0;

		SystemContext<System> system_context;
		shared::TriggerContextImpl<System> trig_context;


		// private constructor since System should only be creatable through build_system(...)
		System(
			const ContainerDecl& container_cfg,
			const container::internal::ContainerCreateInfo & container_info,
			const env::Box& domain_in,
			const std::vector<ParticleRec>& particles,
			const BoundaryTable& boundaries_in,
			const ForceTable& forces_in,
			const ControllerStorage& controllers_in,
			const FieldStorage& fields_in
		)
			: simulation_box(domain_in),
			  boundary_table(boundaries_in),
			  force_table(forces_in),
			  controllers(controllers_in),
			  fields(fields_in),
			  particle_container(Container(container_cfg, container_info)),
			  system_context(*this),
			  trig_context(*this)
		{
			particle_container.invoke_build(particles);
			controllers.for_each_item([&](auto& c) { c.dispatch_init(context()); });
			fields.for_each_item([&](auto& f) { f.dispatch_init(context()); });
		}

		// Friend factory: only entry point for constructing a System
		// Avoids exposing constructor internals publicly.
		template <class Container, env::IsEnvironment EnvT>
		requires container::IsContainerDecl<Container, typename EnvT::traits>
		friend System<Container, typename EnvT::traits>
		build_system(
			 const EnvT & environment,
			 const Container& container,
			 BuildInfo * build_info
		);
	};




	// Default: assume any type is not a System
	template<typename>
	inline constexpr bool is_system_v = false;

	// Specialization: mark all System<C, Env> instantiations as true
	template<class C, env::internal::IsEnvironmentTraits Traits>
	inline constexpr bool is_system_v<System<C, Traits>> = true;

	// Concept: true if T (after removing cv/ref) is a System specialization
	template<typename T>
	concept IsSystem = is_system_v<std::remove_cvref_t<T>>;




	// ---- Implementations -----
	template <class C, env::internal::IsEnvironmentTraits Traits> requires container::IsContainerDecl<C, Traits>
	void System<C, Traits>::update_forces() {

		// batch update lambda. passed into container::for_each_interaction_batch
		auto update_batch = [&]<batching::IsBatch Batch, batching::IsBCP BCP>(const Batch& batch, BCP && apply_bcp) {
			using Par = batching::ParallelPolicy;
			using Upd = batching::UpdatePolicy;
			using Cmp = batching::ComputePolicy;

			auto apply_batch_update =  [&] <force::IsForce ForceT> (const ForceT & force) {
				constexpr env::FieldMask fields = ForceT::fields;
				constexpr env::FieldMask upd_fields = env::Field::force | env::Field::position;
				constexpr env::FieldMask all_fields = upd_fields | fields;

				//---------------
				// PHYSICS KERNEL
				//---------------
				auto kernel = [&](auto & p1, auto & p2) {
					vec3 r;

					if constexpr (std::is_same_v<std::decay_t<BCP>, batching::NoBatchBCP>) {
						r = p2.position - p1.position;
					} else {
						r = apply_bcp(p2.position - p1.position);
					}

					if (r.norm_squared() > force.cutoff2()) {
						return;
					}

					const vec3 f = force(p1.to_view(), p2.to_view(), r);

					if constexpr (Batch::update_policy == Upd::Atomic) {
						static_assert(Batch::update_policy == Upd::Atomic, "atomic force update not implemented yet");
					} else {
						p1.force += f;
						p2.force -= f;
					}
				};

				// Execution
				auto execute_atom = [&](const auto& atom) {
					if constexpr (Batch::compute_policy == Cmp::Scalar) {
						atom.template for_each_pair<all_fields>(kernel);
					} else {
						static_assert(Batch::compute_policy != Cmp::Scalar, "Vectorization not implemented yet");
					}
				};

				auto execute_range_serial = [&](const auto& range) {
					for (const auto& atom : range) {
						execute_atom(atom);
					}
				};

				auto execute_range_parallel = [&](const auto&) {
					std::unreachable();
				};

				// Routing
				if constexpr (batching::IsBatchAtom<Batch>) {
					execute_atom(batch);
				}
				else if constexpr (batching::IsBatchAtomRange<Batch>) {
					if constexpr (Batch::parallel_policy == Par::Inner) {
						static_assert(Batch::parallel_policy == Cmp::Scalar, "Vectorization not implemented yet");
						execute_range_parallel(batch);
					} else {
						execute_range_serial(batch);
					}
				}

			};

			auto [t1, t2] = batch.types;
			force_table.dispatch(t1, t2, apply_batch_update);
		};

		auto reset_force = [](auto && p) {
			p.force = {};
		};

		particle_container.template invoke_for_each_particle<+env::Field::force>(reset_force);
		particle_container.invoke_for_each_interaction_batch(update_batch);


		// handle id interactions
		auto update_global_batch = [&](const auto & batch) {

			auto apply_batch_update = [&] <force::IsForce ForceT> (const ForceT & force) {
				constexpr env::FieldMask fields = ForceT::fields;
				constexpr env::FieldMask upd_fields = env::Field::force | env::Field::position;
				constexpr env::FieldMask all_fields = upd_fields | fields;

				for (const auto & [id1, id2] : batch.pairs) {
					auto && p1 = particle_container.template restricted_at_id<all_fields>(id1);
					auto && p2 = particle_container.template restricted_at_id<all_fields>(id2);

					vec3 r = p2.position - p1.position;
					const vec3 f = force(p1.to_view(), p2.to_view(), r);

					p1.force += f;
					p2.force -= f;
				}
			};

			force_table.dispatch_id(batch.id1, batch.id2, apply_batch_update);
		};

		particle_container.invoke_for_each_topology_batch(update_global_batch);
	}


	template <class C, env::internal::IsEnvironmentTraits Traits> requires container::IsContainerDecl<C, Traits>
	void System<C, Traits>::apply_boundary_conditions() {
		particles_to_update_buffer.clear();

		const env::Box domain_box = this->box();

		for (boundary::Face face : boundary::all_faces) {

			const auto& compiled_boundary = boundary_table[face];

			std::vector<size_t> particle_ids = particle_container.invoke_collect_indices_in_region(compiled_boundary.boundary_region);

			auto boundary_condition_inside = [&]<typename B>(const B & bc) {
				constexpr env::FieldMask M = std::decay_t<B>::fields;

				for (auto p_idx : particle_ids) {
					auto p = particle_container.template at<M>(p_idx);
					bc.apply(p, domain_box, face);

					if (compiled_boundary.topology.may_change_particle_position) {
						particles_to_update_buffer.push_back(p_idx);
					}
				}
			};

			auto boundary_condition_outside = [&]<typename B>(const B & bc) {
				static constexpr env::FieldMask detect_mask = env::Field::position | env::Field::old_position;
				constexpr env::FieldMask M = std::decay_t<B>::fields | detect_mask;

				for (auto p_idx : particle_ids) {
					auto particle = particle_container.template at<M>(p_idx);

					// make sure the particle exited through the current boundary face
					// solve for intersection of the particles path with the boundary face
					// with the equation y = t * diff + p where:
					// diff is the path traveled, p is the particles starting position and y is the face

					const int ax = axis_of_face(face);
					const vec3 diff = particle.position - particle.old_position;
					const double y = diff[ax] < 0 ? domain_box.min[ax] : domain_box.max[ax];
					const double t = (y - particle.old_position[ax]) / diff[ax];

					const vec3 intersection = t * diff + particle.old_position;

					// and check if that point is on the domains surface
					auto [ax1, ax2] = non_face_axis(face);
					if (domain_box.max[ax1] >= intersection[ax1] && domain_box.min[ax1] <= intersection[ax1] &&
						domain_box.max[ax2] >= intersection[ax2] && domain_box.min[ax2] <= intersection[ax2]) {

						bc.apply(particle, domain_box, face);

						if (compiled_boundary.topology.may_change_particle_position) {
							particles_to_update_buffer.push_back(p_idx);
						}
					}
				}
			};

			if (compiled_boundary.topology.boundary_thickness >= 0) {
				compiled_boundary.dispatch(boundary_condition_inside);
			} else {
				compiled_boundary.dispatch(boundary_condition_outside);
			}
		}

		if (!particles_to_update_buffer.empty()) {
			particle_container.invoke_notify_moved(particles_to_update_buffer);
		}
	}

	template <class C, env::internal::IsEnvironmentTraits Traits> requires container::IsContainerDecl<C, Traits>
	void System<C, Traits>::apply_controllers() {
		controllers.for_each_item([this](auto & controller) {
			if (controller.should_trigger(trig_context)) {
				controller.dispatch_apply(system_context);
			}
		});
	}

	template <class C, env::internal::IsEnvironmentTraits Traits> requires container::IsContainerDecl<C, Traits>
	void System<C, Traits>::apply_force_fields() {
		fields.for_each_item([this]<typename F>(F & field) {
			for (size_t i = 0; i < size(); ++i) {
				constexpr env::FieldMask M = F::fields;
				auto restricted = particle_container.template restricted_at<M>(i);
				field.dispatch_apply(restricted);
			}
		});
	}

	template <class C, env::internal::IsEnvironmentTraits Traits> requires container::IsContainerDecl<C, Traits>
	void System<C, Traits>::update_all_components() {
		fields.for_each_item([this](auto & field) {
			field.template dispatch_update<System>(system_context);
		});

		controllers.for_each_item([this](auto & controller) {
			controller.template dispatch_update<System>(system_context);
		});
	}
}


