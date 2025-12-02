#pragma once

#include "april/particle/fields.h"
#include "april/forces/force.h"
#include "april/env/environment.h"
#include "april/containers/container.h"
#include "april/env/domain.h"
#include "april/system/context.h"
#include "april/containers/batch.h"


// TODO: add deleter callbacks to detect if system goes out of scope too early (e.g. prior to integrator.run)


namespace april::core {
	struct BuildInfo;

	template <class C, env::internal::IsEnvironmentTraits Traits>
	requires container::IsContainerDecl<C, Traits>
	class System;


	template <class Container, env::IsEnvironment EnvT>
	requires container::IsContainerDecl<Container, typename EnvT::traits>
	System<Container, typename EnvT::traits> build_system(
		const EnvT & environment,
		const Container& container,
		BuildInfo * build_info = nullptr
	);


	template <class C, env::internal::IsEnvironmentTraits Traits>
	requires container::IsContainerDecl<C, Traits>
	class System final {
		using Controllers    	= typename Traits::controller_storage_t;
		using Fields	     	= typename Traits::field_storage_t;
		using ForceTable     	= typename Traits::force_table_t;
		using BoundaryTable  	= typename Traits::boundary_table_t;

		using Container      	= typename C::template impl<typename Traits::user_data_t>;
		using ContainerFlags 	= container::internal::ContainerFlags;
		using TypeInteraction	= force::internal::TypeInteraction<typename Traits::force_variant_t>;
		using IdInteraction	 	= force::internal::IdInteraction<typename Traits::force_variant_t>;

		using SysContext 		= SystemContext<System>;
		using TrigContext		= shared::TriggerContextImpl<System>;

	public:
		using ParticleRec = typename Traits::particle_record_t;
		using user_data_t = typename Traits::user_data_t;
		template<env::FieldMask M> using ParticleRef    		= typename Traits::template particle_ref_t<M>;
		template<env::FieldMask M> using ParticleView   		= typename Traits::template particle_view_t<M>;
		template<env::FieldMask M> using RestrictedParticleRef	= typename Traits::template restricted_particle_ref_t<M>;

		[[nodiscard]] env::Domain domain() const { return {simulation_box.min, simulation_box.extent}; }
		[[nodiscard]] env::Box box() const { return simulation_box; }

		[[nodiscard]] double time() const noexcept { return time_; }
		[[nodiscard]] size_t step() const noexcept { return step_; }
		void update_time(const double dt) noexcept { time_ += dt; }
		void increment_step() noexcept { ++step_; }
		void reset_time() noexcept { time_ = 0; step_ = 0; }

		[[nodiscard]] size_t size(const env::ParticleState = env::ParticleState::ALL) const noexcept {
			// return container.size(state);
			// TODO implement this method (system::size) properly
			return index_end() - index_start();
		}



		[[nodiscard]] std::vector<size_t> collect_indices_in_region(const env::Box & region) {
			return container.dispatch_collect_indices_in_region(region);
		}

		[[nodiscard]] std::vector<size_t> collect_indices_in_region(const env::Domain & region) {
			return collect_indices_in_region(env::Box(region.min_corner().value(), region.max_corner().value()));
		}

		// call to register particle movements. This may cause container internals to change/be rebuilt
		void register_all_particle_movements() {
			container.dispatch_register_all_particle_movements();
		}

		void register_particle_movement(env::ParticleID id) {
			size_t idx = container.id_to_index(id);
			container.dispatch_register_particle_movement(idx);
		}

		// call to update all pairwise forces between particles
		void update_forces();

		// call to apply boundary conditions to all particles.
		// should not be called before register_particle_movements
		void apply_boundary_conditions();

		void apply_controllers();

		void apply_force_fields();

		void update_all_components();



		[[nodiscard]] SysContext & context() { return system_context; }
		[[nodiscard]] const SysContext & context() const { return system_context; }

		[[nodiscard]] TrigContext & trigger_context() { return trig_context; }
		[[nodiscard]] const TrigContext & trigger_context() const { return trig_context; }


		// get a particle reference by its id. Usually slower than getting it by its index.
		// Useful for stable iterations and accessing a specific particle
		template<env::FieldMask M>
		[[nodiscard]] ParticleRef<M> get_particle_by_id(const env::ParticleID id) noexcept {
			return ParticleRef<M>(container.dispatch_get_fetcher_by_id(id));
		}

		template<env::FieldMask M>
		[[nodiscard]] ParticleView<M> get_particle_by_id(const env::ParticleID id) const noexcept {
			return ParticleView<M>(container.dispatch_get_fetcher_by_id(id));
		}

		// get the lowest particle id
		[[nodiscard]] env::ParticleID id_start() const noexcept{
			return container.dispatch_id_start();
		}

		// get the highest particle id
		[[nodiscard]] env::ParticleID id_end() const noexcept {
			return container.dispatch_id_end();
		}


		// get a particle by its container specific id. useful for non-stable (but fast) iteration over particles
		template<env::FieldMask M>
		[[nodiscard]] ParticleRef<M> get_particle_by_index(const size_t index) noexcept {
			return ParticleRef<M>{container.dispatch_get_fetcher_by_index(index)};
		}

		template<env::FieldMask M>
		[[nodiscard]] ParticleView<M> get_particle_by_index(const size_t index) const noexcept {
			return ParticleView<M>{container.dispatch_get_fetcher_by_index(index)};
		}

		// get the first particle index (usually 0)
		[[nodiscard]] size_t index_start() const noexcept {
			return container.dispatch_index_start();
		}

		// get the last particle index (usually n-1 with n = #particles)
		[[nodiscard]] size_t index_end() const noexcept {
			return container.dispatch_index_end();
		}


	private:
		env::Box simulation_box;
		BoundaryTable boundary_table;
		ForceTable force_table;
		Controllers controllers;
		Fields fields;
		Container container;

		std::vector<size_t> particles_to_update_buffer;

		double time_ = 0;
		size_t step_ = 0;

		SystemContext<System> system_context;
		shared::TriggerContextImpl<System> trig_context;


		// private constructor since System should only be creatable through build_system(...)
		System(
			const C& container_cfg,
			const ContainerFlags& container_flags,
			const env::Box& domain_in,
			const std::vector<ParticleRec>& particles,
			const BoundaryTable& boundaries_in,
			const ForceTable& forces_in,
			const Controllers& controllers_in,
			const Fields& fields_in
		)
			: simulation_box(domain_in),
			  boundary_table(boundaries_in),
			  force_table(forces_in),
			  controllers(controllers_in),
			  fields(fields_in),
			  container(Container(container_cfg, container_flags, domain_in)),
			  system_context(*this),
			  trig_context(*this)
		{
			container.dispatch_build(particles);
			controllers.for_each_item([&](auto& c) { c.dispatch_init(context()); });
			fields.for_each_item([&](auto& f) { f.dispatch_init(context()); });

		}

		// Friend factory: only entry point for constructing a System
		// Avoids exposing constructor internals publicly.
		template <class Cont, env::IsEnvironment Env>
		requires container::IsContainerDecl<Cont, typename Env::traits>
		friend System<Cont, typename Env::traits>
		build_system(
			 const Env & environment,
			 const Cont& container,
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
			auto update_batch = [&]<container::IsBatch Batch, container::IsBCP BCP>(const Batch& batch, BCP && apply_bcp) {
				auto apply_batch_update =  [&] <force::IsForce ForceT> (const ForceT & force) {
					using Sym = container::BatchSymmetry;
					using Par = container::ParallelPolicy;
					using Upd = container::UpdatePolicy;

					// PHYSICS KERNEL
					auto update_force = [&](auto & p1, auto & p2) {
						vec3 diff;

						if constexpr (std::is_same_v<std::decay_t<BCP>, container::NoBatchBCP>) {
							diff = p2.position() - p1.position();
						} else {
							diff = apply_bcp(p2.position() - p1.position());
						}

						if (diff.norm_squared() > force.cutoff2()) {
							return;
						}

						// TODO change to use Restricted Refs instead of fetchers
						const vec3 f = force(
							env::internal::ConstParticleRecordFetcher(p1),
							env::internal::ConstParticleRecordFetcher(p2),
							diff);

						if constexpr (batch.update_policy == Upd::Atomic) {
							throw std::logic_error("atomic force update not implemented yet");
						} else {
							p1.force() += f;
							p2.force() -= f;
						}
					};


					// INNER LOOP
					auto run_symmetric_parallel = [&](const auto&) {
						throw std::logic_error("parallelization not implemented");
					};

					auto run_asymmetric_parallel = [&](const auto&) {
						throw std::logic_error("parallelization not implemented");
					};

					auto run_symmetric_serial = [&](const auto& batch_item) {
						for (size_t i = 0; i < batch_item.indices.size(); ++i) {
							auto && p1 = container.dispatch_get_fetcher_by_index(batch_item.indices[i]);

							for (size_t j = i + 1; j < batch_item.indices.size(); ++j) {
								auto && p2 = container.dispatch_get_fetcher_by_index(batch_item.indices[j]);

								update_force(p1, p2);
							}
						}
					};

					auto run_asymmetric_serial = [&](const auto& batch_item) {
						for (size_t i = 0; i < batch_item.indices1.size(); ++i) {
							auto && p1 = container.dispatch_get_fetcher_by_index(batch_item.indices1[i]);

							for (size_t j = 0; j < batch_item.indices2.size(); ++j) {
								auto && p2 = container.dispatch_get_fetcher_by_index(batch_item.indices2[j]);

								update_force(p1, p2);
							}
						}
					};


					// INNER LOOP DISPATCHING
					auto run_symmetric_inner_loop = [&](const auto& batch_item){
						if constexpr (batch.parallel_policy == Par::InnerLoop) {
							run_symmetric_parallel(batch_item);
						} else {
							run_symmetric_serial(batch_item);
						}
					};

					auto run_asymmetric_inner_loop = [&](const auto& batch_item){
						if constexpr (batch.parallel_policy == Par::InnerLoop) {
							run_asymmetric_parallel(batch_item);
						} else {
							run_asymmetric_serial(batch_item);
						}
					};


					// OUTER LOOP GENERIC STRATEGY
					auto execute_strategy = [&](auto&& run_inner) {
						if constexpr (container::IsChunkedBatch<Batch>) {
							// Chunked Path
							if constexpr (batch.parallel_policy == Par::Chunks) {
								throw std::logic_error("parallelized chunked processing not implemented yet");
							} else {
								for (const auto & chunk : batch.chunks) {
									run_inner(chunk);
								}
							}
						} else {
							run_inner(batch);
						}
					};


					// OUTER LOOP DISPATCHING
					if constexpr (batch.symmetry == Sym::Symmetric) {
						execute_strategy(run_symmetric_inner_loop);
					} else if constexpr (batch.symmetry == Sym::Asymmetric) {
						execute_strategy(run_asymmetric_inner_loop);
					}
				};

				auto [t1, t2] = batch.types;
				force_table.dispatch(t1, t2, apply_batch_update);
			};


			container.dispatch_prepare_force_update();
			// container.prepare_batches();
			container.dispatch_for_each_batch(update_batch);
		}


	template <class C, env::internal::IsEnvironmentTraits Traits> requires container::IsContainerDecl<C, Traits>
	void System<C, Traits>::apply_boundary_conditions() {
		using Boundary = boundary::internal::CompiledBoundary<typename Traits::boundary_variant_t>;

		auto sim_box = box();
		particles_to_update_buffer.clear();

		for (boundary::Face face : boundary::all_faces) {

			const Boundary & boundary = boundary_table.get_boundary(face);
			std::vector<size_t> particle_ids = container.dispatch_collect_indices_in_region(boundary.region);

			if (boundary.topology.boundary_thickness >= 0) {
				for (auto p_idx : particle_ids) {

					auto p = container.dispatch_get_fetcher_by_index(p_idx);
					boundary.apply(p, sim_box, face);

					if (boundary.topology.may_change_particle_position) {
						particles_to_update_buffer.push_back(p_idx);
					}
				}
			} else {
				for (auto p_idx : particle_ids) {
					static constexpr env::FieldMask M = env::Field::position | env::Field::old_position;
					auto particle = get_particle_by_index<M>(p_idx);

					// make sure the particle exited through the current boundary face
					// solve for intersection of the particles path with the boundary face
					// with the equation y = t * diff + p where:
					// diff is the path traveled, p is the particles starting position and y is the face
					const int ax = axis_of_face(face);
					const vec3 diff = particle.position - particle.old_position;
					const double y = diff[ax] < 0 ? sim_box.min[ax] : sim_box.max[ax];
					const double t = (y - particle.old_position[ax]) / diff[ax];

					const vec3 intersection = t * diff + particle.old_position;

					// and check if that point is on the domains surface
					auto [ax1, ax2] = non_face_axis(face);
					if (sim_box.max[ax1] >= intersection[ax1] && sim_box.min[ax1] <= intersection[ax1] &&
						sim_box.max[ax2] >= intersection[ax2] && sim_box.min[ax2] <= intersection[ax2]) {

						auto p = container.dispatch_get_fetcher_by_index(p_idx);
						boundary.apply(p, sim_box, face);

						if (boundary.topology.may_change_particle_position) {
							particles_to_update_buffer.push_back(p_idx);
						}
					}
				}
			}
		}

		container.dispatch_register_particle_movement(particles_to_update_buffer);

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
			for (size_t i = index_start(); i < index_end(); ++i) {
				constexpr env::FieldMask M = F::fields;
				RestrictedParticleRef<M> restricted (container.dispatch_get_fetcher_by_id(i));
				field.template dispatch_apply<user_data_t>(restricted);
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


