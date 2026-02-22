#pragma once

#include "april/env/environment.hpp"
#include "april/env/domain.hpp"
#include "april/containers/container.hpp"
#include "april/containers/batching/common.hpp"
#include "april/containers/kernel_traits.hpp"
#include "april/exec/policy.hpp"
#include "april/exec/particle_kernel.hpp"
#include "april/system/context.hpp"

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
		using ControllerStorage = Traits::controller_storage_t;
		using FieldStorage	    = Traits::field_storage_t;
		using ForceTable     	= Traits::force_table_t;
		using BoundaryTable  	= Traits::boundary_table_t;
		using ParticleData		= Traits::user_data_t;
	public:

		// ----------------
		// PUBLIC API TYPES
		// ----------------
		using SysContext = SystemContext<System>;
		using TrigContext = shared::TriggerContextImpl<System>;
		using Container = ContainerDecl::template impl<typename Traits::user_data_t>;
		using ParticleRec = Traits::particle_record_t;

		template<ParticleField M>
		using ParticleRef = Traits::template particle_ref_t<M>;

		template<ParticleField M>
		using ParticleView = Traits::template particle_view_t<M>;

		template<ParticleField M>
		using RestrictedParticleRef = Traits::template restricted_particle_ref_t<M>;


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
		template<ParticleField M>
		[[nodiscard]] env::ScalarParticleRef<M, ParticleData> at(const size_t index) {
			return particle_container.template at<M>(index);
		}

		template<ParticleField M>
		[[nodiscard]] env::ScalarParticleView<M, ParticleData> view(const size_t index) const {
			return particle_container.template view<M>(index);
		}

		template<ParticleField M>
		[[nodiscard]] env::ScalarRestrictedParticleRef<M, ParticleData> restricted_at(const size_t index) {
			return particle_container.template restricted_at<M>(index);
		}

		// ID ACCESSORS (stable)
		template<ParticleField M>
		[[nodiscard]] env::ScalarParticleRef<M, ParticleData> at_id(const ParticleID id) {
			return particle_container.template at_id<M>(id);
		}

		template<ParticleField M>
		[[nodiscard]] env::ScalarParticleView<M, ParticleData> view_id(const ParticleID id) const {
			return particle_container.template view_id<M>(id);
		}

		template<ParticleField M>
		[[nodiscard]] env::ScalarRestrictedParticleRef<M, ParticleData> restricted_at_id(const ParticleID id) {
			return particle_container.template restricted_at_id<M>(id);
		}


		// -----------
		// ID INDEXING
		// -----------
		// get the lowest particle id
		[[nodiscard]] ParticleID min_id() const noexcept{
			return particle_container.min_id();
		}

		// get the largest particle id
		[[nodiscard]] ParticleID max_id() const noexcept {
			return particle_container.max_id();
		}

		// check if particle id is valid
		[[nodiscard]] bool contains_id(const ParticleID id) const noexcept {
			// with min_id anc max_id can write a stable for loop:
			//	for (auto id = min_id(); id < max_id; id++) { if !contains_id(id) continue; ...}
			return particle_container.invoke_contains_id(id);
		}

		// convert id to index
		[[nodiscard]] size_t id_to_index(const ParticleID id) const noexcept {
			return particle_container.invoke_id_to_index(id);
		}


		// -------
		// QUERIES
		// -------
		[[nodiscard]] size_t size(const ParticleState = ParticleState::ALL) const noexcept {
			// TODO implement this method (system::size) properly
			return particle_container.invoke_particle_count();
		}
		[[nodiscard]] size_t capacity() const noexcept {
			return particle_container.capacity();
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
		template<
			ParticleField M,
			ParallelPolicy P=ParallelPolicy::Serial,
			VectorPolicy V=VectorPolicy::Auto,
			exec::IsKernel Kernel>
		void for_each_particle(Kernel && func, ParticleState state = ParticleState::ALL) {
			particle_container.template for_each_particle<M, P, V, Kernel>(std::forward<Kernel>(func), state);
		}

		template<
			ParticleField M,
			ParallelPolicy P=ParallelPolicy::Serial,
			VectorPolicy V=VectorPolicy::Auto,
			exec::IsKernel Kernel>
		void for_each_particle_view(Kernel && func, ParticleState state = ParticleState::ALL) const {
			particle_container.template for_each_particle_view<M, P, V, Kernel>(std::forward<Kernel>(func), state);
		}

		template<ParticleField M, typename Mapper, typename T, typename Reducer = std::plus<T>>
		[[nodiscard]] T reduce(
			T initial_value,
			Mapper&& map_func, // map particle -> T
			Reducer&& reduce_func = {},
			ParticleState state = ParticleState::ALIVE
		) const {
			return particle_container.template invoke_reduce<M>(initial_value, map_func, reduce_func, state);
		}



		template<typename Func>
		void for_each_interaction_batch(Func && func) { // func(batch, bcp)
			particle_container.invoke_for_each_interaction_batch(std::forward<Func>(func));
		}

		template<ParticleField M, typename Func>
		void for_each_interaction_pair(Func && func) { // func(particle, particle, dist)
			auto update_batch = [&]<container::IsBatch Batch, /*TODO container::IsBCP*/ typename BCP>(const Batch& batch, BCP && apply_bcp) {
				auto kernel = [&](auto && p1, auto && p2) {
					vec3 r = {};
					if constexpr (env::has_field_v<M, ParticleField::position> &&
						std::is_same_v<std::decay_t<BCP>, container::NoBatchBCP>) {
						r = p2.position - p1.position;
						} else {
							r = apply_bcp(p2.position - p1.position);
						}

					func(p1, p2, r);
				};

				execute_batch_kernel<M, ParallelPolicy::Serial, VectorPolicy::Scalar>(batch, kernel);
			};

			for_each_interaction_batch(update_batch);
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

		void notify_moved_id(const std::vector<ParticleID> & ids) {
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
		[[nodiscard]] SysContext & context() {
			return system_context;
		}
		[[nodiscard]] TrigContext & trigger_context() {
			return trig_context;
		}

		[[nodiscard]] const SysContext & context() const {
			return system_context;
		}
		[[nodiscard]] const TrigContext & trigger_context() const {
			return trig_context;
		}


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

		template<ParticleField M, ParallelPolicy P, VectorPolicy V, container::IsBatch Batch, typename Kernel>
		void execute_batch_kernel(const Batch& batch, Kernel&& kernel);
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

} // namespace april::core

#include "april/system/internal/system_impl.hpp"

