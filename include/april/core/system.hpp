#pragma once

#include "april/core/environment.hpp"
#include "april/core/domain.hpp"
#include "april/containers/container.hpp"
#include "april/containers/batching/common.hpp"
#include "april/exec/policy.hpp"
#include "april/exec/particle_kernel.hpp"
#include "april/core/context.hpp"
#include "april/exec/executor.hpp"

namespace april {
	struct BuildInfo;

	template <class C, core::internal::IsEnvironmentTraits Traits>
	requires container::IsContainerDecl<C, Traits>
	class System;


	template <class Container, core::IsEnvironment EnvT>
	requires container::IsContainerDecl<Container, typename EnvT::traits>
	System<Container, typename EnvT::traits> build_system(
		const EnvT & environment,
		const Container & container_config,
		BuildInfo * build_info = nullptr
	);

	template <class ContainerDecl, core::internal::IsEnvironmentTraits Traits>
	requires container::IsContainerDecl<ContainerDecl, Traits>
	class System final {
		// --------------
		// INTERNAL TYPES
		// --------------
		using ControllerStorage = Traits::controller_storage_t;
		using FieldStorage	    = Traits::field_storage_t;
		using ForceTable     	= Traits::force_table_t;
		using BoundaryTable  	= Traits::boundary_table_t;
		using ParticleAttributes= Traits::particle_attributes_t;
	public:

		// ----------------
		// PUBLIC API TYPES
		// ----------------
		using SysContext = core::SystemContext<System>;
		using TrigContext = utility::internal::TriggerContextImpl<System>;
		using Container = ContainerDecl::template impl<ParticleAttributes>;
		using ParticleRec = Traits::particle_record_t;


		// -----------------
		// LIFECYCLE & STATE
		// -----------------
		[[nodiscard]] double time() const noexcept { return time_; }
		[[nodiscard]] size_t step() const noexcept { return step_; }
		[[nodiscard]] Domain domain() const { return {simulation_box.min, simulation_box.extent}; }
		[[nodiscard]] core::Box box() const { return simulation_box; }

		void update_time(const double dt) noexcept { time_ += dt; }
		void increment_step() noexcept { ++step_; }
		void reset_time() noexcept { time_ = 0; step_ = 0; }


		// ------------------
		// PARTICLE ACCESSORS
		// ------------------
		// "at" implies mutable access, "view" implies read-only,

		// INDEX ACCESSORS (fast)
		template<ParticleField Read, ParticleField Write = Read>
		[[nodiscard]] auto at(const size_t index) {
			// type: particle::internal::ScalarParticleRef<Read, Write, ParticleData>
			return particle_container.template at<Read, Write>(index);
		}

		template<ParticleField Read>
		[[nodiscard]] auto view(const size_t index) const {
			// type: particle::internal::ScalarParticleRef<Read, ParticleField::none, ParticleData>
			return particle_container.template view<Read>(index);
		}

		// ID ACCESSORS (stable)
		template<ParticleField Read, ParticleField Write = Read>
		[[nodiscard]] auto at_id(const ParticleID id) {
			// type: particle::internal::ScalarParticleRef<Read, Write, ParticleData>
			return particle_container.template at_id<Read, Write>(id);
		}

		template<ParticleField Read>
		[[nodiscard]]auto view_id(const ParticleID id) const {
			// type: particle::internal::ScalarParticleRef<Read, ParticleField::none, ParticleData>
			return particle_container.template view_id<Read>(id);
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
		[[nodiscard]] std::vector<size_t> query_region(const core::Box & region) const {
			return particle_container.invoke_collect_indices_in_region(region);
		}

		[[nodiscard]] std::vector<size_t> query_region(const Domain & region) const {
			return query_region(core::Box(region.min_corner().value(), region.max_corner().value()));
		}


		// --------------
		// FUNCTIONAL OPS
		// --------------
		template<
			ParallelPolicy P=ParallelPolicy::Serial,
			VectorPolicy V=VectorPolicy::Auto,
			exec::IsKernel Kernel>
		void for_each_particle(Kernel && func, ParticleState state = ParticleState::ALL) {
			particle_container.template for_each_particle<P, V, Kernel>(std::forward<Kernel>(func), state);
		}

		template<
			ParallelPolicy P=ParallelPolicy::Serial,
			VectorPolicy V=VectorPolicy::Auto,
			exec::IsKernel Kernel>
		void for_each_particle_view(Kernel && func, ParticleState state = ParticleState::ALL) const {
			static_assert(std::remove_cvref_t<Kernel>::Write == ParticleField::none,
				"[APRIL] System: Kernel in for_each_particle_view must not have "
				"mutable fields (Kernel::Write must equal ParticleField::none"
			);
			particle_container.template for_each_particle<P, V, Kernel>(std::forward<Kernel>(func), state);
		}

		// TODO implement and fix
		template<ParticleField M, typename Mapper, typename T, typename Reducer = std::plus<T>>
		[[nodiscard]] T reduce(
			T initial_value,
			Mapper&& map_func, // map particle -> T
			Reducer&& reduce_func = {},
			ParticleState state = ParticleState::ALIVE
		) const {
			return particle_container.invoke_iterate_state(initial_value, map_func, reduce_func, state);
		}



		template<typename Func>
		void for_each_interaction_batch(Func && func) { // func(batch, bcp)
			particle_container.invoke_for_each_interaction_batch(std::forward<Func>(func));
		}

		template<exec::IsKernel Kernel>
		void for_each_interaction_pair(Kernel && func) { // func(particle, particle, dist)
			auto update_batch = [&]<container::batching::IsBCP BCP>(const container::batching::IsBatch auto& batch, BCP && apply_bcp) {
				using K = std::remove_cvref_t<Kernel>;

				auto bridge = [&](auto&& p1, auto&& p2) {
					// If it doesn't accept 2 arguments, assume it requires the distance vector (3 arguments)
					constexpr bool takes_2_args = std::is_invocable_v<K, decltype(p1), decltype(p2)>;
					constexpr bool requires_r = !takes_2_args;

					constexpr bool position_requested = particle::internal::has_field_v<K::Read, ParticleField::position>;

					if constexpr (requires_r && position_requested) {
						auto diff = p2.position - p1.position;

						const auto r = [&] {
							if constexpr (std::is_same_v<std::remove_cvref_t<BCP>, container::batching::NoBatchBCP>) {
								return diff;
							} else {
								return apply_bcp(diff);
							}
						}();

						func(p1, p2, r);
					} else if constexpr (requires_r && !position_requested) {
						static_assert(false,
							"[APRIL] ERROR: Kernel requires distance vector 'r', but ParticleField::position is missing from the ReadMask.");
					} else {
						func(p1, p2); // Fallback to 2-argument invocation
					}
				};

				auto kernel = exec::internal::KernelWrapper<K::Read, K::Write, K::Mode, decltype(bridge)>{bridge};
				execute_batch_kernel<ParallelPolicy::Serial, VectorPolicy::Scalar>(batch, kernel);
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
		core::Box simulation_box; // TODO rename to domain
		exec::Executor thread_executor;
		BoundaryTable boundary_table;
		ForceTable force_table;
		ControllerStorage controllers;
		FieldStorage fields;
		Container particle_container;

		struct alignas(64) PaddedThreadBuffer { std::vector<size_t> buffer; };
		std::vector<PaddedThreadBuffer> thread_update_buffers;
		std::vector<size_t> particles_to_update_buffer;

		double time_ = 0;
		size_t step_ = 0;

		core::SystemContext<System> system_context;
		utility::internal::TriggerContextImpl<System> trig_context;


		// private constructor since System should only be creatable through build_system(...)
		System(
			const ContainerDecl& container_cfg,
			const container::internal::ContainerCreateInfo & container_info,
			const core::Box& domain_in,
			const std::vector<ParticleRec>& particles,
			const BoundaryTable& boundaries_in,
			const ForceTable& forces_in,
			const ControllerStorage& controllers_in,
			const FieldStorage& fields_in
		)
			: simulation_box(domain_in),
			  thread_executor(exec::Executor()), // init stub incase we want to pass args in the future
			  boundary_table(boundaries_in),
			  force_table(forces_in),
			  controllers(controllers_in),
			  fields(fields_in),
			  particle_container(Container(container_cfg, container_info, thread_executor)),
			  system_context(*this),
			  trig_context(*this)
		{
			particle_container.invoke_build(particles);
			controllers.for_each_item([&](auto& c) { c.dispatch_init(context()); });
			fields.for_each_item([&](auto& f) { f.dispatch_init(context()); });
		}

		// Friend factory: only entry point for constructing a System
		// Avoids exposing constructor internals publicly.
		template <class Container, core::IsEnvironment EnvT>
		requires container::IsContainerDecl<Container, typename EnvT::traits>
		friend System<Container, typename EnvT::traits>
		build_system(
			 const EnvT & environment,
			 const Container& container,
			 BuildInfo * build_info
		);

		template<ParallelPolicy P, VectorPolicy V, container::batching::IsBatch Batch, exec::IsKernel Kernel>
		void execute_batch_kernel(const Batch& batch, Kernel&& kernel);
	};


	namespace core {
		namespace internal {
			template<typename>
			inline constexpr bool is_system_v = false;

			// template specialization: mark all System<C, Env> instantiations as true
			template<class C, IsEnvironmentTraits Traits>
			inline constexpr bool is_system_v<System<C, Traits>> = true;
		}

		template<typename T>
		concept IsSystem = internal::is_system_v<std::remove_cvref_t<T>>;
	}

} // namespace april::core

#include "april/core/internal/system_impl.hpp"












