#pragma once

#include "april/core/environment.hpp"
#include "april/core/domain.hpp"
#include "april/containers/container.hpp"
#include "april/containers/batching/common.hpp"
#include "april/exec/policy.hpp"
#include "april/exec/particle_kernel.hpp"
#include "april/core/context.hpp"
#include "april/exec/config.hpp"

namespace april {
	struct BuildInfo;

	template <class SystemConfig>
	class System;


	template<
		class ParticleContainer,
		core::internal::IsEnvironmentTraits EnvTraits,
		exec::IsExecutionConfig ExecCfg>
	struct SystemConfig {
		using Container = ParticleContainer;
		using ExecutionConfig = ExecCfg;
		using BoundaryTable = EnvTraits::boundary_table_t;
		using InteractionTable = EnvTraits::force_table_t;
		using ParticleRecord = EnvTraits::particle_record_t;
		using ParticleAttributes = EnvTraits::particle_attributes_t;
		using Controllers = EnvTraits::controller_storage_t;
		using Fields = EnvTraits::field_storage_t;

		Container container; // Keep as value
		ExecutionConfig execution_config; // Keep as value
		std::vector<ParticleRecord> particles; // Change to value!
		BoundaryTable boundaries; // Change to value!
		InteractionTable interactions; // Change to value!
		Controllers controllers; // Change to value!
		Fields fields; // Change to value!
	};


	template <class Container, core::IsEnvironment Env, exec::IsExecutionConfig ExecCfg>
	requires container::IsContainerDecl<Container, typename Env::traits>
	auto build_system(
		const Env & environment,
		const Container & container_config,
		const ExecCfg & execution_config,
		BuildInfo * build_info
	);




	template <class SystemConfig>
	class System final {
		// --------------
		// INTERNAL TYPES
		// --------------
		using Controllers		= SystemConfig::Controllers;
		using Fields			= SystemConfig::Fields;
		using InteractionTable 	= SystemConfig::InteractionTable;
		using BoundaryTable  	= SystemConfig::BoundaryTable;
	public:
		// ----------------
		// PUBLIC API TYPES
		// ----------------
		using SysContext = core::SystemContext<System>;
		using TrigContext = utility::internal::TriggerContextImpl<System>;
		using ParticleAttributes= SystemConfig::ParticleAttributes;
		using ParticleRec		= SystemConfig::ParticleRecord;
		using Container			= SystemConfig::Container;
		using ExecutionConfig   = SystemConfig::ExecutionConfig;

		// convinience aliases
		static constexpr auto parallel_policy = ExecutionConfig::parallel_policy;
		static constexpr auto vector_policy = ExecutionConfig::vector_policy;

	private:
		// ---------------
		// PRIVATE MEMBERS
		// ---------------
		exec::Executor thread_executor;
		BoundaryTable boundary_table;
		InteractionTable force_table;
		Controllers controllers;
		Fields fields;
		Container particle_container;
		ExecutionConfig exec_config;

		struct alignas(64) PaddedThreadBuffer { std::vector<size_t> buffer; };
		std::vector<PaddedThreadBuffer> thread_update_buffers;
		std::vector<size_t> particles_to_update_buffer;

		double time_ = 0;
		size_t step_ = 0;

		core::SystemContext<System> system_context;
		utility::internal::TriggerContextImpl<System> trig_context;


		// private constructor since System should only be creatable through build_system(...)
		explicit System(SystemConfig&& config)
		  : thread_executor(typename ExecutionConfig::ThreadExecutor(config.execution_config.executer_config)),
			boundary_table(std::move(config.boundaries)),     // Zero-cost pointer swap!
			force_table(std::move(config.interactions)),      // Zero-cost pointer swap!
			controllers(std::move(config.controllers)),       // Zero-cost pointer swap!
			fields(std::move(config.fields)),                 // Zero-cost pointer swap!
			particle_container(std::move(config.container)),  // Zero-cost pointer swap!
			exec_config(std::move(config.execution_config)),
			system_context(*this),
			trig_context(*this)
		{
			particle_container.bind_executor(&thread_executor);
			particle_container.invoke_build(config.particles);
			controllers.for_each_item([&](auto& c) { c.dispatch_init(context()); });
			fields.for_each_item([&](auto& f) { f.dispatch_init(context()); });
		}

		// Friend factory: only entry point for constructing a System
		// Avoids exposing constructor internals publicly.
		template <class C, core::IsEnvironment E, exec::IsExecutionConfig EC>
	    requires container::IsContainerDecl<C, typename E::traits>
	    friend auto build_system(
		    const E&, const C&, const EC&, BuildInfo*
	    );

		template<VectorPolicy V, container::batching::IsBatch Batch, exec::IsKernel Kernel>
		void execute_batch_kernel(const Batch& batch, Kernel&& kernel);
	public:
		// make non copyable
		System(const System&) = delete;
		System& operator=(const System&) = delete;
		System(System&&) = delete;
		System& operator=(System&&) = delete;

		// -----------------
		// LIFECYCLE & STATE
		// -----------------
		[[nodiscard]] double time() const noexcept { return time_; }
		[[nodiscard]] size_t step() const noexcept { return step_; }
		// TODO check if domain() is used anywhere and if it can be deleted
		[[nodiscard]] Domain domain() const { return { particle_container.simulation_domain().min,  particle_container.simulation_domain().extent}; }
		[[nodiscard]] core::Box box() const { return particle_container.simulation_domain(); }

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
			VectorPolicy V=vector_policy,
			exec::IsKernel Kernel>
		void for_each_particle(Kernel && func, ParticleState state = ParticleState::ALL) {
			particle_container.template for_each_particle<P, V, Kernel>(std::forward<Kernel>(func), state);
		}

		template<
			ParallelPolicy P=ParallelPolicy::Serial,
			VectorPolicy V=vector_policy,
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



		template<ParallelPolicy P, typename Func>
		void for_each_interaction_batch(Func && func) { // func signature: func(batch, bcp)
			particle_container.template invoke_for_each_interaction_batch<P>(std::forward<Func>(func));
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
				execute_batch_kernel<ParallelPolicy::Serial, vector_policy>(batch, kernel);
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
		// note: these methods will be eventually abstracted into compute stages which will be composable in a compute pipeline
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
	};


	namespace core {
		namespace internal {
			template<typename>
			inline constexpr bool is_system_v = false;

			// mark all System<C, Env, ExecCfg> instantiations as true (template specialization
			template<class Container, IsEnvironmentTraits EnvTraits, exec::IsExecutionConfig ExecCfg>
			inline constexpr bool is_system_v<System<SystemConfig<Container, EnvTraits, ExecCfg>>> = true;
		}

		template<typename T>
		concept IsSystem = internal::is_system_v<std::remove_cvref_t<T>>;
	}

} // namespace april::core

#include "april/core/internal/system_impl.hpp"

