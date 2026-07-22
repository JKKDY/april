#pragma once

#include <vector>

#include "april/exec/policy.hpp"
#include "april/exec/kernel.hpp"
#include "april/exec/threading/executor_config.hpp"
#include "april/exec/threading/executor_reference.hpp"
#include "april/exec/config.hpp"

#include "april/math/range.hpp"

#include "april/interactions/interaction_table.hpp"

#include "april/core/domain.hpp"
#include "april/core/internal/environment_traits.hpp"

#include "april/particle/access/scalar_access.hpp"
#include "april/particle/access/packed_access.hpp"



namespace april::container {
	struct ContainerFlags {
		bool periodic_x;			// domain is periodic along x-axis
		bool periodic_y;			// domain is periodic along y-axis
		bool periodic_z;			// domain is periodic along z-axis
		bool infinite_domain;		// particles outside of domain still interact normally (time complexity may go to O(n^2)
		bool particle_addable;		// particles can be added during run time
		bool particle_deletable;	// particles can be deleted during run time
	};

	struct ContainerHints {
		// TODO add regions that will be queried in the future so the container can keep track of particles better
		std::vector<ParticleID> interacting_particles;
		std::vector<core::Box> query_regions;
	};

	template<class ContainerCfg, class ExecutionCfg, particle::IsParticleAttributes Attributes>
	struct ContainerBuildConfig {
		using ContainerConfig = ContainerCfg;
		using ExecutionConfig = ExecutionCfg;
		using ParticleAttributes = Attributes;
		using ThreadExecutor = ExecutionConfig::ThreadExecutor;
		using Container = ContainerCfg::template impl<ContainerBuildConfig>;

		ExecutionConfig exec;
		ContainerConfig config;
		ContainerFlags flags {};
		ContainerHints hints {};
		interactions::internal::InteractionMap interaction_map {};
		core::Box domain {};
	};



	namespace internal {
		// check if T is a ContainerBuildConfig
		template<typename T>
		struct is_container_build_context : std::false_type {};

		template<typename Cfg, typename Exec, typename Attributes>
		struct is_container_build_context<
			ContainerBuildConfig<Cfg, Exec, Attributes>
		> : std::true_type {};
	}

	template<typename T>
	concept IsContainerBuildConfig =
		internal::is_container_build_context<std::remove_cvref_t<T>>::value;



	/**
	 * @brief CRTP-style base class for APRIL containers.
	 *
	 * Container implementations derive from this class to provide particle storage,
	 * traversal, ID lookup, and optional spatial acceleration structures. The base
	 * class supplies common particle accessors, kernel adaptation, fallback
	 * dispatch, and integration with APRIL's execution policies.
	 *
	 * Most simulation users interact with concrete containers such as DirectSum or
	 * LinkedCells rather than this base class directly.
	 */
	template<IsContainerBuildConfig BuildConfiguration>
	class Container {
	public:
		using ParticleAttributes = BuildConfiguration::ParticleAttributes;
		using ParticleRecord = particle::ParticleRecord<ParticleAttributes>;
		using Config = BuildConfiguration::ContainerConfig;
		using ExecutionConfig = BuildConfiguration::ExecutionConfig;
		using ThreadExecutor = BuildConfiguration::ThreadExecutor;

		static constexpr auto parallel_policy = ExecutionConfig::parallel_policy;
		static constexpr auto vector_policy = ExecutionConfig::vector_policy;

		explicit Container(const BuildConfiguration & context):
			config(context.config),
			flags(context.flags),
			hints(context.hints),
			interaction_map(context.interaction_map),
			domain(context.domain)
		{}

		void invoke_build(this auto&& self, const std::vector<ParticleRecord>& particles) {
			self.build(particles);
		}


		// ------------------
		// PARTICLE ACCESSORS
		// ------------------
		// INDEX ACCESSORS
		template<ParticleField Read, ParticleField Write>
		[[nodiscard]] auto at(this auto&& self, size_t index) {
			return particle::internal::ScalarParticleRef<Read, Write, ParticleAttributes> {
				self.template access_particle<Read, Write>(index)
			};
		}

		template<ParticleField Read>
		[[nodiscard]] auto view(this const auto& self, size_t index) {
			return particle::internal::ScalarParticleRef<Read, ParticleField::none, ParticleAttributes> {
				self.template access_particle<Read, ParticleField::none>(index)
			};
		}

		template<ParticleField Read, ParticleField Write>
		[[nodiscard]] auto at_packed(this auto&& self, size_t index) {
			return particle::internal::PackedParticleRef<Read, Write, ParticleAttributes> {
				self.template access_particle<Read, Write>(index)
			};
		}

		template<ParticleField Read>
		[[nodiscard]] auto view_packed(this const auto& self, size_t index) {
			return particle::internal::PackedParticleRef<Read, ParticleField::none, ParticleAttributes> {
				self.template access_particle<Read, ParticleField::none>(index)
			};
		}


		// ID ACCESSORS
		template<ParticleField Read, ParticleField Write>
		[[nodiscard]] auto at_id(this auto&& self, ParticleID id) {
			return particle::internal::ScalarParticleRef<Read, Write, ParticleAttributes> {
				self.template access_particle_id<Read, Write>(id)
			};
		}

		template<ParticleField Read>
		[[nodiscard]] auto view_id(this const auto & self, ParticleID id) {
			return particle::internal::ScalarParticleRef<Read, ParticleField::none, ParticleAttributes> {
				self.template access_particle_id<Read, ParticleField::none>(id)
			};
		}


		// ----------------
		// PREFETCHING
		// ----------------
		template <ParticleField Mask>
	    APRIL_FORCE_INLINE void prefetch_particle(this const auto& self, auto... args);

	    template <ParticleField Mask>
	    APRIL_FORCE_INLINE void prefetch_particle_nta(this const auto& self, auto... args);


		// ------------------
		// PARTICLE ITERATION
		// ------------------
		// filter by state (safe, performs checks to skip garbage data)
		template<
			ParallelPolicy P = ParallelPolicy::Serial,
			VectorPolicy V = vector_policy,
			exec::IsKernel Kernel>
		void for_each_particle(this auto&& self, Kernel && func, ParticleState state = ParticleState::ALL) {
			self.template invoke_iterate_state<P, V, false>(func, state);
		} // TODO add shortcuiting (if kernel returns a bool, stop when a true is encountered)

		// direct range based access (fast & branchless but unsafe; will not perform any checks)
		template<
			ParallelPolicy P = ParallelPolicy::Serial,
			VectorPolicy V = vector_policy,
			exec::IsKernel Kernel>
		void for_each_particle(this auto&& self, size_t start, size_t stop, Kernel && func) {
			APRIL_ASSERT(start <= self.capacity(), "Start index out of bounds: " + std::to_string(start));
			APRIL_ASSERT(stop <= self.capacity(), "Stop index out of bounds: " + std::to_string(stop));
			APRIL_ASSERT(start <= stop, "Invalid range: start > stop");

			self.template invoke_iterate_range<P, V, false>(func, start, stop);
		}


		template<ParticleField M, typename T, typename Mapper, typename Reducer = std::plus<T>>
		[[nodiscard]] T invoke_reduce( // TODO restrict callable Mapper, Reducer (invoke_reduce)
			this const auto& self,
			T initial_value,
			Mapper&& map_func,
			Reducer&& reduce_func = {},
			ParticleState state = ParticleState::ALIVE
		) {
			if constexpr (requires { self.reduce(initial_value, map_func, reduce_func, state); }) {
				// custom/optimized reducer
				return self.reduce(initial_value, std::forward<Mapper>(map_func), std::forward<Reducer>(reduce_func), state);
			} else {
				// default implementation using iterate function
				T curr = initial_value;
				auto kernel = [&](size_t, auto&& p) {
					T val = map_func(p);
					curr = reduce_func(curr, val);
				};

				self.template invoke_iterate<M, ParallelPolicy::Serial, vector_policy, true>(kernel, state);
				return curr;
			}
		}


		// --------
		// INDEXING
		// --------
		[[nodiscard]] size_t invoke_id_to_index(this const auto& self, ParticleID id) {
			return self.id_to_index(id);
		}

		[[nodiscard]] ParticleID invoke_min_id(this const auto& self) {
			return self.min_id();
		}

		[[nodiscard]] ParticleID invoke_max_id(this const auto& self) {
			return self.max_id();
		}
		[[nodiscard]] std::vector<math::Range> iteration_ranges() const {
			return {}; // TODO implement safe_iteration_ranges
		}


		// ---------------
		// BATCH ITERATION
		// ---------------
		template<ParallelPolicy P, typename Func>
		void invoke_for_each_interaction_batch(this auto&& self, Func && func) {
			self.template for_each_interaction_batch<P>(std::forward<Func>(func));
		}

		template<ParallelPolicy P, typename Func>
		void invoke_for_each_topology_batch(this auto&& self, Func&& func) {
			self.template for_each_topology_batch<P>(std::forward<Func>(func));
		}


		// -----------------
		// STRUCTURE UPDATES
		// -----------------
		// update container internals by scanning for all particle movements
		void invoke_rebuild_structure(this auto&& self) {
			self.rebuild_structure();
		}

		// perform partial container update given a list of particle indices
		void invoke_notify_moved(this auto&& self, const std::vector<size_t>& indices) {
			if constexpr (requires { self.notify_moved(indices); }) {
				self.notify_moved(indices);
			} else {
				// fallback: rebuild entire structure
				// useful in cases where the container cant do partial updates
				self.invoke_rebuild_structure();
			}
		}


		// ---------
		// MODIFIERS
		// ---------
		void invoke_add_particle(this auto&& self, const ParticleRecord & record) {
			// TODO implement add particle
			APRIL_ASSERT(false, "add_particle not supported yet");
			self.add_particle(record);
		}
		void invoke_remove_particle(this auto&& self, const ParticleID id) {
			// TODO implement remove particle
			APRIL_ASSERT(false, "remove_particle not supported yet");
			self.remove_particle(id);
		}
		void invoke_resize_domain(this auto&& self, const core::Box & new_domain) {
			// TODO implement resize domain
			APRIL_ASSERT(false, "resize_domain not supported yet");
			self.resize_domain(new_domain);
		}


		// -------
		// QUERIES
		// -------
		[[nodiscard]] bool invoke_index_is_valid(this const auto& self, size_t index) {
			return self.index_is_valid(index);
		}
		[[nodiscard]] bool invoke_contains_id(this const auto& self, ParticleID id) {
			return self.contains_id(id);
		}
		[[nodiscard]] size_t invoke_capacity(this const auto& self) {
			return self.capacity();
		}
		[[nodiscard]] size_t invoke_particle_count(this const auto& self) {
			return self.particle_count();
		}
		[[nodiscard]] std::vector<size_t> invoke_collect_indices_in_region(this const auto& self, const core::Box & region) {
			// TODO add tests for collect_indices_in_region
			if constexpr (requires { self.collect_indices_in_region(region); }) {
				return self.collect_indices_in_region(region);
			} else {
				std::vector<size_t> buffer;
				self.collect_indices_in_region(region, buffer);
				return buffer;
			}
		}
		[[nodiscard]] auto simulation_domain() const noexcept {
			return domain;
		}

		// TODO implement on containers, use this instead of non buffer version in system
		void invoke_collect_indices_in_region(this const auto& self, const core::Box & region, std::vector<size_t> & buffer) {
			return self.collect_indices_in_region(region, buffer);
		}

		void bind_executor(ThreadExecutor* raw_executor_ptr) {
			thread_executor.bind(raw_executor_ptr);
		}

	protected:
		const Config config;
		const ContainerFlags flags;
		const ContainerHints hints;
		const interactions::internal::InteractionMap interaction_map;
		const core::Box domain; // Note: in the future this may be adjustable during run time
		exec::ThreadExecutorRef<ThreadExecutor> thread_executor;


		template<ParallelPolicy P, VectorPolicy V, bool is_const, exec::IsKernel Kernel>
		void invoke_iterate_range(this auto&& self, Kernel && func, size_t start, size_t end);

		template<ParallelPolicy P, VectorPolicy V, bool is_const, exec::IsKernel Kernel>
		void invoke_iterate_state(this auto&& self, Kernel && func, ParticleState state);

		template<ParticleField F>
		auto invoke_get_field_ptr(this auto&& self, auto ... args);

		template<ParticleField F>
		auto invoke_get_field_ptr_id(this auto&& self, ParticleID id);

		//------------------------
		// PARTICLE DATA ACCESSORS
		//------------------------
		template<ParticleField Read, ParticleField Write>
		[[nodiscard]] auto access_particle(this auto&& self, const auto ... args);

		template<ParticleField Read, ParticleField Write>
		[[nodiscard]] auto access_particle_id(this auto&& self, const ParticleID id);
	};



	template<typename C>
	concept HasContainerOps = requires (
	    C c,
	    const C cc,
	    ParticleID id,
	    size_t index,
	    const core::Box& region,
	    const particle::ParticleRecord<typename C::ParticleAttributes>& p,
	    const std::vector<particle::ParticleRecord<typename C::ParticleAttributes>>& particles
	) {
		// minimal implemented interface (except get_field_ptr since that is not part of the public API)
	    { c.build(particles) };
	    { c.rebuild_structure() };

		{ cc.capacity() } -> std::convertible_to<size_t>;
	    { cc.particle_count() } -> std::convertible_to<size_t>;
		{ cc.min_id() } -> std::convertible_to<ParticleID>;
	    { cc.max_id() } -> std::convertible_to<ParticleID>;

	    { cc.id_to_index(id) } -> std::convertible_to<size_t>;
		{ cc.contains_id(id) } -> std::convertible_to<bool>;
		{ cc.index_is_valid(id) } -> std::convertible_to<bool>;

	    { c.collect_indices_in_region(region) } -> std::convertible_to<std::vector<size_t>>;

	    { c.template for_each_interaction_batch<ParallelPolicy::Serial>([](auto&&){}) };
		{ c.template for_each_topology_batch<ParallelPolicy::Serial>([](auto&&){}) };
	};


	template<typename C>
	concept IsContainer =
		requires {
			typename C::Config;
			typename C::ExecutionConfig;
			typename C::ParticleAttributes;
		}
		&& HasContainerOps<C> // has implemented minimal interface
		&& std::derived_from<C, // is derived from container base class
			Container<ContainerBuildConfig<typename C::Config, typename C::ExecutionConfig, typename C::ParticleAttributes>>
		>;


	template<typename ContainerDecl, typename Traits, typename ExecCfg>
	concept IsContainerDecl =
		core::internal::IsEnvironmentTraits<Traits> &&
		exec::IsExecutionConfig<ExecCfg> &&
		requires { // check if ContainerDecl has an impl member that is template-able on a build configuration
			typename ContainerBuildConfig<ContainerDecl, ExecCfg, typename Traits::particle_attributes_t>;

			typename ContainerDecl::template impl<
				ContainerBuildConfig<ContainerDecl, ExecCfg, typename Traits::particle_attributes_t>
			>;
		} &&
		IsContainer< // and has represents a valid container
			typename ContainerDecl::template impl<
				ContainerBuildConfig<ContainerDecl, ExecCfg, typename Traits::particle_attributes_t>
			>
		>;
} // namespace april::container

#include "april/containers/internal/container_impl.hpp"
