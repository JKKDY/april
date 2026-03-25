#pragma once

#include <bit>
#include <vector>

#include "april/exec/policy.hpp"
#include "april/exec/particle_kernel.hpp"
#include "april/exec/config.hpp"

#include "april/math/range.hpp"

#include "april/interactions/interaction_table.hpp"
#include "april/core/domain.hpp"
#include "april/core/internal/environment_traits.hpp"

#include "april/particle/scalar_access.hpp"
#include "april/particle/packed_access.hpp"



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
			return particle::internal::PackedParticleRef<Read, ParticleField::none, ParticleAttributes> { // Assuming you kept the View alias
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
			AP_ASSERT(start <= self.capacity(), "Start index out of bounds: " + std::to_string(start));
			AP_ASSERT(stop <= self.capacity(), "Stop index out of bounds: " + std::to_string(stop));
			AP_ASSERT(start <= stop, "Invalid range: start > stop");

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
			AP_ASSERT(false, "add_particle not supported yet");
			self.add_particle(record);
		}
		void invoke_remove_particle(this auto&& self, const ParticleID id) {
			// TODO implement remove particle
			AP_ASSERT(false, "remove_particle not supported yet");
			self.remove_particle(id);
		}
		void invoke_resize_domain(this auto&& self, const core::Box & new_domain) {
			// TODO implement resize domain
			AP_ASSERT(false, "resize_domain not supported yet");
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
		exec::ExecutorRef<ThreadExecutor> thread_executor;

		template<exec::IsKernel Kernel>
		static auto adapt_indexed_kernel(Kernel && kernel) {
			return exec::make_kernel_wrapper<Kernel>(
				[kernel = std::forward<Kernel>(kernel)]<bool is_packed>(size_t i, auto && p) {
					if constexpr (requires { kernel(i, p); }) {
						return kernel(i, p); // user wants index
					} else if constexpr (requires { kernel(p); }) {
						return kernel(p); // user only wants particle
					} else {
						// TODO in C++26 use std::format and introspection to print out received signature
						// TODO print kernel name by implementing a name demangler
						static_assert(false, "[APRIL] Kernel is malformed! it must either have signature (size_t, auto && p) or (auto && p)");
					}
				}
			);
		}

		template<exec::IsKernel Kernel>
		static auto adapt_buffered_kernel(Kernel && kernel) {
			return exec::make_kernel_wrapper<Kernel>(
				[kernel = std::forward<Kernel>(kernel)]<bool is_packed>(size_t i, auto && p) {
					if constexpr (particle::IsPackedParticleRef<decltype(p)>) {
						static_assert(particle::IsPackedParticleAccessor<decltype(p)>);
						auto buffer = p.load_buffer();
						auto view = buffer.to_view();
						kernel(i, view);
						buffer.update_into(p);
					} else {
						static_assert(particle::IsScalarParticleAccessor<decltype(p)>);
						kernel(i, p);
					}
				}
			);
		}

		template<exec::IsKernel Kernel>
		static auto adapt_iterator_kernel(Kernel && kernel) {
			auto indexed_kernel = adapt_indexed_kernel(std::forward<Kernel>(kernel));
			auto buffered_kernel = adapt_buffered_kernel(std::move(indexed_kernel));
			return buffered_kernel;
		}



		template<ParallelPolicy P, VectorPolicy V, bool is_const, exec::IsKernel Kernel>
		void invoke_iterate_range(this auto&& self, Kernel && func, size_t start, size_t end) {
			constexpr auto mode = exec::internal::valid_execution_modes<V, std::remove_cvref_t<Kernel>::Mode>();
			auto kernel = adapt_iterator_kernel(func);
			self.template iterate_range<P, mode, is_const>(kernel, start, end);
		}

		template<ParallelPolicy P, VectorPolicy V, bool is_const, exec::IsKernel Kernel>
		void invoke_iterate_state(this auto&& self, Kernel && func, const ParticleState state) {
			using K = std::remove_cvref_t<Kernel>;
			auto kernel = adapt_iterator_kernel(func);
			constexpr auto mode = exec::internal::valid_execution_modes<V, K::Mode>();

			// try optimized implementation else fallback to default. Default assumes valid data for the entire iteration range
			if constexpr (requires {self.template iterate<P, V, is_const>(kernel, state);}) {
				self.template iterate<P, mode, is_const>(kernel, state);
			} else {
				// note: iterate_range makes no checks so if it encounters memory that cannot be interpreted as
				// particle data or memory it is not allowed to access it can crash
				// meaning the following default implementation is only safe if the container can guarantee that
				// only valid memory will be accessed (The build in AoS, SoA, AoSoA layouts support this)
				auto state_filter = [&]<typename Part>(size_t i, Part && p) {
					if constexpr (particle::IsPackedParticleAccessor<Part>) {
						static_assert(particle::IsPackedParticleRef<Part>);

						// const auto mask = (p.state.load() & +state) != 0;
						// if (!any(mask)) return; // if no particle is in requested state, skip this execution
						auto temp_buf = p.load_buffer();
						const auto mask = (temp_buf.state & +state) != 0;

						kernel(i, p.mask_with(mask));
					} else {
						static_assert(particle::IsScalarParticleAccessor<decltype(p)>);
						if (self.index_is_valid(i) && static_cast<int>(p.state & state)) {
							kernel(i, p);
						}

					}
				};

				self.template iterate_range<P, mode, is_const>(
					exec::internal::KernelWrapper<K::Read | ParticleField::state, K::Write, mode, decltype(state_filter)>(state_filter),
					0, self.capacity());
			}
		}

		template<ParticleField F>
		auto invoke_get_field_ptr(this auto&& self, size_t i) {
			AP_ASSERT(i < self.capacity(), "Index lies outside of capacity: " + std::to_string(i));
			return self.template get_field_ptr<F>(i);
		}

		template<ParticleField F>
		auto invoke_get_field_ptr_id(this auto&& self, ParticleID id) {
			AP_ASSERT(self.contains_id(id), "Got invalid Id: " + std::to_string(id));
			return self.template get_field_ptr_id<F>(id);
		}

		//------------------------
		// PARTICLE DATA ACCESSORS
		//------------------------
		template<ParticleField Read, ParticleField Write>
		[[nodiscard]] auto access_particle(this auto&& self, const size_t i) {

			constexpr bool is_const = std::is_const_v<std::remove_reference_t<decltype(self)>>;

			static_assert(!(is_const && Write != ParticleField::none),
				"APRIL ERROR: Cannot request write permissions (WriteMask != none) on a const Container. "
				"Either drop the write mask or ensure the container is mutable.");

			particle::internal::ParticleSource<Read, Write, ParticleAttributes> src;
			constexpr auto Mask = Read | Write;

			if constexpr (particle::internal::has_field_v<Mask, ParticleField::force>)
				src.force = self.template invoke_get_field_ptr<ParticleField::force>(i);
			if constexpr (particle::internal::has_field_v<Mask, ParticleField::position>)
				src.position = self.template invoke_get_field_ptr<ParticleField::position>(i);
			if constexpr (particle::internal::has_field_v<Mask, ParticleField::velocity>)
				src.velocity = self.template invoke_get_field_ptr<ParticleField::velocity>(i);
			if constexpr (particle::internal::has_field_v<Mask, ParticleField::old_position>)
				src.old_position = self.template invoke_get_field_ptr<ParticleField::old_position>(i);
			if constexpr (particle::internal::has_field_v<Mask, ParticleField::mass>)
				src.mass = self.template invoke_get_field_ptr<ParticleField::mass>(i);
			if constexpr (particle::internal::has_field_v<Mask, ParticleField::state>)
				src.state = self.template invoke_get_field_ptr<ParticleField::state>(i);
			if constexpr (particle::internal::has_field_v<Mask, ParticleField::type>)
				src.type = self.template invoke_get_field_ptr<ParticleField::type>(i);
			if constexpr (particle::internal::has_field_v<Mask, ParticleField::id>)
				src.id = self.template invoke_get_field_ptr<ParticleField::id>(i);
			if constexpr (particle::internal::has_field_v<Mask, ParticleField::attributes>)
				src.attributes = self.template invoke_get_field_ptr<ParticleField::attributes>(i);

			return src;
		}

		template<ParticleField Read, ParticleField Write>
		[[nodiscard]] auto access_particle_id(this auto&& self, const ParticleID id) {
			// it's optional to implement get_field_ptr_id. The fallback is to use id -> index and access_particle

			constexpr auto Mask = Read | Write;
			if constexpr (Mask == ParticleField::none) { // guard against none because 1 << std::countr_zero would produce UB
				return particle::internal::ParticleSource<Read, Write, ParticleAttributes>{};
			}

			// We pick the first active field in the Mask to test if 'get_field_ptr_id' exists.
			// We cannot check the function "in general" because it is a template.
			[[maybe_unused]] constexpr auto TestMask = static_cast<ParticleField>(1 << std::countr_zero(+(Read|Write)));

		    // does get_field_ptr_id<TestF>(id) compile?
		    if constexpr (requires { self.template get_field_ptr_id<TestMask>(id); }) {

		        // specialized path (direct ID access)
				constexpr bool is_const = std::is_const_v<std::remove_reference_t<decltype(self)>>;

		    	static_assert(!(is_const && Write != ParticleField::none),
		    		"APRIL ERROR: Cannot request write permissions (WriteMask != none) on a const Container. "
					"Either drop the write mask or ensure the container is mutable.");

		    	particle::internal::ParticleSource<Read, Write, ParticleAttributes> src;

		        if constexpr (particle::internal::has_field_v<Mask, ParticleField::force>)
        			src.force = self.template invoke_get_field_ptr_id<ParticleField::force>(id);
		        if constexpr (particle::internal::has_field_v<Mask, ParticleField::position>)
        			src.position = self.template invoke_get_field_ptr_id<ParticleField::position>(id);
		        if constexpr (particle::internal::has_field_v<Mask, ParticleField::velocity>)
        			src.velocity = self.template invoke_get_field_ptr_id<ParticleField::velocity>(id);
		        if constexpr (particle::internal::has_field_v<Mask, ParticleField::old_position>)
        			src.old_position = self.template invoke_get_field_ptr_id<ParticleField::old_position>(id);
		        if constexpr (particle::internal::has_field_v<Mask, ParticleField::mass>)
        			src.mass = self.template invoke_get_field_ptr_id<ParticleField::mass>(id);
		        if constexpr (particle::internal::has_field_v<Mask, ParticleField::state>)
        			src.state = self.template invoke_get_field_ptr_id<ParticleField::state>(id);
		        if constexpr (particle::internal::has_field_v<Mask, ParticleField::type>)
        			src.type = self.template invoke_get_field_ptr_id<ParticleField::type>(id);
		        if constexpr (particle::internal::has_field_v<Mask, ParticleField::id>)
        			src.id = self.template invoke_get_field_ptr_id<ParticleField::id>(id);
		        if constexpr (particle::internal::has_field_v<Mask, ParticleField::attributes>)
        			src.attributes = self.template invoke_get_field_ptr_id<ParticleField::attributes>(id);

		        return src;

		    } else {

		        // fallback path (ID -> Index -> Access)
		        return self.template access_particle<Read, Write>(self.invoke_id_to_index(id));
		    }
		}
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
		// 1. Must expose the standard traits (which your base Container class now does)
		requires {
			typename C::Config;
			typename C::ExecutionConfig;
			typename C::ParticleAttributes;
		} &&
		// 2. Must physically implement the operations
		HasContainerOps<C>
		&& std::derived_from<C, Container< ContainerBuildConfig<typename C::Config, typename C::ExecutionConfig, typename C::ParticleAttributes> >>;


	template<typename ContainerDecl, typename Traits, typename ExecCfg>
	concept IsContainerDecl =
		core::internal::IsEnvironmentTraits<Traits> &&
		exec::IsExecutionConfig<ExecCfg> &&
		requires {
			// 1. Synthesize the Build Configuration
			typename ContainerBuildConfig<ContainerDecl, ExecCfg, typename Traits::particle_attributes_t>;

			// 2. Check if the user's tag has an `impl` alias that accepts it
			typename ContainerDecl::template impl<
				ContainerBuildConfig<ContainerDecl, ExecCfg, typename Traits::particle_attributes_t>
			>;
		} &&
		// 3. Check if the resulting type satisfies the physical container requirements
		IsContainer<
			typename ContainerDecl::template impl<
				ContainerBuildConfig<ContainerDecl, ExecCfg, typename Traits::particle_attributes_t>
			>
		>;
} // namespace april::container


