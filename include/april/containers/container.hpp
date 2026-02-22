#pragma once

#include <bit>
#include <vector>

#include "april/exec/policy.hpp"
#include "april/exec/particle_kernel.hpp"

#include "april/math/range.hpp"

#include "april/forces/force_table.hpp"
#include "april/env/domain.hpp"
#include "april/env/traits.hpp"

#include "april/particle/scalar_access.hpp"
#include "april/particle/packed_access.hpp"





namespace april::container {

	namespace internal {
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
			std::vector<env::Box> query_regions;
		};

		struct ContainerCreateInfo {
			ContainerFlags flags {};
			ContainerHints hints {};
			force::internal::InteractionSchema force_schema {};
			env::Box domain {};
		};
	}



	template<class C, env::IsUserData U>
	class Container {
	public:
		using ParticleRecord = env::internal::ParticleRecord<U>;
		using UserData = U;
		using Config = C;

		Container(const Config & config, const internal::ContainerCreateInfo & info):
			config(config), flags(info.flags), hints(info.hints), force_schema(info.force_schema), domain(info.domain)
		{}

		void invoke_build(this auto&& self, const std::vector<ParticleRecord>& particles) {
			self.build(particles);
		}


		// ------------------
		// PARTICLE ACCESSORS
		// ------------------
		// INDEX ACCESSORS
		template<ParticleField M>
		[[nodiscard]] auto at(this auto&& self, size_t index) {
			return env::ScalarParticleRef<M, U>{ self.template access_particle<M>(index) };
		}

		template<ParticleField M>
		[[nodiscard]] auto view(this const auto& self, size_t index) {
			return env::ScalarParticleView<M, U>{ self.template access_particle<M>(index) };
		}

		template<ParticleField M>
		[[nodiscard]] auto restricted_at(this auto&& self, size_t index) {
			return env::ScalarRestrictedParticleRef<M, U>{ self.template access_particle<M>(index) };
		}

		template<ParticleField M>
		[[nodiscard]] auto at_packed(this auto&& self, size_t index) {
			return env::PackedParticleRef<M, U>{ self.template access_particle<M>(index) };
		}

		template<ParticleField M>
		[[nodiscard]] auto view_packed(this const auto& self, size_t index) {
			return env::PackedParticleView<M, U>{ self.template access_particle<M>(index) };
		}

		template<ParticleField M>
		[[nodiscard]] auto restricted_at_packed(this auto&& self, size_t index) {
			return env::PackedRestrictedParticleRef<M, U>{ self.template access_particle<M>(index) };
		}


		// ID ACCESSORS
		template<ParticleField M>
		[[nodiscard]] auto at_id(this auto&& self, ParticleID id) {
			return env::ScalarParticleRef<M, U>{ self.template access_particle_id<M>(id) };
		}

		template<ParticleField M>
		[[nodiscard]] auto view_id(this const auto & self, ParticleID id) {
			return env::ScalarParticleView<M, U>{ self.template access_particle_id<M>(id) };
		}

		template<ParticleField M>
		[[nodiscard]] auto restricted_at_id(this auto&& self, ParticleID id) {
			return env::ScalarRestrictedParticleRef<M, U>{ self.template access_particle_id<M>(id) };
		}



		// ------------------
		// PARTICLE ITERATION
		// ------------------
		// filter by state (safe, performs checks to skip garbage data)
		template<
			ParticleField M,
			ParallelPolicy P = ParallelPolicy::Serial,
			VectorPolicy V = VectorPolicy::Auto,
			exec::IsKernel Kernel>
		void for_each_particle(this auto&& self, Kernel && func, ParticleState state = ParticleState::ALL) {
			self.template invoke_iterate_state<M, P, V, false>(func, state);
		}

		// view variant (const)
		template<
			ParticleField M,
			ParallelPolicy P = ParallelPolicy::Serial,
			VectorPolicy V = VectorPolicy::Auto,
			exec::IsKernel Kernel>
		void for_each_particle_view(this const auto& self, Kernel && func, ParticleState state = ParticleState::ALL) {
			self.template invoke_iterate_state<M, P, V, true>(func, state);
		}

		// direct range based access (fast & branchless but unsafe; will not perform any checks)
		template<
			ParticleField M,
			ParallelPolicy P = ParallelPolicy::Serial,
			VectorPolicy V = VectorPolicy::Auto,
			exec::IsKernel Kernel>
		void for_each_particle(this auto&& self, size_t start, size_t stop, Kernel && func) {
			AP_ASSERT(start <= self.capacity(), "Start index out of bounds: " + std::to_string(start));
			AP_ASSERT(stop <= self.capacity(), "Stop index out of bounds: " + std::to_string(stop));
			AP_ASSERT(start <= stop, "Invalid range: start > stop");

			self.template invoke_iterate_range<M,  P, V, false>(func, start, stop);
		}

		// view variant (const)
		template<
			ParticleField M,
			ParallelPolicy P = ParallelPolicy::Serial,
			VectorPolicy V = VectorPolicy::Auto,
			exec::IsKernel Kernel>
		void for_each_particle_view(this const auto& self, size_t start, size_t stop, Kernel && func) {
			AP_ASSERT(start <= self.capacity(), "Start index out of bounds: " + std::to_string(start));
			AP_ASSERT(stop <= self.capacity(), "Stop index out of bounds: " + std::to_string(stop));
			AP_ASSERT(start <= stop, "Invalid range");

			self.template invoke_iterate_range<M,  P, V, true>(func, start, stop);
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

				self.template invoke_iterate<M, ParallelPolicy::Serial, VectorPolicy::Auto, true>(kernel, state);
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
		template<typename Func> // TODO restrict callable Func (invoke_for_each_batch)
		void invoke_for_each_interaction_batch(this auto&& self, Func && func) {
			self.for_each_interaction_batch(std::forward<Func>(func));
		}

		template<typename Func>
		void invoke_for_each_topology_batch(this auto&& self, Func&& func) {
			self.for_each_topology_batch(std::forward<Func>(func));
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
			AP_ASSERT(false, "add_particle not supported yet");
			self.add_particle(record);
		}
		void invoke_remove_particle(this auto&& self, const ParticleID id) {
			AP_ASSERT(false, "remove_particle not supported yet");
			self.remove_particle(id);
		}
		void invoke_resize_domain(this auto&& self, const env::Box & new_domain) {
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
		[[nodiscard]] std::vector<size_t> invoke_collect_indices_in_region(this const auto& self, const env::Box & region) {
			if constexpr (requires { self.collect_indices_in_region(region); }) {
				return self.collect_indices_in_region(region);
			} else {
				std::vector<size_t> buffer;
				self.collect_indices_in_region(region, buffer);
				return buffer;
			}
		}

		// TODO implement on containers, use this instead of non buffer version in system
		void invoke_collect_indices_in_region(this const auto& self, const env::Box & region, std::vector<size_t> & buffer) {
			return self.collect_indices_in_region(region, buffer);
		}


	protected:
		const Config config;
		const internal::ContainerFlags flags;
		const internal::ContainerHints hints;
		const force::internal::InteractionSchema force_schema;
		const env::Box domain; // Note: in the future this may be adjustable during run time

		template<ParticleField M, ParallelPolicy P, VectorPolicy V, bool is_const, exec::IsKernel Kernel>
		void invoke_iterate_range(this auto&& self, Kernel && func, size_t start, size_t end) {
			auto bridge = [&](size_t i, auto && p) {
				if constexpr (requires { func(i, p); }) {
					return func(i, p); // user wants index
				} else {
					return func(p); // user only wants particle
				}
			};

			// auto kernel = exec::internal::KernelWrapper<Kernel::Mode, decltype(bridge)>{bridge};
			auto kernel = exec::internal::KernelWrapper<std::remove_cvref_t<Kernel>::Mode, decltype(bridge)>{bridge};
			self.template iterate_range<M, P, V, is_const>(kernel, start, end);
		}

		template<ParticleField M, ParallelPolicy P, VectorPolicy V, bool is_const, exec::IsKernel Kernel>
		void invoke_iterate_state(this auto&& self, Kernel && func, const ParticleState state) {
			auto kernel = universal_kernel([&](size_t i, auto && p) {
				if constexpr (requires { func(i, p); }) {
					func(i, p); // user wants index
				} else {
					func(p); // user only wants particle
				}
			});

			// try optimized implementation. Else fallback to default. Default assumes valid data for the entire iteration range
			if constexpr (requires {self.template iterate<M, P, V, is_const>(kernel, state);}) {
				self.template iterate<M, P, V, is_const>(kernel, state);
			} else {
				// note: iterate_range makes no checks so if it encounters memory that cannot be interpreted as
				// particle data or memory it is not allowed to access it can crash.
				self.template iterate_range<M | ParticleField::state, P, V, is_const>(universal_kernel(
					[&](size_t i, auto && p) {
					if (self.index_is_valid(i) && static_cast<int>(p.state & state)) {
						kernel(i, p);
					}
				}), 0, self.capacity());
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
		template<ParticleField M>
		[[nodiscard]] auto access_particle(this auto&& self, const size_t i) {

			constexpr bool IsConst = std::is_const_v<std::remove_reference_t<decltype(self)>>;
			env::internal::ParticleSource<M, U, IsConst> src;

			if constexpr (env::has_field_v<M, ParticleField::force>)
				src.force = self.template invoke_get_field_ptr<ParticleField::force>(i);
			if constexpr (env::has_field_v<M, ParticleField::position>)
				src.position = self.template invoke_get_field_ptr<ParticleField::position>(i);
			if constexpr (env::has_field_v<M, ParticleField::velocity>)
				src.velocity = self.template invoke_get_field_ptr<ParticleField::velocity>(i);
			if constexpr (env::has_field_v<M, ParticleField::old_position>)
				src.old_position = self.template invoke_get_field_ptr<ParticleField::old_position>(i);
			if constexpr (env::has_field_v<M, ParticleField::mass>)
				src.mass = self.template invoke_get_field_ptr<ParticleField::mass>(i);
			if constexpr (env::has_field_v<M, ParticleField::state>)
				src.state = self.template invoke_get_field_ptr<ParticleField::state>(i);
			if constexpr (env::has_field_v<M, ParticleField::type>)
				src.type = self.template invoke_get_field_ptr<ParticleField::type>(i);
			if constexpr (env::has_field_v<M, ParticleField::id>)
				src.id = self.template invoke_get_field_ptr<ParticleField::id>(i);
			if constexpr (env::has_field_v<M, ParticleField::user_data>)
				src.user_data = self.template invoke_get_field_ptr<ParticleField::user_data>(i);

			return src;
		}

		template<ParticleField M>
		[[nodiscard]] auto access_particle_id(this auto&& self, const ParticleID id) {
			// its optional to implement get_field_ptr_id. The fallback is to use id -> index and access_particle

			// We pick the first active field in the Mask to test if 'get_field_ptr_id' exists.
			// We cannot check the function "in general" because it is a template.
			[[maybe_unused]] constexpr auto TestF = static_cast<ParticleField>(1 << std::countr_zero(+M));

		    // does 'get_field_ptr_id<TestF>(id)' compile?
		    if constexpr (requires { self.template get_field_ptr_id<TestF>(id); }) {

		        // specialized path (direct ID access)
		        constexpr bool IsConst = std::is_const_v<std::remove_reference_t<decltype(self)>>;
		        env::internal::ParticleSource<M, U, IsConst> src;

		        if constexpr (env::has_field_v<M, ParticleField::force>)
        			src.force = self.template invoke_get_field_ptr_id<ParticleField::force>(id);
		        if constexpr (env::has_field_v<M, ParticleField::position>)
        			src.position = self.template invoke_get_field_ptr_id<ParticleField::position>(id);
		        if constexpr (env::has_field_v<M, ParticleField::velocity>)
        			src.velocity = self.template invoke_get_field_ptr_id<ParticleField::velocity>(id);
		        if constexpr (env::has_field_v<M, ParticleField::old_position>)
        			src.old_position = self.template invoke_get_field_ptr_id<ParticleField::old_position>(id);
		        if constexpr (env::has_field_v<M, ParticleField::mass>)
        			src.mass = self.template invoke_get_field_ptr_id<ParticleField::mass>(id);
		        if constexpr (env::has_field_v<M, ParticleField::state>)
        			src.state = self.template invoke_get_field_ptr_id<ParticleField::state>(id);
		        if constexpr (env::has_field_v<M, ParticleField::type>)
        			src.type = self.template invoke_get_field_ptr_id<ParticleField::type>(id);
		        if constexpr (env::has_field_v<M, ParticleField::id>)
        			src.id = self.template invoke_get_field_ptr_id<ParticleField::id>(id);
		        if constexpr (env::has_field_v<M, ParticleField::user_data>)
        			src.user_data = self.template invoke_get_field_ptr_id<ParticleField::user_data>(id);

		        return src;

		    } else {

		        // fallback path (ID -> Index -> Access)
		        return self.template access_particle<M>(self.invoke_id_to_index(id));
		    }
		}
	};



	template<typename C>
	concept HasContainerOps = requires (
	    C c,
	    const C cc,
	    ParticleID id,
	    size_t index,
	    const env::Box& region,
	    const env::internal::ParticleRecord<typename C::UserData>& p,
	    const std::vector<env::internal::ParticleRecord<typename C::UserData>>& particles
	) {
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

	    { c.for_each_interaction_batch([](auto&&){}) };
		{ c.for_each_topology_batch([](auto&&){}) };

		// for future use
		// { c.add_particle(p) };
		// { c.remove_particle(id) };
	};


	template<typename C> concept IsContainer =
		// must define types (Config, UserData)
		requires {
			typename C::Config;
			typename C::UserData;
		} &&
		// Config must have impl typename pointing to Container type
		// container must only depend on user data as template argument
		requires {
			typename C::Config::template impl<typename C::UserData>;
			requires std::same_as<C, typename C::Config::template impl<typename C::UserData>>;
		} &&
		// Must inherit from the Container
		std::derived_from<C, Container<
			typename C::Config,
			typename C::UserData
		>> &&
		// Must implement the Structural Contract
		HasContainerOps<C>;


	template<typename ContainerDecl, typename Traits> concept IsContainerDecl =
		env::internal::IsEnvironmentTraits<Traits>
		&& requires { typename ContainerDecl::template impl<typename Traits::user_data_t>; }
		&& IsContainer<typename ContainerDecl::template impl<typename Traits::user_data_t>>;

} // namespace april::container





