#pragma once

#include <bit>
#include <vector>

#include "april/forces/force_table.hpp"
#include "april/particle/fields.hpp"
#include "april/particle/access.hpp"
#include "april/env/domain.hpp"
#include "april/env/traits.hpp"



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
			std::vector<env::ParticleID> interacting_particles;
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
		template<env::FieldMask M>
		[[nodiscard]] auto at(this auto&& self, size_t index) {
			return env::ParticleRef<M, U>{ self.template access_particle<M>(index) };
		}

		template<env::FieldMask M>
		[[nodiscard]] auto view(this const auto& self, size_t index) {
			return env::ParticleView<M, U>{ self.template access_particle<M>(index) };
		}

		template<env::FieldMask M>
		[[nodiscard]] auto restricted_at(this auto&& self, size_t index) {
			return env::RestrictedParticleRef<M, U>{ self.template access_particle<M>(index) };
		}

		// ID ACCESSORS
		template<env::FieldMask M>
		[[nodiscard]] auto at_id(this auto&& self, env::ParticleID id) {
			return env::ParticleRef<M, U>{ self.template access_particle_id<M>(id) };
		}

		template<env::FieldMask M>
		[[nodiscard]] auto view_id(this const auto & self, env::ParticleID id) {
			return env::ParticleView<M, U>{ self.template access_particle_id<M>(id) };
		}

		template<env::FieldMask M>
		[[nodiscard]] auto restricted_at_id(this auto&& self, env::ParticleID id) {
			return env::RestrictedParticleRef<M, U>{ self.template access_particle_id<M>(id) };
		}



		// ------------------
		// PARTICLE ITERATION
		// ------------------
		// filter by state (safe, performs checks to skip garbage data)
		template<env::FieldMask M, ExecutionPolicy Policy = ExecutionPolicy::Seq, typename Func>
		void for_each_particle(this auto&& self, Func && func, env::ParticleState state = env::ParticleState::ALL) {
			self.template invoke_iterate_state<M, Policy, false>(func, state);
		}
		template<env::FieldMask M, ExecutionPolicy Policy = ExecutionPolicy::Seq, typename Func>
		void for_each_particle_view(this const auto& self, Func && func, env::ParticleState state = env::ParticleState::ALL) {
			self.template invoke_iterate_state<M,Policy, true>(func, state);
		}

		// direct range based access (fast & branchless but unsafe will not perform any checks)
		template<env::FieldMask M, ExecutionPolicy Policy = ExecutionPolicy::Seq, typename Func>
		void for_each_particle(this auto&& self, size_t start, size_t stop, Func && func) {
			AP_ASSERT(start <= self.capacity(), "Start index out of bounds: " + std::to_string(start));
			AP_ASSERT(stop <= self.capacity(), "Stop index out of bounds: " + std::to_string(stop));
			AP_ASSERT(start <= stop, "Invalid range: start > stop");

			self.template invoke_iterate_range<M, Policy, false>(func, start, stop);
		}
		template<env::FieldMask M, ExecutionPolicy Policy = ExecutionPolicy::Seq, typename Func>
		void for_each_particle_view(this const auto& self, size_t start, size_t stop, Func && func) {
			AP_ASSERT(start <= self.capacity(), "Start index out of bounds: " + std::to_string(start));
			AP_ASSERT(stop <= self.capacity(), "Stop index out of bounds: " + std::to_string(stop));
			AP_ASSERT(start <= stop, "Invalid range");

			self.template invoke_iterate_range<M, Policy, true>(func, start, stop);
		}

		template<env::FieldMask M, typename T, typename Mapper, typename Reducer = std::plus<T>>
		[[nodiscard]] T invoke_reduce( // TODO restrict callable Mapper, Reducer (invoke_reduce)
			this const auto& self,
			T initial_value,
			Mapper&& map_func,
			Reducer&& reduce_func = {},
			env::ParticleState state = env::ParticleState::ALIVE
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

				self.template invoke_iterate<M, ExecutionPolicy::Seq, true>(kernel, state);
				return curr;
			}
		}


		// --------
		// INDEXING
		// --------
		[[nodiscard]] size_t invoke_id_to_index(this const auto& self, env::ParticleID id) {
			return self.id_to_index(id);
		}

		[[nodiscard]] env::ParticleID invoke_min_id(this const auto& self) {
			return self.min_id();
		}

		[[nodiscard]] env::ParticleID invoke_max_id(this const auto& self) {
			return self.max_id();
		}
		[[nodiscard]] std::vector<std::pair<size_t, size_t>> safe_iteration_ranges() const {
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
		void invoke_remove_particle(this auto&& self, const env::ParticleID id) {
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
		[[nodiscard]] bool invoke_contains_id(this const auto& self, env::ParticleID id) {
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

		template<env::FieldMask M, ExecutionPolicy Policy, bool is_const, typename Func>
		void invoke_iterate_range(this auto&& self, Func && func, size_t start, size_t end) {
			auto kernel = [&](size_t i, auto && p) {
				if constexpr (requires { func(i, p); }) {
					func(i, p); // user wants index
				} else {
					func(p); // user only wants particle
				}
			};

			self.template iterate_range<M, Policy, is_const>(kernel, start, end);
		}

		template<env::FieldMask M, ExecutionPolicy Policy, bool is_const, typename Func>
		void invoke_iterate_state(this auto&& self, Func && func, const env::ParticleState state) {
			auto kernel = [&](size_t i, auto && p) {
				if constexpr (requires { func(i, p); }) {
					func(i, p); // user wants index
				} else {
					func(p); // user only wants particle
				}
			};

			// try optimized implementation. Else fallback to default. Default assumes valid data for the entire iteration range
			if constexpr (requires {self.template iterate<M, Policy, is_const>(kernel, state);}) {
				self.template iterate<M, Policy, is_const>(kernel, state);
			} else {
				// note: iterate_range makes no checks so if it encounters memory that cannot be interpreted as
				// particle data or memory it is not allowed to access it can crash.
				self.template iterate_range<M | env::Field::state, Policy, is_const>(
					[&](size_t i, auto && p) {
					if (self.index_is_valid(i) && static_cast<int>(p.state & state)) {
						kernel(i, p);
					}
				}, 0, self.capacity());
			}
		}

		template<env::Field F>
		auto invoke_get_field_ptr(this auto&& self, size_t i) {
			return self.template get_field_ptr<F>(i);
		}

		template<env::Field F>
		auto invoke_get_field_ptr_id(this auto&& self, env::ParticleID id) {
			return self.template get_field_ptr_id<F>(id);
		}

		//------------------------
		// PARTICLE DATA ACCESSORS
		//------------------------
		template<env::FieldMask M>
		[[nodiscard]] auto access_particle(this auto&& self, const size_t i) {
			// Safety check: If Mask is empty, return empty source immediately.

			if constexpr (M == 0) return env::ParticleSource<0, U, false>{};

			constexpr bool IsConst = std::is_const_v<std::remove_reference_t<decltype(self)>>;
			env::ParticleSource<M, U, IsConst> src;

			if constexpr (env::has_field_v<M, env::Field::force>)
				src.force = self.template invoke_get_field_ptr<env::Field::force>(i);
			if constexpr (env::has_field_v<M, env::Field::position>)
				src.position = self.template invoke_get_field_ptr<env::Field::position>(i);
			if constexpr (env::has_field_v<M, env::Field::velocity>)
				src.velocity = self.template invoke_get_field_ptr<env::Field::velocity>(i);
			if constexpr (env::has_field_v<M, env::Field::old_position>)
				src.old_position = self.template invoke_get_field_ptr<env::Field::old_position>(i);
			if constexpr (env::has_field_v<M, env::Field::mass>)
				src.mass = self.template invoke_get_field_ptr<env::Field::mass>(i);
			if constexpr (env::has_field_v<M, env::Field::state>)
				src.state = self.template invoke_get_field_ptr<env::Field::state>(i);
			if constexpr (env::has_field_v<M, env::Field::type>)
				src.type = self.template invoke_get_field_ptr<env::Field::type>(i);
			if constexpr (env::has_field_v<M, env::Field::id>)
				src.id = self.template invoke_get_field_ptr<env::Field::id>(i);
			if constexpr (env::has_field_v<M, env::Field::user_data>)
				src.user_data = self.template invoke_get_field_ptr<env::Field::user_data>(i);

			return src;
		}

		template<env::FieldMask M>
		[[nodiscard]] auto access_particle_id(this auto&& self, const env::ParticleID id) {
			// its optional to implement get_field_ptr_id. The fallback is to use id -> index and access_particle

			// Safety check: If Mask is empty, return empty source immediately.
			if constexpr (M == 0) return env::ParticleSource<0, U, false>{};

			// Strategy: We pick the first active field in the Mask to test if 'get_field_ptr_id' exists.
			// We cannot check the function "in general" because it is a template.
			[[maybe_unused]] constexpr auto TestF = static_cast<env::Field>(1 << std::countr_zero(M));

		    // does 'get_field_ptr_id<TestF>(id)' compile?
		    if constexpr (requires { self.template get_field_ptr_id<TestF>(id); }) {

		        // specialized path (direct ID access)
		        constexpr bool IsConst = std::is_const_v<std::remove_reference_t<decltype(self)>>;
		        env::ParticleSource<M, U, IsConst> src;

		        if constexpr (env::has_field_v<M, env::Field::force>)
        			src.force = self.template invoke_get_field_ptr_id<env::Field::force>(id);
		        if constexpr (env::has_field_v<M, env::Field::position>)
        			src.position = self.template invoke_get_field_ptr_id<env::Field::position>(id);
		        if constexpr (env::has_field_v<M, env::Field::velocity>)
        			src.velocity = self.template invoke_get_field_ptr_id<env::Field::velocity>(id);
		        if constexpr (env::has_field_v<M, env::Field::old_position>)
        			src.old_position = self.template invoke_get_field_ptr_id<env::Field::old_position>(id);
		        if constexpr (env::has_field_v<M, env::Field::mass>)
        			src.mass = self.template invoke_get_field_ptr_id<env::Field::mass>(id);
		        if constexpr (env::has_field_v<M, env::Field::state>)
        			src.state = self.template invoke_get_field_ptr_id<env::Field::state>(id);
		        if constexpr (env::has_field_v<M, env::Field::type>)
        			src.type = self.template invoke_get_field_ptr_id<env::Field::type>(id);
		        if constexpr (env::has_field_v<M, env::Field::id>)
        			src.id = self.template invoke_get_field_ptr_id<env::Field::id>(id);
		        if constexpr (env::has_field_v<M, env::Field::user_data>)
        			src.user_data = self.template invoke_get_field_ptr_id<env::Field::user_data>(id);

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
	    env::ParticleID id,
	    size_t index,
	    const env::Box& region,
	    const env::internal::ParticleRecord<typename C::UserData>& p,
	    const std::vector<env::internal::ParticleRecord<typename C::UserData>>& particles
	) {
	    { c.build(particles) };
	    { c.rebuild_structure() };

		{ cc.capacity() } -> std::convertible_to<size_t>;
	    { cc.particle_count() } -> std::convertible_to<size_t>;
		{ cc.min_id() } -> std::convertible_to<env::ParticleID>;
	    { cc.max_id() } -> std::convertible_to<env::ParticleID>;

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
