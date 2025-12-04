#pragma once
#include <vector>

#include "april/particle/fields.h"
#include "april/particle/access.h"
#include "april/env/domain.h"



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
	}

	/*
	// -----------------------------------------------------------------------------
	// REFERENCE CONTAINER IMPLEMENTATION
	// Copy-paste this to start a new container.
	// -----------------------------------------------------------------------------
	class ReferenceContainer {
	public:
		using config_type_t     = ...;
		using user_type_t       = ...;
		using mutable_fetcher_t = ...;
		using const_fetcher_t   = ...;

		// --- MANDATORY (Checked by HasContainerOps) ---
		void build(const std::vector<ParticleRecord>& p);
		void rebuild_structure();

		size_t particle_count() const;
		env::ParticleID min_id() const;
		env::ParticleID max_id() const;

		size_t id_to_index(env::ParticleID id) const;
		bool contains(env::ParticleID id) const; // Optional but recommended

		// The Fetcher Hook
		MutableFetcher get_fetcher(size_t index);
		ConstFetcher   get_fetcher(size_t index) const;

		void add_particle(const ParticleRecord& p);
		void remove_particle(env::ParticleID id);

		template<typename F> void for_each_particle(F&& func);
		template<typename F> void for_each_interaction_batch(F&& func);

		std::vector<size_t> collect_indices_in_region(const env::Box& region);

		// --- OPTIONAL OPTIMIZATIONS (Checked via SFINAE / if constexpr) ---

		// 1. Partial Updates (If missing, System falls back to rebuild_structure)
		void notify_moved(const std::vector<size_t>& indices);

		// 2. Fast ID Fetch (If missing, System uses id_to_index + get_fetcher)
		MutableFetcher get_fetcher_by_id(env::ParticleID id);
		ConstFetcher   get_fetcher_by_id(env::ParticleID id) const;
	};
	*/

	template<class C, env::IsUserData U, env::IsFetcher F>
	class Container {
		static_assert(std::same_as<U, typename F::UserDataT>, "U must match MF::user_data_t");
	public:
		using ParticleRecord = env::internal::ParticleRecord<U>;
		using FetcherT = F;
		using UserDataT = U;
		using ConfigT = C;

		Container(const ConfigT & config,
			const internal::ContainerFlags & flags,
			const internal::ContainerHints & hints,
			const env::Box & box)
			: config(config), flags(flags), hints(hints), domain(box)
		{}

		void invoke_build(this auto&& self, const std::vector<ParticleRecord>& particles) {
			self.build(particles);
		}



		// ---------
		// ITERATION
		// ---------
		template<env::FieldMask M, typename Func, bool parallelize=false> // TODO restrict callable (invoke_for_each_particle)
		void invoke_for_each_particle(this auto&& self, Func && func) {
			if constexpr (requires {self.for_each_particle<M>(func); }) {
				self.for_each_particle<M>(func);
			} else if constexpr (parallelize) {
				throw std::logic_error("parallelization not implemented yet");
			} else {
				for (size_t i = 0; i < self.particle_count(); i++) {
					func(self.template at<M>(i));
				}
			}
		}

		template<typename Func> // TODO restrict callable (invoke_for_each_batch)
		void invoke_for_each_interaction_batch(this auto&& self, Func && func) {
			self.for_each_interaction_batch(func);
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
				self.rebuild_structure();
			}
		}


		// ---------
		// MODIFIERS
		// ---------

		void invoke_add_particle(this auto&&, const ParticleRecord &) {
			throw std::logic_error("dispatch_add_particle not implemented yet");
		}

		void invoke_remove_particle(this auto&&, const env::ParticleID) {
			throw std::logic_error("dispatch_remove_particle not implemented yet");
		}


		// ------------------
		// PARTICLE ACCESSORS
		// ------------------
		// INDEX ACCESSORS
		template<env::FieldMask M>
		[[nodiscard]] auto at(this auto&& self, size_t index) {
			return env::ParticleRef<M, UserDataT>{ self.access_fetcher(index) };
		}

		template<env::FieldMask M>
		[[nodiscard]] auto view(this const auto& self, size_t index) {
			return env::ParticleView<M, UserDataT>{ self.access_fetcher(index) };
		}

		template<env::FieldMask M>
		[[nodiscard]] auto restricted_at(this auto&& self, size_t index) {
			return env::RestrictedParticleRef<M, UserDataT>{ self.access_fetcher(index) };
		}

		// ID ACCESSORS
		template<env::FieldMask M>
		[[nodiscard]] auto at_id(this auto&& self, env::ParticleID id) {
			return env::ParticleRef<M, UserDataT>{ self.access_fetcher_by_id(id) };
		}

		template<env::FieldMask M>
		[[nodiscard]] auto view_id(this const auto & self, env::ParticleID id) {
			return env::ParticleView<M, UserDataT>{ self.access_fetcher_by_id(id) };
		}

		template<env::FieldMask M>
		[[nodiscard]] auto restricted_at_id(this auto&& self, env::ParticleID id) {
			return env::RestrictedParticleRef<M, UserDataT>{ self.access_fetcher_by_id(id) };
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


		// -------
		// QUERIES
		// -------
		[[nodiscard]] bool invoke_contains(this const auto& self, env::ParticleID id) {
			return self.contains(id);
		}

		[[nodiscard]] size_t invoke_particle_count(this const auto& self) {
			return self.particle_count();
		}

		[[nodiscard]] std::vector<size_t> invoke_collect_indices_in_region(this const auto& self, const env::Box & region) {
			return self.collect_indices_in_region(region);
		}

		FetcherT access_fetcher(this auto&& self, size_t index) noexcept {
			static_assert(
				requires {
					{ self.get_fetcher(index) } -> std::convertible_to<FetcherT>;
				},
				"get_fetcher(index) not implemented or invalid return type"
			);
			return self.get_fetcher(index);
		}

		FetcherT access_fetcher_by_id(this auto&& self, env::ParticleID id) noexcept{
			if constexpr ( requires { { self.get_fetcher_by_id(id) } -> std::convertible_to<FetcherT>; }) {
				return self.get_fetcher_by_id(id);
			} else {
				size_t idx = self.invoke_id_to_index(id);
				return self.get_fetcher(idx);
			}
		}

	protected:
		ConfigT config;
		const internal::ContainerFlags flags;
		const internal::ContainerHints hints;
		env::Box domain;


		//--------------
		// DATA FETCHERS
		//--------------
	};



	template<typename C>
	concept HasContainerOps = requires (
	    C c,
	    const C cc,
	    env::ParticleID id,
	    size_t index,
	    const env::Box& region,
	    const env::internal::ParticleRecord<typename C::UserDataT>& p,
	    const std::vector<env::internal::ParticleRecord<typename C::UserDataT>>& particles
	) {
	    { c.build(particles) };
	    { c.rebuild_structure() };

	    { cc.particle_count() } -> std::convertible_to<size_t>;
		{ cc.min_id() } -> std::convertible_to<env::ParticleID>;
	    { cc.max_id() } -> std::convertible_to<env::ParticleID>;

	    { cc.id_to_index(id) } -> std::convertible_to<size_t>;
		{ cc.contains(id) } -> std::convertible_to<bool>;

	    { c.collect_indices_in_region(region) } -> std::convertible_to<std::vector<size_t>>;

	    { c.for_each_interaction_batch([](auto&&){}) };

		// for future use
		// { c.add_particle(p) };
		// { c.remove_particle(id) };
	};


	template<typename C> concept IsContainer =
		// 1. must define types (Config, UserData, Fetcher)
		requires {
			typename C::ConfigT;
			typename C::UserDataT;
			typename C::FetcherT;
		} &&
		// 2. Config must have impl typename pointing to Container type
		// container must only depend on user data as template argument
		requires {
			typename C::ConfigT::template impl<typename C::UserDataT>;
			requires std::same_as<C, typename C::ConfigT::template impl<typename C::UserDataT>>;
		} &&
		// 3. Must inherit from the Middleware (Container)
		std::derived_from<C, Container<
			typename C::ConfigT,
			typename C::UserDataT,
			typename C::FetcherT
		>> &&
		// 4. Must implement the Structural Contract
		HasContainerOps<C>;



	template<typename ContainerDecl, typename Traits> concept IsContainerDecl =
		env::internal::IsEnvironmentTraits<Traits>
	&&
	requires {
		typename ContainerDecl::template impl<typename Traits::user_data_t>;
	} &&
		IsContainer<typename ContainerDecl::template impl<typename Traits::user_data_t>>;

} // namespace april::container
