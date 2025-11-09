#pragma once
#include <vector>

#include "april/particle/fields.h"
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

		/// Containers must implement the following functions:
		///   build();
		///   register_all_particle_movements()
		///   register_particle_movement()
		///   calculate_forces()
		///   get_particle_by_id, id_start, id_end
		///   id_to_index
		///   get_particle_by_index (optional), index_start, index_end
		///   particle_count()
		///   collect_indices_in_region()
		/// in the future:
		///   remove_particle
		///   add_particle
		template<env::IsUserData U, env::IsMutableFetcher MF, env::IsConstFetcher CF>
		class ContainerInterface {
		public:
			static_assert(std::same_as<U, typename MF::user_data_t>, "U must match MF::user_data_t");
			static_assert(std::same_as<U, typename CF::user_data_t>, "U must match CF::user_data_t");

			using ParticleID = env::ParticleID;
			using ParticleType = env::ParticleType;
			using ParticleRecord = env::internal::ParticleRecord<U>;
			using MutableFetcher = MF;
			using ConstFetcher = CF;

			// template<env::FieldMask M> using ParticleView = env::ParticleView<M, U>;
			// template<env::FieldMask M> using ParticleRef = env::ParticleRef<M, U>;

			ContainerInterface() = default;

			// user hook to initialize the container
			// TODO add regions that will be queried in the future so the container can keep track of particles better
			void dispatch_build(this auto&& self, const std::vector<ParticleRecord>& particles) {
		        static_assert(
		            requires { self.build(particles); },
		            "Container subclass must implement: void build(const vector<Particle>&)"
		        );
		        self.build(particles);
		    }

			// update container internals by scanning for all particle movements
			void dispatch_register_all_particle_movements(this auto&& self) {
				static_assert(
					requires { self.register_all_particle_movements(); },
					"Container subclass must implement: void register_particle_movements()"
				);
				self.register_all_particle_movements();
			}

			// update container internals by scanning one particle for its movement
			void dispatch_register_particle_movement(this auto&& self, size_t idx) {
				static_assert(
					requires { self.register_particle_movement(idx); },
					"Container subclass must implement: void register_particle_movements()"
				);
				self.register_particle_movement(idx);
			}

			// calculate inter-particle forces
			void dispatch_calculate_forces(this auto&& self) {
		        static_assert(
		            requires { self.calculate_forces(); },
		            "Container subclass must implement: void calculate_forces()"
		        );
		        self.calculate_forces();
		    }


			size_t dispatch_id_to_index(this auto&& self, ParticleID id) {
				static_assert(
				   requires { { self.id_to_index(id) } -> std::same_as<size_t>; },
				   "Container subclass must implement: void id_to_index(id)"
			   );
				return self.id_to_index(id);
			}

			MutableFetcher dispatch_get_fetcher_by_id(this auto&& self, ParticleID id) noexcept{
				if constexpr ( requires { { self.get_fetcher_by_id(id) } -> std::same_as<MutableFetcher>; }) {
					return self.get_fetcher_by_id(id);
				} else {
					size_t idx = self.dispatch_id_to_index(id);
					return self.dispatch_get_fetcher_by_index(idx);
				}
			}

			ConstFetcher dispatch_get_fetcher_by_id(this const auto & self, ParticleID id) noexcept{
				if constexpr ( requires { { self.get_fetcher_by_id(id) } -> std::same_as<ConstFetcher>; }) {
					return self.get_fetcher_by_id(id);
				} else {
					size_t idx = self.dispatch_id_to_index(id);
					return self.dispatch_get_fetcher_by_index(idx);
				}
			}

			// ids are always dense in [0, N-1]
			// use ids for stable iteration
		    ParticleID dispatch_id_start(this auto&& self) {
		        static_assert(
		            requires { { self.id_start() } -> std::same_as<ParticleID>; },
		            "Container subclass must implement: ParticleID id_start()"
		        );
		        return self.id_start();
		    }

		    ParticleID dispatch_id_end(this auto&& self) {
		        static_assert(
		            requires { { self.id_end() } -> std::same_as<ParticleID>; },
		            "Container subclass must implement: ParticleID id_end()"
		        );
		        return self.id_end();
		    }


			MutableFetcher dispatch_get_fetcher_by_index(this auto&& self, size_t index) noexcept {
		        static_assert(
		            requires { { self.get_fetcher_by_index(index) } -> std::same_as<MutableFetcher>; },
		            "Container subclass must implement: MutableFetcher get_fetcher_by_index(size_t)"
		        );
		        return self.get_fetcher_by_index(index);
		    }

			ConstFetcher dispatch_get_fetcher_by_index(this const auto & self, size_t index) noexcept {
				static_assert(
					requires { { self.get_fetcher_by_index(index) } -> std::same_as<ConstFetcher>; },
					"Container subclass must implement: ConstFetcher get_fetcher_by_index(size_t) const"
				);
				return self.get_fetcher_by_index(index);
			}

			ConstFetcher get_const_fetcher_by_index(size_t index) noexcept {
				return std::as_const(*this).dispatch_get_fetcher_by_index(index);
			}

		    size_t dispatch_index_start(this auto&& self) {
		        static_assert(
		            requires { { self.index_start() } -> std::same_as<size_t>; },
		            "Container subclass must implement: size_t index_start()"
		        );
		        return self.index_start();
		    }

		    size_t dispatch_index_end(this auto&& self) {
		        static_assert(
		            requires { { self.index_end() } -> std::same_as<size_t>; },
		            "Container subclass must implement: size_t index_end()"
		        );
		        return self.index_end();
		    }



		    size_t dispatch_particle_count(this auto&& self) {
		        static_assert(
		            requires { { self.particle_count() } -> std::same_as<size_t>; },
		            "Container subclass must implement: size_t particle_count()"
		        );
		        return self.particle_count();
		    }

			// returns a list of indices to the particles container in region
			std::vector<size_t> dispatch_collect_indices_in_region(this auto&& self, const env::Box & region) {
				static_assert(
				   requires { { self.collect_indices_in_region(region) } -> std::same_as<std::vector<size_t>>; },
				   "Container subclass must implement: size_t particles_in_domain()"
			   );

				return self.collect_indices_in_region(region);
			}

			void dispatch_add_particle(this auto&&, const ParticleRecord &) {
				throw std::logic_error("dispatch_add_particle not implemented yet");
			}

			void dispatch_remove_particle(this auto&&, ParticleID) {
				throw std::logic_error("dispatch_remove_particle not implemented yet");
			}
		};
	} // namespace internal


	template<typename CFG, class ForceTable, env::IsUserData U, env::IsMutableFetcher MF, env::IsConstFetcher CF>
	class Container : public internal::ContainerInterface<U, MF, CF> {
	public:
		using config_type_t = CFG;
		using user_type_t = U;
		using force_table_t = ForceTable;
		using mutable_fetcher_t = MF;
		using const_fetcher = CF;

		Container(const CFG & config,
		  const internal::ContainerFlags & flags,
		  const env::Box & box,
		  ForceTable * force_table)
		: cfg(config), flags(flags), domain(box), force_table(force_table) {}

	protected:
		CFG cfg;
		const internal::ContainerFlags flags;
		env::Box domain;
		ForceTable * force_table{};
	};


	template<typename C> concept IsContainer =
	requires {
		typename C::config_type_t;
		typename C::force_table_t;
		typename C::user_type_t;
		typename C::mutable_fetcher_t;
		typename C::const_fetcher;
	} && requires {
		// Check config_type_t defines an impl template taking (force_table_t, user_type_t)
		typename C::config_type_t::template impl<
			typename C::force_table_t,
			typename C::user_type_t>;
	} && requires {
		// Check impl<force_table_t, user_type_t> resolves to this container type
		requires std::same_as<C, typename C::config_type_t::template impl<typename C::force_table_t, typename C::user_type_t>>;
	} && requires {
		// check that C is a derivative of Container
		requires std::derived_from<C, Container<
			typename C::config_type_t,
			typename C::force_table_t,
			typename C::user_type_t,
			typename C::mutable_fetcher_t,
			typename C::const_fetcher
		>>;
	};


	template<typename ContainerDecl, typename Traits> concept IsContainerDecl =
		env::internal::IsEnvironmentTraits<Traits> &&
	requires {
		typename ContainerDecl::template impl<typename Traits::force_table_t, typename Traits::user_data_t>;
	} &&
		IsContainer<typename ContainerDecl::template impl<typename Traits::force_table_t, typename Traits::user_data_t>>;

} // namespace april::container
