#pragma once
#include <vector>

#include "april/env/particle.h"
#include "april/forces/force_table.h"
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

		// Containers must implement the following functions:
		//   build();
		//   register_all_particle_movements()
		//   register_particle_movement()
		//   calculate_forces()
		//   get_particle_by_id, id_start, id_end
		//   id_to_index
		//   get_particle_by_index (optional), index_start, index_end
		//   particle_count()
		//   collect_indices_in_region()
		// in the future:
		//   remove_particle
		//   add_particle
		template <class Env> class ContainerInterface {
		public:
			using ForceTable = force::internal::ForceTable<Env>;
			using Domain = env::Domain;
			using Particle = env::internal::Particle;
			using ParticleID = env::internal::ParticleID;

			ContainerInterface() = default;
			void init(ForceTable & force_table, const Domain & dom) {
				interactions = &force_table;
				domain = dom;
			}
			virtual ~ContainerInterface() = default;

			// user hook to initialize the container
			void dispatch_build(this auto&& self, const std::vector<Particle>& particles) {
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
				return self.id_to_index();
			}

			// ids are always dense in [0, N-1]
			// use ids for stable iteration
		    Particle& dispatch_get_particle_by_id(this auto&& self, ParticleID id) {
				if constexpr (
					requires { { self.get_particle_by_id(id) } -> std::same_as<Particle&>; }
					) {
					return self.get_particle_by_id(id);
				} else {
					size_t idx = self.dispatch_id_to_index(id);
					return self.dispatch_get_particle_by_index(idx);
				}
		    }

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


		    Particle& dispatch_get_particle_by_index(this auto&& self, size_t index) noexcept {
		        static_assert(
		            requires { { self.get_particle_by_index(index) } -> std::same_as<Particle&>; },
		            "Container subclass must implement: Particle& get_particle_by_index(size_t)"
		        );
		        return self.get_particle_by_index(index);
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

			void dispatch_add_particle(this auto&&, const Particle &) {
				throw std::logic_error("dispatch_add_particle not implemented yet");
			}

			void dispatch_remove_particle(this auto&&, ParticleID) {
				throw std::logic_error("dispatch_remove_particle not implemented yet");
			}

		protected:
			ForceTable * interactions{};
			Domain domain;
		};
	} // namespace internal


	template<typename Config, typename Env> class Container : public internal::ContainerInterface<Env> {
	public:
		using CFG = Config;
		explicit Container(Config config, const internal::ContainerFlags & flags): cfg(config), flags(flags) {}

		using typename internal::ContainerInterface<Env>::Particle;
		using typename internal::ContainerInterface<Env>::ParticleID;
	protected:
		Config cfg;
		const internal::ContainerFlags flags;
	};


	template<typename A> concept IsContainer =
		// std::derived_from<A, Container<typename A::CFG, >> &&
			requires {
		typename A::CFG::impl;
	};


	// TODO implement Container declaration concept
	template<typename A> concept IsContDecl = std::true_type::value;
	// 	requires {
	// 	typename A::impl;
	// }  && IsContainer<typename A::impl>;

} // namespace april::container
