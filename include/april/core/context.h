#pragma once

#include <vector>
#include <cstddef>
#include "april/env/particle.h"
#include "april/env/domain.h"


namespace april::core {

	template<class System>
	class SystemContext {
		using ParticleID = env::ParticleID;
	public:
		using user_data_t = typename System::user_data_t;
		template<env::FieldMask M> using ParticleRef    		= typename System::template ParticleRef<M>;
		template<env::FieldMask M> using ParticleView   		= typename System::template ParticleView<M>;
		template<env::FieldMask M> using RestrictedParticleRef	= typename System::template RestrictedParticleRef<M>;

		explicit SystemContext(System & sys): system(sys) {}

		// ---- Core information ----
		[[nodiscard]] env::Box box() const noexcept { return system.box(); }
		[[nodiscard]] double time() const noexcept { return system.time(); }
		[[nodiscard]] size_t step() const noexcept { return system.step(); }
		[[nodiscard]] size_t size(env::ParticleState state = env::ParticleState::ALL) const noexcept {
			return system.size(state);
		}


		// ---- Particle modification ----
		void register_particle_movement(ParticleID id) { system.register_particle_movement(id); }
		void register_all_particle_movements() { system.register_all_particle_movements(); }


		// ---- Region Query ----
		[[nodiscard]] std::vector<size_t> collect_indices_in_region(const env::Box& region) const {
			return system.collect_indices_in_region(region);
		}


		// ---- Particle access by ID ----
		template<env::FieldMask M>
		[[nodiscard]]  ParticleRef<M> get_particle_by_id(ParticleID id) noexcept {
			return system.template get_particle_by_id<M>(id);
		}

		template<env::FieldMask M>
		[[nodiscard]] ParticleView<M> get_particle_by_id(ParticleID id) const noexcept {
			return std::as_const(system).template get_particle_by_id<M>(id);
		}

		[[nodiscard]] std::vector<env::ParticleID> particle_id_list() const noexcept {
			AP_ASSERT(false, "particle_id_list not implemented yet");
			return {};
		}


		// ---- Particle access by index ----
		template<env::FieldMask M>
		[[nodiscard]] ParticleRef<M> get_particle_by_index(size_t index) noexcept {
			return system.template get_particle_by_index<M>(index);
		}

		template<env::FieldMask M>
		[[nodiscard]] ParticleView<M> get_particle_by_index(size_t index) const noexcept {
			return std::as_const(system).template get_particle_by_index<M>(index);
		}

		[[nodiscard]] size_t index_start() const noexcept {
			return system.index_start();
		}

		[[nodiscard]] size_t index_end() const noexcept {
			return system.index_end();
		}

	private:
		System& system;
	};

}