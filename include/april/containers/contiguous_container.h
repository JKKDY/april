#pragma once

#include "april/particle/particle.h"
#include "april/containers/container.h"

namespace april::container {

	template<typename Config, env::IsUserData U>
	class ContiguousContainer : public Container<Config, U>{
	public:
		using Base = Container<Config, U>;
		friend Base;

		using Particle = env::internal::ParticleRecord<U>;
		using Base::Base;

		void build_storage(const std::vector<Particle>& particles) {
			AP_ASSERT(!is_built, "storage has already been built");

			this->particles = std::vector(particles);
			id_to_index_map.resize(particles.size());
			for (size_t i = 0; i < particles.size(); i++) {
				const auto id = static_cast<size_t>(particles[i].id);
				id_to_index_map[id] = i;
			}
			this->is_built = true;
		}

		// INDEXING
		[[nodiscard]] size_t id_to_index(const env::ParticleID id) const {
			return id_to_index_map[static_cast<size_t>(id)];
		}

		[[nodiscard]] env::ParticleID min_id() const {
			return 0;
		}

		[[nodiscard]] env::ParticleID max_id() const {
			return static_cast<env::ParticleID>(particles.size());
		}


		// QUERIES
		[[nodiscard]] bool contains(const env::ParticleID id) const {
			return id <= max_id();
		}

		[[nodiscard]] size_t particle_count() const {
			return  particles.size();
		}

	protected:
		void swap_particles (uint32_t idx1, uint32_t idx2) {
			auto id1 = particles[idx1].id, id2 = particles[idx2].id;
			std::swap(particles[idx1], particles[idx2]);
			std::swap(id_to_index_map[id1], id_to_index_map[id2]);
		}

		// Deducing 'this' automatically propagates constness to the return type
		template<env::Field F>
		auto get_field_ptr(this auto&& self, size_t i) {
			if constexpr (F == env::Field::force)				return &self.particles[i].force;
			else if constexpr (F == env::Field::position)	  	return &self.particles[i].position;
			else if constexpr (F == env::Field::velocity)	  	return &self.particles[i].velocity;
			else if constexpr (F == env::Field::old_position) 	return &self.particles[i].old_position;
			else if constexpr (F == env::Field::mass)			return &self.particles[i].mass;
			else if constexpr (F == env::Field::state)			return &self.particles[i].state;
			else if constexpr (F == env::Field::type)			return &self.particles[i].type;
			else if constexpr (F == env::Field::id)				return &self.particles[i].id;
			else if constexpr (F == env::Field::user_data)		return &self.particles[i].user_data;
		}

		std::vector<Particle> particles = {};
		std::vector<uint32_t> id_to_index_map; // map id to index
	private:
		bool is_built = false;
	};
}