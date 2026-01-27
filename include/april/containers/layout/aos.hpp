#pragma once

#include "april/particle/particle.hpp"
#include "april/containers/container.hpp"
#include "april/containers/batching.hpp"

namespace april::container::layout {

	template<typename Config, env::IsUserData U>
	class AoS : public Container<Config, U>{
	public:
		using Base = Container<Config, U>;
		using Base::force_schema;
		using Base::Base;
		friend Base;

		using Particle = env::internal::ParticleRecord<U>;

		AoS(const Config & config, const internal::ContainerCreateInfo & info):
			Container<Config, U>(config, info)
		{
			// precompute topology batches (id based batches)
			for (size_t i = 0; i < force_schema.interactions.size(); ++i) {
				const auto& prop = force_schema.interactions[i];

				if (!prop.used_by_ids.empty() && prop.is_active) {
					batching::TopologyBatch batch;
					batch.id1 = prop.used_by_ids[0].first;
					batch.id2 = prop.used_by_ids[0].second;
					batch.pairs = prop.used_by_ids;

					topology_batches.push_back(std::move(batch));
				}
			}
		}

		template<typename Func>
		void for_each_topology_batch(Func && func) {
			for (const auto & batch : topology_batches) {
				func(batch);
			}
		}

		template<env::FieldMask M, ExecutionPolicy Policy, bool is_const, typename Kernel>
		void iterate(this auto&& self, Kernel && kernel, const env::ParticleState state) {
			constexpr env::FieldMask fields = M | env::Field::state;

			for (size_t i = 0; i < self.capacity(); i++) {
				auto raw_state = self.particles[i].state;
				if (static_cast<int>(raw_state & (state & ~env::ParticleState::INVALID))) {
					if constexpr (is_const) {
						kernel(i, self.template view<fields>(i));
					} else {
						kernel(i, self.template at<fields>(i));
					}
				}
			}
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
		[[nodiscard]] bool contains_id(const env::ParticleID id) const {
			return id <= max_id();
		}
		[[nodiscard]] bool index_is_valid(const size_t index) const {
			return index < particle_count();
		}


		// QUERIES
		[[nodiscard]] size_t capacity() const {
			return particle_count();
		}
		[[nodiscard]] size_t particle_count() const {
			return  particles.size();
		}

	protected:
		std::vector<Particle> tmp = {};
		std::vector<Particle> particles = {};
		std::vector<uint32_t> id_to_index_map; // map id to index

		void build_storage(const std::vector<Particle>& particles_in) {
			particles = std::vector(particles_in);
			id_to_index_map.resize(particles.size());
			for (size_t i = 0; i < particles.size(); i++) {
				const auto id = static_cast<size_t>(particles[i].id);
				id_to_index_map[id] = i;
			}

			tmp.resize(particles.size());
		}

		void reorder_storage(const std::vector<std::vector<size_t>> & bins, const bool = false) {
			// scatter particles into bins
			size_t current_idx = 0;
			for (const auto& bin : bins) {
				for (size_t old_idx : bin) {
					tmp[current_idx++] = particles[old_idx];
				}
			}
			// and swap storages
			std::swap(particles, tmp);

			// rebuild id map: we must update the look-up table because every particle moved
			for (size_t i = 0; i < particles.size(); i++) {
				const auto id = particles[i].id;
				id_to_index_map[id] = i;
			}
		}

		void swap_particles (uint32_t i, uint32_t j) {
			if (i == j) return;

			auto id1 = particles[i].id, id2 = particles[j].id;
			std::swap(particles[i], particles[j]);
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

	private:
		std::vector<batching::TopologyBatch> topology_batches;
	};
}