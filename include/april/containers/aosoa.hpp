#pragma once


#include <array>
#include <cstddef>
#include <bit>

#include "batching.hpp"
#include "container.hpp"
#include "april/common.hpp"
#include "april/particle/defs.hpp"



namespace april::container::internal {
	template <typename U, size_t Size = 8>
	struct alignas(64) ParticleChunk {
		// Enforce alignment requirements (Size 8 * double 8 bytes = 64 bytes = 1 AVX-512 Register)
		static_assert(std::has_single_bit(Size),
			"Chunk Size must be a Power of 2 (e.g., 8, 16) for bitwise indexing optimizations.");
		static_assert(Size >= 8,
			"Chunk Size must be at least 8 to fill a standard 64-byte Cache Line / AVX-512 register.");

		static constexpr size_t size = Size;

		// Position
		alignas(64) std::array<vec3::type, Size> pos_x;
		alignas(64) std::array<vec3::type, Size> pos_y;
		alignas(64) std::array<vec3::type, Size> pos_z;

		// Velocity
		alignas(64) std::array<vec3::type, Size> vel_x;
		alignas(64) std::array<vec3::type, Size> vel_y;
		alignas(64) std::array<vec3::type, Size> vel_z;

		// Force
		alignas(64) std::array<vec3::type, Size> frc_x;
		alignas(64) std::array<vec3::type, Size> frc_y;
		alignas(64) std::array<vec3::type, Size> frc_z;

		// Old Position (for Verlet)
		alignas(64) std::array<vec3::type, Size> old_x;
		alignas(64) std::array<vec3::type, Size> old_y;
		alignas(64) std::array<vec3::type, Size> old_z;

		// Scalars
		// Note: For types smaller than double (like ParticleType), padding will occur
		// between arrays to maintain 64-byte alignment of the start of the next array.
		alignas(64) std::array<double, Size>             mass;
		alignas(64) std::array<env::ParticleState, Size> state;
		alignas(64) std::array<env::ParticleType, Size>  type;
		alignas(64) std::array<env::ParticleID, Size>    id;
		alignas(64) std::array<U, Size>                  user_data;
	};


	template<typename Chunk>
	struct ChunkedStorage {
		static constexpr size_t CHUNK_SIZE = Chunk::size;
		std::vector<Chunk> chunks;
		size_t n_particles = 0;

		// resize the storage to fit n particles
		void resize(const size_t n) {
			n_particles = n;
			size_t n_chunks = (n + CHUNK_SIZE - 1) / CHUNK_SIZE; // ceiling division: ceil(n / size)
			chunks.resize(n_chunks);
		}

		[[nodiscard]] std::pair<size_t, size_t> locate(const size_t index) const {
			return {index / CHUNK_SIZE, index & (CHUNK_SIZE - 1)};
		}

		// Swap particle i and j (Handles cross-chunk swapping)
		void swap(const size_t i, const size_t j) {
			if (i == j) return;

			// Decode Linear Index -> (Chunk Index, Lane Index)
			auto [c1_i, l1] = locate(i);
			auto [c2_i, l2] = locate(j);

			auto& c1 = chunks[c1_i];
			auto& c2 = chunks[c2_i];

			// swap values
			std::swap(c1.pos_x[l1], c2.pos_x[l2]); std::swap(c1.pos_y[l1], c2.pos_y[l2]); std::swap(c1.pos_z[l1], c2.pos_z[l2]);
			std::swap(c1.vel_x[l1], c2.vel_x[l2]); std::swap(c1.vel_y[l1], c2.vel_y[l2]); std::swap(c1.vel_z[l1], c2.vel_z[l2]);
			std::swap(c1.frc_x[l1], c2.frc_x[l2]); std::swap(c1.frc_y[l1], c2.frc_y[l2]); std::swap(c1.frc_z[l1], c2.frc_z[l2]);
			std::swap(c1.old_x[l1], c2.old_x[l2]); std::swap(c1.old_y[l1], c2.old_y[l2]); std::swap(c1.old_z[l1], c2.old_z[l2]);

			std::swap(c1.mass[l1],      c2.mass[l2]);
			std::swap(c1.state[l1],     c2.state[l2]);
			std::swap(c1.type[l1],      c2.type[l2]);
			std::swap(c1.id[l1],        c2.id[l2]);
			std::swap(c1.user_data[l1], c2.user_data[l2]);
		}
	};
}

namespace april::container {
	template<size_t chunk_size, typename Config, env::IsUserData U>
	class AoSoAContainer : public Container<Config, U> {
	public:
		using Base = Container<Config, U>;
		using Base::force_schema;
		using Base::Base; // Inherit constructors
		friend Base;

		using ChunkType = internal::ParticleChunk<U, chunk_size>;


		AoSoAContainer(const Config & config, const internal::ContainerCreateInfo & info)
			: Base(config, info)
		{
			// Topology batch initialization (identical to SoA/AoS)
			for (size_t i = 0; i < force_schema.interactions.size(); ++i) {
				const auto& prop = force_schema.interactions[i];
				if (!prop.used_by_ids.empty() && prop.is_active) {
					batching::TopologyBatch batch;
					batch.id1 = prop.used_by_ids[0].first;
					batch.id2 = prop.used_by_ids[0].second;
					batch.pairs = prop.used_by_ids;
					topology_batches.push_back(batch);
				}
			}
		}


		template<typename Func>
		void for_each_topology_batch(Func && func) {
			for (const auto & batch : topology_batches) {
				func(batch);
			}
		}

		void build_storage(const std::vector<env::internal::ParticleRecord<U>>& particles) {
			const size_t n = particles.size();
			data.resize(n);
			id_to_index_map.resize(n);

			for (size_t i = 0; i < n; ++i) {
				const auto& p = particles[i];

				// locate destination
				const auto [c_idx, l_idx] = data.locate(i);
				auto& chunk = data.chunks[c_idx];

				// fill chunk
				chunk.pos_x[l_idx] = p.position.x;
				chunk.pos_y[l_idx] = p.position.y;
				chunk.pos_z[l_idx] = p.position.z;

				chunk.vel_x[l_idx] = p.velocity.x;
				chunk.vel_y[l_idx] = p.velocity.y;
				chunk.vel_z[l_idx] = p.velocity.z;

				chunk.frc_x[l_idx] = p.force.x;
				chunk.frc_y[l_idx] = p.force.y;
				chunk.frc_z[l_idx] = p.force.z;

				chunk.old_x[l_idx] = p.old_position.x;
				chunk.old_y[l_idx] = p.old_position.y;
				chunk.old_z[l_idx] = p.old_position.z;

				chunk.mass[l_idx]  = p.mass;
				chunk.state[l_idx] = p.state;
				chunk.type[l_idx]  = p.type;
				chunk.id[l_idx]    = p.id;
				chunk.user_data[l_idx] = p.user_data;

				// Map ID
				id_to_index_map[static_cast<size_t>(p.id)] = i;
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
			return static_cast<env::ParticleID>(id_to_index_map.size());
		}


		// QUERIES
		[[nodiscard]] bool contains(const env::ParticleID id) const {
			return id <= max_id();
		}
		[[nodiscard]] size_t particle_count() const {
			return data.n_particles;
		}

	protected:
		internal::ChunkedStorage<ChunkType> data;
		std::vector<uint32_t> id_to_index_map;
		void swap_particles(size_t i, size_t j) {
			if (i == j) return;

			// get IDs before swap to update map
			const auto [c1, l1] = data.locate(i);
			const auto [c2, l2] = data.locate(j);

			const auto id1 = data.chunks[c1].id[l1];
			const auto id2 = data.chunks[c2].id[l2];

			// swap data & ids
			data.swap(i, j);
			std::swap(id_to_index_map[static_cast<size_t>(id1)],
					  id_to_index_map[static_cast<size_t>(id2)]);
		}

		template<env::Field F>
		auto get_field_ptr(this auto&& self, size_t i) {

			// locate data
			const auto [chunk_idx, lane_idx] = self.data.locate(i);
			auto& chunk = self.data.chunks[chunk_idx];

			// return vector pointer
			if constexpr (F == env::Field::position)
				return utils::Vec3Ptr { &chunk.pos_x[lane_idx], &chunk.pos_y[lane_idx], &chunk.pos_z[lane_idx] };
			else if constexpr (F == env::Field::velocity)
				return utils::Vec3Ptr { &chunk.vel_x[lane_idx], &chunk.vel_y[lane_idx], &chunk.vel_z[lane_idx] };
			else if constexpr (F == env::Field::force)
				return utils::Vec3Ptr { &chunk.frc_x[lane_idx], &chunk.frc_y[lane_idx], &chunk.frc_z[lane_idx] };
			else if constexpr (F == env::Field::old_position)
				return utils::Vec3Ptr { &chunk.old_x[lane_idx], &chunk.old_y[lane_idx], &chunk.old_z[lane_idx] };

			// return scalar pointer
			else if constexpr (F == env::Field::mass)      return &chunk.mass[lane_idx];
			else if constexpr (F == env::Field::state)     return &chunk.state[lane_idx];
			else if constexpr (F == env::Field::type)      return &chunk.type[lane_idx];
			else if constexpr (F == env::Field::id)        return &chunk.id[lane_idx];
			else if constexpr (F == env::Field::user_data) return &chunk.user_data[lane_idx];
		}

	private:
		std::vector<batching::TopologyBatch> topology_batches;
	};
}
