#pragma once


#include <array>
#include <cstddef>
#include <bit>

#include "april/containers/batching.hpp"
#include "april/containers/container.hpp"
#include "april/common.hpp"
#include "april/particle/defs.hpp"



namespace april::container::layout {
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


	template<size_t size, typename Config, env::IsUserData U>
	class AoSoA : public Container<Config, U> {
	public:
		static constexpr size_t chunk_size = size;
		using Chunk = ParticleChunk<U, chunk_size>;

		using Base = Container<Config, U>;
		using Base::force_schema;
		using Base::Base; // Inherit constructors
		friend Base;

		// make inherited accessors explicit.
		// Otherwise the compiler cant find them due to existing overrides in this class
		using Base::view;
		using Base::at;
		using Base::restricted_at;
		using Base::access_particle;



		AoSoA(const Config & config, const internal::ContainerCreateInfo & info)
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


		template<env::FieldMask M, ExecutionPolicy Policy, bool is_const, typename Kernel>
		void iterate(this auto&& self, Kernel && kernel, const env::ParticleState state) {
			constexpr env::FieldMask fields = M | env::Field::state;
			const size_t n_chunks = self.data.size();
			if (n_chunks == 0) return;

			auto iter_chunk = [&](size_t c) {
				for (size_t i = 0; i < self.chunk_size; i++) {

					// peek at state to check if the data is valid or garbage
					auto raw_state = self.data[c].state[i];

					if (static_cast<int>(raw_state & (state & ~env::ParticleState::INVALID))) {
						size_t phys_idx = c * self.chunk_size + i;

						if constexpr (is_const) {// read only
							kernel(phys_idx, self.template view<fields>(c, i));
						} else { // read-write
							kernel(phys_idx, self.template at<fields>(c, i));
						}
					}
				}
			};

			AP_PREFETCH(&self.data[0]);
			for (size_t c = 0; c < n_chunks - 1; c++) {
				AP_PREFETCH(&self.data[c + 1]);
				iter_chunk(c);
			}

			iter_chunk(n_chunks - 1);
		}

		// ACCESSORS (chunk based)
		template<env::FieldMask M>
		[[nodiscard]] auto at(this auto&& self, size_t chunk_idx, size_t lane_idx) {
			return env::ParticleRef<M, U>{ self.template access_particle<M>(chunk_idx, lane_idx) };
		}

		template<env::FieldMask M>
		[[nodiscard]] auto view(this const auto& self, size_t chunk_idx, size_t lane_idx) {
			return env::ParticleView<M, U>{ self.template access_particle<M>(chunk_idx, lane_idx) };
		}

		template<env::FieldMask M>
		[[nodiscard]] auto restricted_at(this auto&& self, size_t chunk_idx, size_t lane_idx) {
			return env::RestrictedParticleRef<M, U>{ self.template access_particle<M>(chunk_idx, lane_idx) };
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
		[[nodiscard]] size_t capacity() const {
			return particle_capacity;
		}


		// QUERIES
		[[nodiscard]] bool index_is_valid(const size_t index) const {
			if (index >= particle_capacity) return false;

			auto [c, l] = locate(index);
			return data[c].state[l] != env::ParticleState::INVALID;
		}
		[[nodiscard]] bool contains_id(const env::ParticleID id) const {
			if (static_cast<size_t>(id) >= id_to_index_map.size()) return false;
			return id_to_index_map[static_cast<size_t>(id)] != ID_NOT_FOUND;
		}
		[[nodiscard]] size_t particle_count() const {
			return n_particles;
		}

	protected:
		static constexpr size_t chunk_mask = chunk_size - 1; // chunk_size is a power of 2
		static constexpr size_t chunk_shift = std::countr_zero(chunk_size); // log_2(chunk_size)
		static constexpr uint32_t ID_NOT_FOUND = std::numeric_limits<uint32_t>::max();

		// hoist data pointer outside and restrict
		Chunk* AP_RESTRICT ptr_chunks = nullptr;

		size_t particle_capacity{};
		size_t n_particles{};
		std::vector<Chunk> data;
		std::vector<Chunk> tmp;
		std::vector<size_t> bin_starts;
		std::vector<uint32_t> id_to_index_map;

		void update_cache() {
			ptr_chunks = data.data();
		}

		void build_storage(const std::vector<env::internal::ParticleRecord<U>>& particles) {
			const size_t n = particles.size();

			const size_t n_chunks = (n + chunk_size - 1) / chunk_size;
			particle_capacity = n_chunks * chunk_size;

			data.resize(n_chunks);
			id_to_index_map.resize(n);

			update_cache();

			for (size_t i = 0; i < n; ++i) {
				const auto& p = particles[i];

				// locate destination
				const auto [c_idx, l_idx] = locate(i);
				auto& chunk = data[c_idx];

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

			for (size_t i = n; i < particle_capacity; ++i) {
				const auto [c_idx, l_idx] = locate(i);
				data[c_idx].state[l_idx] = env::ParticleState::INVALID;
				data[c_idx].id[l_idx] = std::numeric_limits<env::ParticleID>::max();
			}
		}


		void reorder_storage(const std::vector<std::vector<size_t>> & bins, const bool sentinel_pad=true) {
			tmp.clear();
			bin_starts.clear();

			// allocate space, calculate chunk index of first chunk in each bin
			bin_starts.reserve(bins.size() + 1);

			size_t n_chunks = 0;
			for (const auto & bin : bins) {
				bin_starts.push_back(n_chunks);
				if (bin.empty()) continue;
				const size_t bin_size = bin.size();
				n_chunks +=  (bin_size + chunk_size - 1) / chunk_size;
			}

			bin_starts.push_back(n_chunks);
			particle_capacity = n_chunks * chunk_size;
			tmp.resize(n_chunks);
			id_to_index_map.assign(particle_capacity, ID_NOT_FOUND);

			// traversal cursors
			size_t dst_c = 0; // destination chunk index
			size_t dst_l = 0; // destination lane index

			for (const auto & bin : bins) {
				if (bin.empty()) continue;

				for (const size_t src_idx : bin) {
					// locate source
					auto [src_c, src_l] = locate(src_idx);

					auto& src_chunk = data[src_c];
					auto& dst_chunk = tmp[dst_c];

					// copy fields
					dst_chunk.pos_x[dst_l] = src_chunk.pos_x[src_l];
					dst_chunk.pos_y[dst_l] = src_chunk.pos_y[src_l];
					dst_chunk.pos_z[dst_l] = src_chunk.pos_z[src_l];

					dst_chunk.vel_x[dst_l] = src_chunk.vel_x[src_l];
					dst_chunk.vel_y[dst_l] = src_chunk.vel_y[src_l];
					dst_chunk.vel_z[dst_l] = src_chunk.vel_z[src_l];

					dst_chunk.frc_x[dst_l] = src_chunk.frc_x[src_l];
					dst_chunk.frc_y[dst_l] = src_chunk.frc_y[src_l];
					dst_chunk.frc_z[dst_l] = src_chunk.frc_z[src_l];

					dst_chunk.old_x[dst_l] = src_chunk.old_x[src_l];
					dst_chunk.old_y[dst_l] = src_chunk.old_y[src_l];
					dst_chunk.old_z[dst_l] = src_chunk.old_z[src_l];

					dst_chunk.mass[dst_l]  = src_chunk.mass[src_l];
					dst_chunk.state[dst_l] = src_chunk.state[src_l];
					dst_chunk.type[dst_l]  = src_chunk.type[src_l];
					dst_chunk.id[dst_l]    = src_chunk.id[src_l];
					dst_chunk.user_data[dst_l] = src_chunk.user_data[src_l];

					const size_t new_physical_idx = dst_c * chunk_size + dst_l;
					const env::ParticleID id = dst_chunk.id[dst_l];
					id_to_index_map[id] = new_physical_idx;

					dst_l++;
					if (dst_l == chunk_size) {
						dst_l = 0;
						dst_c++;
					}
				}

				// If the bin ended mid-chunk, fill the rest with safe garbage and skip to next chunk.
				if (sentinel_pad && dst_l > 0) {
					auto& dst_chunk = tmp[dst_c];
					while (dst_l < chunk_size) {
						// mark as dead so physics kernels ignore it
						dst_chunk.state[dst_l] = env::ParticleState::INVALID;

						// move far away to be safe against distance checks
						dst_chunk.pos_x[dst_l] = std::numeric_limits<double>::max();
						dst_chunk.pos_y[dst_l] = std::numeric_limits<double>::max();
						dst_chunk.pos_z[dst_l] = std::numeric_limits<double>::max();

						// set ID to max to avoid map lookups
						dst_chunk.id[dst_l] = std::numeric_limits<env::ParticleID>::max();

						dst_l++;
					}
					// chunk is now "full"
					dst_l = 0;
					dst_c++;
				}
			}
			std::swap(data, tmp);
			update_cache();
		}



		[[nodiscard]] std::pair<size_t, size_t> bin_range(const size_t bin_index) const {
			return {bin_starts[bin_index], bin_starts[bin_index+1]};
		}

		[[nodiscard]] std::pair<size_t, size_t> locate(const size_t physical_index) const {
			return {
				physical_index >> chunk_shift, // physical_index / chunk_size
				physical_index & chunk_mask    // physical_index % chunk_size
			 };
		}

		template<env::Field F>
		auto get_field_ptr(this auto&& self, const size_t i) {

			// locate data
			const auto [chunk_idx, lane_idx] = self.locate(i);
			auto& chunk = self.ptr_chunks[chunk_idx];

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

		template<env::FieldMask M>
		[[nodiscard]] auto access_particle(this auto&& self, size_t chunk_idx, size_t lane_idx) {

			constexpr bool IsConst = std::is_const_v<std::remove_reference_t<decltype(self)>>;
			env::ParticleSource<M, U, IsConst> src;

			auto& chunk = self.ptr_chunks[chunk_idx];

			if constexpr (env::has_field_v<M, env::Field::force>)
				src.force = utils::Vec3Ptr { &chunk.frc_x[lane_idx], &chunk.frc_y[lane_idx], &chunk.frc_z[lane_idx] };
			if constexpr (env::has_field_v<M, env::Field::position>)
				src.position = utils::Vec3Ptr { &chunk.pos_x[lane_idx], &chunk.pos_y[lane_idx], &chunk.pos_z[lane_idx] };
			if constexpr (env::has_field_v<M, env::Field::velocity>)
				src.velocity = utils::Vec3Ptr { &chunk.vel_x[lane_idx], &chunk.vel_y[lane_idx], &chunk.vel_z[lane_idx] };
			if constexpr (env::has_field_v<M, env::Field::old_position>)
				src.old_position = utils::Vec3Ptr { &chunk.old_x[lane_idx], &chunk.old_y[lane_idx], &chunk.old_z[lane_idx] };
			if constexpr (env::has_field_v<M, env::Field::mass>)
				src.mass = &chunk.mass[lane_idx];
			if constexpr (env::has_field_v<M, env::Field::state>)
				src.state = &chunk.state[lane_idx];
			if constexpr (env::has_field_v<M, env::Field::type>)
				src.type = &chunk.type[lane_idx];
			if constexpr (env::has_field_v<M, env::Field::id>)
				src.id = &chunk.id[lane_idx];
			if constexpr (env::has_field_v<M, env::Field::user_data>)
				src.user_data = &chunk.user_data[lane_idx];

			return src;
		}
	};
}
