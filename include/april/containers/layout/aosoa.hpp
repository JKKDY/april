#pragma once


#include <array>
#include <cstddef>
#include <bit>

#include "april/containers/container.hpp"
#include "april/base/types.hpp"
#include "april/particle/properties.hpp"
#include "april/exec/policy.hpp"
#include "../../exec/threading/scheduling.hpp"


#include "april/containers/layout/internal/soa_chunk.hpp"

namespace april::container::layout {

    template <typename ContainerConfig, size_t ChunkSize>
    class AoSoA : public Container<ContainerConfig> {
    public:
        using Base = Container<ContainerConfig>;
        using Base::interaction_map;
        using ParticleAttributes = Base::ParticleAttributes;
        using Base::Base;
        using Base::thread_executor;
        friend Base;

        static constexpr size_t chunk_size = ChunkSize;
        using ChunkT = ParticleChunk<ParticleAttributes, chunk_size>;


        // make inherited accessors explicit otherwise the compiler cant find them due to existing overrides in this class
        using Base::view;
        using Base::at;
        using Base::access_particle;

        explicit AoSoA(const ContainerConfig & config): Base(config) {
            for (size_t k = 0; k < packed::size(); ++k) idx_arr[k] = static_cast<double>(k);
        }

        void bind_executor(Base::ThreadExecutor* raw_executor_ptr) {
            thread_executor.bind(raw_executor_ptr);
            pair_schedule_config = exec::BlockConfig(thread_executor.num_threads(), 2);
            linear_schedule_config = exec::BlockConfig(thread_executor.num_threads(), 8);
            pair_schedule_config.alignment = chunk_size;
            linear_schedule_config.alignment = chunk_size;
        }


        // ACCESSORS (chunk based)
        template <ParticleField Read, ParticleField Write>
        [[nodiscard]] auto at(this auto&& self, size_t chunk_idx, size_t lane_idx) {
            return particle::internal::ScalarParticleRef<Read, Write, ParticleAttributes>{
                self.template access_particle<Read, Write>(chunk_idx, lane_idx)
            };
        }

        template <ParticleField Read>
        [[nodiscard]] auto view(this const auto& self, size_t chunk_idx, size_t lane_idx) {
            return particle::internal::ScalarParticleRef<Read, ParticleField::none, ParticleAttributes>{
                self.template access_particle<Read, ParticleField::none>(chunk_idx, lane_idx)
            };
        }

        template <ParticleField Read, ParticleField Write>
        [[nodiscard]] auto at_packed(this auto&& self, size_t chunk_idx, size_t lane_idx) {
            return particle::internal::PackedParticleRef<Read, Write, ParticleAttributes>{
                self.template access_particle<Read, Write>(chunk_idx, lane_idx)
            };
        }

        template <ParticleField Read>
        [[nodiscard]] auto view_packed(this const auto& self, size_t chunk_idx, size_t lane_idx) {
            return particle::internal::PackedParticleRef<Read, ParticleField::none, ParticleAttributes>{
                self.template access_particle<Read, ParticleField::none>(chunk_idx, lane_idx)
            };
        }


        // INDEXING
        [[nodiscard]] size_t id_to_index(const ParticleID id) const {
            return id_to_index_map[static_cast<size_t>(id)];
        }

        [[nodiscard]] ParticleID min_id() const {
            return 0;
        }

        [[nodiscard]] ParticleID max_id() const {
            return static_cast<ParticleID>(n_particles);
        }

        [[nodiscard]] size_t capacity() const {
            return particle_capacity;
        }


        // PREFETCHING
      // ----------------
       // AOSOA PREFETCHING
       // ----------------

       template <ParticleField Mask>
       APRIL_FORCE_INLINE void prefetch_packed(this const auto& self, size_t c, size_t lane_idx = 0) {
           // Standard prefetch (Temporal Locality = 3) - For outer loops
           const auto& chunk = self.chunks[c]; // Adjust this based on how your container stores chunks

           if constexpr (particle::internal::has_field_v<Mask, ParticleField::force>) {
               APRIL_PREFETCH(&chunk.frc_x[lane_idx]);
               APRIL_PREFETCH(&chunk.frc_y[lane_idx]);
               APRIL_PREFETCH(&chunk.frc_z[lane_idx]);
           }
           if constexpr (particle::internal::has_field_v<Mask, ParticleField::position>) {
               APRIL_PREFETCH(&chunk.pos_x[lane_idx]);
               APRIL_PREFETCH(&chunk.pos_y[lane_idx]);
               APRIL_PREFETCH(&chunk.pos_z[lane_idx]);
           }
           if constexpr (particle::internal::has_field_v<Mask, ParticleField::velocity>) {
               APRIL_PREFETCH(&chunk.vel_x[lane_idx]);
               APRIL_PREFETCH(&chunk.vel_y[lane_idx]);
               APRIL_PREFETCH(&chunk.vel_z[lane_idx]);
           }
           if constexpr (particle::internal::has_field_v<Mask, ParticleField::old_position>) {
               APRIL_PREFETCH(&chunk.old_x[lane_idx]);
               APRIL_PREFETCH(&chunk.old_y[lane_idx]);
               APRIL_PREFETCH(&chunk.old_z[lane_idx]);
           }
           if constexpr (particle::internal::has_field_v<Mask, ParticleField::mass>) {
               APRIL_PREFETCH(&chunk.mass[lane_idx]);
           }
           if constexpr (particle::internal::has_field_v<Mask, ParticleField::state>) {
               APRIL_PREFETCH(&chunk.state[lane_idx]);
           }
           if constexpr (particle::internal::has_field_v<Mask, ParticleField::type>) {
               APRIL_PREFETCH(&chunk.type[lane_idx]);
           }
           if constexpr (particle::internal::has_field_v<Mask, ParticleField::id>) {
               APRIL_PREFETCH(&chunk.id[lane_idx]);
           }
           if constexpr (particle::internal::has_field_v<Mask, ParticleField::attributes>) {
               APRIL_PREFETCH(&chunk.attributes[lane_idx]);
           }
       }

       template <ParticleField Mask>
       APRIL_FORCE_INLINE void prefetch_packed_nta(this const auto& self, size_t c, size_t lane_idx = 0) {
           // NTA prefetch (Locality = 0) - For inner streaming loops
           const auto& chunk = self.chunks[c];

           if constexpr (particle::internal::has_field_v<Mask, ParticleField::force>) {
               APRIL_PREFETCH_NTA(&chunk.frc_x[lane_idx]);
               APRIL_PREFETCH_NTA(&chunk.frc_y[lane_idx]);
               APRIL_PREFETCH_NTA(&chunk.frc_z[lane_idx]);
           }
           if constexpr (particle::internal::has_field_v<Mask, ParticleField::position>) {
               APRIL_PREFETCH_NTA(&chunk.pos_x[lane_idx]);
               APRIL_PREFETCH_NTA(&chunk.pos_y[lane_idx]);
               APRIL_PREFETCH_NTA(&chunk.pos_z[lane_idx]);
           }
           if constexpr (particle::internal::has_field_v<Mask, ParticleField::velocity>) {
               APRIL_PREFETCH_NTA(&chunk.vel_x[lane_idx]);
               APRIL_PREFETCH_NTA(&chunk.vel_y[lane_idx]);
               APRIL_PREFETCH_NTA(&chunk.vel_z[lane_idx]);
           }
           if constexpr (particle::internal::has_field_v<Mask, ParticleField::old_position>) {
               APRIL_PREFETCH_NTA(&chunk.old_x[lane_idx]);
               APRIL_PREFETCH_NTA(&chunk.old_y[lane_idx]);
               APRIL_PREFETCH_NTA(&chunk.old_z[lane_idx]);
           }
           if constexpr (particle::internal::has_field_v<Mask, ParticleField::mass>) {
               APRIL_PREFETCH_NTA(&chunk.mass[lane_idx]);
           }
           if constexpr (particle::internal::has_field_v<Mask, ParticleField::state>) {
               APRIL_PREFETCH_NTA(&chunk.state[lane_idx]);
           }
           if constexpr (particle::internal::has_field_v<Mask, ParticleField::type>) {
               APRIL_PREFETCH_NTA(&chunk.type[lane_idx]);
           }
           if constexpr (particle::internal::has_field_v<Mask, ParticleField::id>) {
               APRIL_PREFETCH_NTA(&chunk.id[lane_idx]);
           }
           if constexpr (particle::internal::has_field_v<Mask, ParticleField::attributes>) {
               APRIL_PREFETCH_NTA(&chunk.attributes[lane_idx]);
           }
       }


        // QUERIES
        [[nodiscard]] bool index_is_valid(const size_t index) const {
            if (index >= particle_capacity) return false;

            auto [c, l] = locate(index);
            return data[c].state[l] != ParticleState::INVALID;
        }

        [[nodiscard]] bool contains_id(const ParticleID id) const {
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
        ChunkT* APRIL_RESTRICT ptr_chunks = nullptr;

        size_t particle_capacity{};
        size_t n_particles{};
        std::vector<ChunkT> data;
        std::vector<ChunkT> tmp;
        std::vector<size_t> bin_starts; // first chunk index of each bin
        std::vector<size_t> bin_sizes; // number of particles in each bin
        std::vector<uint32_t> id_to_index_map;

        exec::BlockConfig pair_schedule_config;
        exec::BlockConfig linear_schedule_config;

        void update_cache() {
            ptr_chunks = data.data();
        }

        void build_storage(const std::vector<particle::ParticleRecord<ParticleAttributes>>& particles) {
            n_particles = particles.size();

            const size_t n_chunks = (n_particles + chunk_size - 1) / chunk_size;
            particle_capacity = n_chunks * chunk_size;

            data.resize(n_chunks);
            id_to_index_map.resize(n_particles);

            bin_sizes.resize(1);
            bin_starts.resize(1);
            bin_sizes[0] = particles.size();
            bin_starts[0] = 0;

            update_cache();

            for (size_t i = 0; i < n_particles; ++i) {
                const auto& p = particles[i];

                // locate destination
                const auto [c_idx, l_idx] = locate(i);
                auto& chunk = data[c_idx];

                // fill chunk
                chunk.insert_particle(l_idx, p);

                // Map ID
                id_to_index_map[static_cast<size_t>(p.id)] = i;
            }

            for (size_t i = n_particles; i < particle_capacity; ++i) {
                const auto [c_idx, l_idx] = locate(i);
                data[c_idx].state[l_idx] = ParticleState::INVALID;
                data[c_idx].id[l_idx] = std::numeric_limits<ParticleID>::max();
                data[c_idx].type[l_idx] = std::numeric_limits<ParticleType>::max();
                data[c_idx].pos_x[l_idx] = 1e50;
                data[c_idx].pos_y[l_idx] = 1e50;
                data[c_idx].pos_z[l_idx] = 1e50;
                data[c_idx].mass[l_idx]  = 1.0;
            }
        }


        std::vector<size_t> bin_counts;
        std::vector<size_t> cached_bins;
        std::vector<size_t> bin_counts_tls_buffers;
        std::vector<size_t> bin_particles;

        template <typename HashFunc>
        requires std::invocable<HashFunc, size_t> &&
            std::unsigned_integral<std::invoke_result_t<HashFunc, size_t>>
        void reorder_storage(const size_t n_bins, HashFunc&& calc_bin) {
            const size_t n_threads = this->thread_executor.num_threads();
            const size_t old_capacity = this->particle_capacity;

            if (n_particles == 0 || n_bins == 0) return;

            // resize buffers
            bin_starts.resize(n_bins);
            bin_sizes.resize(n_bins);
            bin_counts.resize(n_bins);
            cached_bins.resize(old_capacity);
            bin_counts_tls_buffers.assign(n_bins * n_threads, 0);

            // schedule over the capacity (including holes)
            auto capacity_blocks = exec::make_linear_schedule(math::Range{0, old_capacity}, this->linear_schedule_config);
            auto bin_blocks = exec::make_linear_schedule(math::Range{0, n_bins}, this->linear_schedule_config);

            // each thread counts the number of particles per bin in its assigned particle blocks
            this->thread_executor.execute(capacity_blocks.size(), [&](const size_t t_idx) {
                APRIL_ASSERT(exec::thread_index() < this->thread_executor.num_threads(),
                      "[APRIL] exec::thread_index() must be in range [0, executor.num_threads()]. "
                      "Verify that the executor sets ScopedThreadContext correctly");

                const auto& block = capacity_blocks[t_idx];
                auto* bin_counts_tls = &bin_counts_tls_buffers[exec::thread_index() * n_bins];

                for (size_t i = block.start; i < block.stop; ++i) {
                    // skip sentinel/invalid particles inline
                    const auto [c, l] = locate(i);
                    if (data[c].state[l] == ParticleState::INVALID) continue;

                    const size_t bin = calc_bin(i);
                    cached_bins[i] = bin;
                    ++bin_counts_tls[bin];
                }
            });

            // reduce to a single buffer
            this->thread_executor.execute(bin_blocks.size(), [&](const size_t t_idx) {
                const auto& bins = bin_blocks[t_idx];

                for (const size_t bin : bins) bin_sizes[bin] = 0;

                for (unsigned i = 0; i < n_threads; i++) {
                    for (const size_t bin : bins) {
                        bin_sizes[bin] += bin_counts_tls_buffers[n_bins * i + bin];
                    }
                }
            });

            // compute chunk aligned offsets (prefix sum) (serial)
            size_t current_chunk_offset = 0;
            size_t total_chunks = 0;

            for (size_t i = 0; i < n_bins; ++i) {
                bin_starts[i] = current_chunk_offset;
                bin_counts[i] = current_chunk_offset;

                if (bin_sizes[i] > 0) {
                    // Calculate chunks required for this bin
                    const size_t chunks_for_bin = (bin_sizes[i] + chunk_size - 1) / chunk_size;

                    total_chunks += chunks_for_bin;
                    current_chunk_offset += chunks_for_bin * chunk_size;
                }
            }

            // calcualte new capacity and update buffers
            const size_t new_capacity = total_chunks * chunk_size;
            tmp.resize(total_chunks);
            bin_particles.assign(new_capacity, std::numeric_limits<size_t>::max());

            // Build the mapping: destination_idx -> source_idx (serial)
            for (size_t i = 0; i < old_capacity; ++i) {
                const auto [c, l] = locate(i);
                if (data[c].state[l] == ParticleState::INVALID) continue;

                const size_t bin = cached_bins[i];
                const size_t dest_idx = bin_counts[bin]++;
                bin_particles[dest_idx] = i;
            }

            // gather copy into ping pong buffer (chunk by chunk)
            auto chunk_blocks = exec::make_linear_schedule(math::Range{0, total_chunks}, this->linear_schedule_config);

            this->thread_executor.execute(chunk_blocks.size(), [&](const size_t t_idx) {
                 const auto& block = chunk_blocks[t_idx];

                 // Loop over destination CHUNKS
                 for (size_t dst_c = block.start; dst_c < block.stop; ++dst_c) {
                     auto& dst_chunk = tmp[dst_c];

                     // Loop over destination LANES
                     for (size_t dst_l = 0; dst_l < chunk_size; ++dst_l) {
                         const size_t dest_idx = dst_c * chunk_size + dst_l;
                         const size_t src_idx = bin_particles[dest_idx];

                         if (src_idx == std::numeric_limits<size_t>::max()) {
                             // Write Sentinel Padding
                             dst_chunk.state[dst_l] = ParticleState::INVALID;
                             dst_chunk.pos_x[dst_l] = 1e50;
                             dst_chunk.pos_y[dst_l] = 1e50;
                             dst_chunk.pos_z[dst_l] = 1e50;
                             dst_chunk.id[dst_l]    = std::numeric_limits<ParticleID>::max();
                             dst_chunk.type[dst_l]  = std::numeric_limits<ParticleType>::max();
                             dst_chunk.mass[dst_l]  = 1.0;
                         } else {
                             // Gather Valid Particle
                             const auto [src_c, src_l] = locate(src_idx);
                             dst_chunk.copy_from(dst_l, src_l, data[src_c]);

                             // Update ID map
                             const auto id = static_cast<size_t>(dst_chunk.id[dst_l]);
                             id_to_index_map[id] = static_cast<uint32_t>(dest_idx);
                         }
                     }
                 }
             });

            // Swap and update structural tracking
            std::swap(data, tmp);
            this->particle_capacity = new_capacity;
            update_cache();
        }


        // return physical index range
        [[nodiscard]] math::Range get_physical_bin_range(const size_t type) const {
            size_t start = bin_starts[type];
            size_t end = start + bin_sizes[type]; // end is exact start + count (excludes padding)

            return {start, end};
        }

        [[nodiscard]] std::pair<size_t, size_t> locate(const size_t physical_index) const {
            return {
                physical_index >> chunk_shift, // chunk_index = physical_index / chunk_size
                physical_index & chunk_mask // lane_index = physical_index % chunk_size
            };
        }

        template <ParticleField F>
        APRIL_FORCE_INLINE auto get_field_ptr(this auto&& self, size_t i) {
            // locate data, then forward to the fast path
            const auto [chunk_idx, lane_idx] = self.locate(i);
            return self.template get_field_ptr<F>(chunk_idx, lane_idx);
        }

        template <ParticleField F>
        APRIL_FORCE_INLINE auto get_field_ptr(this auto&& self, size_t chunk_idx, size_t lane_idx) {            // locate data
            auto& chunk = self.ptr_chunks[chunk_idx];

            // return vector pointer
            if constexpr (F == ParticleField::position)
                return math::Vec3Ptr{&chunk.pos_x[lane_idx], &chunk.pos_y[lane_idx], &chunk.pos_z[lane_idx]};
            else if constexpr (F == ParticleField::velocity)
                return math::Vec3Ptr{&chunk.vel_x[lane_idx], &chunk.vel_y[lane_idx], &chunk.vel_z[lane_idx]};
            else if constexpr (F == ParticleField::force)
                return math::Vec3Ptr{&chunk.frc_x[lane_idx], &chunk.frc_y[lane_idx], &chunk.frc_z[lane_idx]};
            else if constexpr (F == ParticleField::old_position)
                return math::Vec3Ptr{&chunk.old_x[lane_idx], &chunk.old_y[lane_idx], &chunk.old_z[lane_idx]};

            // return scalar pointer
            else if constexpr (F == ParticleField::mass) return &chunk.mass[lane_idx];
            else if constexpr (F == ParticleField::state) return &chunk.state[lane_idx];
            else if constexpr (F == ParticleField::type) return &chunk.type[lane_idx];
            else if constexpr (F == ParticleField::id) return &chunk.id[lane_idx];
            else if constexpr (F == ParticleField::attributes) return &chunk.attributes[lane_idx];
        }

    private:
        alignas(64) packed::value_type idx_arr[packed::size()]{}; // for creating packed masks quickly

        template <ParallelPolicy P, exec::ExecutionMode V, bool is_const, exec::IsKernel Kernel>
        void iterate_range(this auto&& self, Kernel&& kernel, const size_t start, const size_t end) {
            // route scalar/vector execution
            auto process_sub_range = [&](const size_t r_start, const size_t r_end) APRIL_FORCE_INLINE {
                if constexpr (V == exec::ExecutionMode::Scalar) {
                    self.template iterate_range_scalar<P, is_const>(kernel, r_start, r_end);
                } else if constexpr (V == exec::ExecutionMode::Packed ||
                    V == (exec::ExecutionMode::Scalar | exec::ExecutionMode::Packed)) {
                    self.template iterate_range_vector<P, is_const>(kernel, r_start, r_end);
                } else {
                    static_assert(false,"[APRIL] invalid ExecutionMode in AoSoA::iterate_range");
                }
            };

            if constexpr (P == ParallelPolicy::Serial) {
                process_sub_range(start, end);
            }
            else if constexpr (P == ParallelPolicy::Threaded) {
                math::Range range = {start, end};
                const auto blocks = exec::make_linear_schedule(range, self.linear_schedule_config);

                self.thread_executor.execute(blocks.size(), [&](const size_t i) {
                    process_sub_range(blocks[i].start, blocks[i].stop);
                });
            }
        }

        template <ParallelPolicy P, bool is_const, exec::IsKernel Kernel>
        APRIL_FORCE_INLINE void iterate_range_scalar(this auto&& self, Kernel&& kernel, const size_t start, const size_t end) {
            using K = std::remove_cvref_t<Kernel>;
            if (start >= end) return;

            const auto [start_chunk, start_idx] = self.locate(start);
            const auto [end_chunk, end_idx] = self.locate(end);

            auto* APRIL_RESTRICT chunks = self.ptr_chunks;
            auto exec_scalar = [&](size_t c, size_t i) APRIL_FORCE_INLINE {
                const size_t physical_idx = (c << chunk_shift) | i;
                if constexpr (is_const) kernel(physical_idx, self.template view<K::Read>(c, i));
                else kernel(physical_idx, self.template at<K::Read, K::Write>(c, i));
            };

            if (start_chunk == end_chunk) {
                APRIL_PREFETCH(chunks + start_chunk);
                for (size_t i = start_idx; i < end_idx; ++i) {
                    exec_scalar(start_chunk, i);
                }
            }
            else {
                APRIL_PREFETCH(chunks + start_chunk);
                if (start_chunk + 1 < end_chunk || end_idx > 0) APRIL_PREFETCH(chunks + start_chunk + 1);

                // head
                for (size_t i = start_idx; i < self.chunk_size; ++i) {
                    exec_scalar(start_chunk, i);
                }

                // body
                for (size_t c = start_chunk + 1; c < end_chunk; ++c) {
                    APRIL_PREFETCH(chunks + c + 1);
                    for (size_t i = 0; i < self.chunk_size; ++i) exec_scalar(c, i);
                }

                // tail
                if (end_idx > 0) {
                    for (size_t i = 0; i < end_idx; ++i) {
                        exec_scalar(end_chunk, i);
                    }
                }
            }
        }

        template <ParallelPolicy P, bool is_const, exec::IsKernel Kernel>
        APRIL_FORCE_INLINE void iterate_range_vector(this auto&& self, Kernel&& kernel, const size_t start, const size_t end) {
            using K = std::remove_cvref_t<Kernel>;
            if (start >= end) return;

            // get the first and last chunks touched in the iteration
            const auto [start_chunk, start_idx] = self.locate(start);
            const auto [end_chunk, end_idx] = self.locate(end);

            constexpr size_t simd_width = packed::size();
            auto* APRIL_RESTRICT chunks = self.ptr_chunks;

            // Align backwards to the nearest vector boundary for head masking
            const size_t head_offset = start_idx % simd_width;
            const size_t aligned_start_idx = start_idx - head_offset;

            // kernel expects the logical base index of the vector register
            const auto lane_indices = packed::load_aligned(self.idx_arr);

            auto exec_vector = [&](size_t c, size_t i) APRIL_FORCE_INLINE {
                const size_t physical_idx = (c << chunk_shift) | i;
                if constexpr (is_const) kernel(physical_idx, self.template view_packed<K::Read>(c, i));
                else kernel(physical_idx, self.template at_packed<K::Read, K::Write>(c, i));
            };

            auto exec_vector_masked = [&](size_t c, size_t i, auto mask) APRIL_FORCE_INLINE {
                const size_t physical_idx = (c << chunk_shift) | i;
                if constexpr (is_const) {
                    auto ref = self.template view_packed<K::Read>(c, i);
                    kernel(physical_idx, ref.mask_with(mask));
                } else {
                    auto ref = self.template at_packed<K::Read, K::Write>(c, i);
                    kernel(physical_idx, ref.mask_with(mask));
                }
            };

            // Single Chunk iteration (start and end chunk are the same => only one chunk is touched)
            if (start_chunk == end_chunk) {
                APRIL_PREFETCH(chunks + start_chunk);

                // Edge Case: The entire range fits inside a single vector register
                if (aligned_start_idx + simd_width >= end_idx) {
                    const auto mask = (lane_indices >= static_cast<packed::value_type>(head_offset)) &&
                                      (lane_indices < static_cast<packed::value_type>(end_idx - aligned_start_idx));
                    exec_vector_masked(start_chunk, aligned_start_idx, mask);
                    return;
                }

                // Head
                if (head_offset > 0) {
                    const auto mask = lane_indices >= static_cast<packed::value_type>(head_offset);
                    exec_vector_masked(start_chunk, aligned_start_idx, mask);
                }

                // Body
                const size_t body_start = head_offset > 0 ? aligned_start_idx + simd_width : aligned_start_idx;
                const size_t body_end = end_idx - (end_idx % simd_width);
                for (size_t i = body_start; i < body_end; i += simd_width) {
                    exec_vector(start_chunk, i);
                }

                // Tail
                const size_t tail_count = end_idx % simd_width;
                if (tail_count > 0) {
                    const auto mask = lane_indices < static_cast<packed::value_type>(tail_count);
                    exec_vector_masked(start_chunk, body_end, mask);
                }
            }
            // process multiple chunks
            else {
                APRIL_PREFETCH(chunks + start_chunk);
                if (start_chunk + 1 < end_chunk || end_idx > 0) APRIL_PREFETCH(chunks + start_chunk + 1);

                // Head Chunk
                if (head_offset > 0) {
                    const auto mask = lane_indices >= static_cast<packed::value_type>(head_offset);
                    exec_vector_masked(start_chunk, aligned_start_idx, mask);
                }
                const size_t head_body_start = head_offset > 0 ? aligned_start_idx + simd_width : aligned_start_idx;
                for (size_t i = head_body_start; i < self.chunk_size; i += simd_width) {
                    exec_vector(start_chunk, i);
                }

                // Body Chunks (Guaranteed Aligned)
                for (size_t c = start_chunk + 1; c < end_chunk; ++c) {
                    APRIL_PREFETCH(chunks + c + 1);
                    for (size_t i = 0; i < self.chunk_size; i += simd_width) {
                        exec_vector(c, i);
                    }
                }

                // Tail Chunk
                if (end_idx > 0) {
                    const size_t tail_body_end = end_idx - (end_idx % simd_width);
                    for (size_t i = 0; i < tail_body_end; i += simd_width) {
                        exec_vector(end_chunk, i);
                    }
                    const size_t tail_count = end_idx % simd_width;
                    if (tail_count > 0) {
                        const auto mask = lane_indices < static_cast<packed::value_type>(tail_count);
                        exec_vector_masked(end_chunk, tail_body_end, mask);
                    }
                }
            }
        }
    };
}


