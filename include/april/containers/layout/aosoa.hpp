#pragma once


#include <array>
#include <cstddef>
#include <bit>

#include "april/containers/container.hpp"
#include "april/base/types.hpp"
#include "april/particle/particle_types.hpp"
#include "april/exec/policy.hpp"

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
        ChunkT* AP_RESTRICT ptr_chunks = nullptr;

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

        void reorder_storage(const std::vector<std::vector<size_t>>& new_bins, const bool sentinel_pad = true) {
            bin_starts.clear();
            bin_sizes.clear();

            std::vector<size_t> offsets(new_bins.size());
            size_t n_chunks = 0;

            for (size_t i = 0; i < new_bins.size(); ++i) {
                const size_t start_idx = n_chunks * chunk_size;
                bin_starts.push_back(start_idx);
                bin_sizes.push_back(new_bins[i].size());
                offsets[i] = start_idx;

                if (!new_bins[i].empty()) {
                    n_chunks += (new_bins[i].size() + chunk_size - 1) / chunk_size;
                }
            }

            particle_capacity = n_chunks * chunk_size;
            bin_starts.push_back(particle_capacity);
            tmp.resize(n_chunks);

            auto bin_blocks = exec::make_linear_schedule(math::Range{0, new_bins.size()}, this->linear_schedule_config);

            this->thread_executor.execute(bin_blocks.size(), [&](const size_t b_idx) {
                const auto& block = bin_blocks[b_idx];

                for (size_t bin_idx = block.start; bin_idx < block.stop; ++bin_idx) {
                    const auto& bin = new_bins[bin_idx];
                    if (bin.empty()) continue;

                    const size_t start_offset = offsets[bin_idx];

                    size_t current_new_idx = start_offset;
                    size_t dst_c = start_offset / chunk_size;
                    size_t dst_l = 0;

                    for (const size_t old_idx : bin) {
                        auto [src_c, src_l] = locate(old_idx);

                        auto& src_chunk = data[src_c];
                        auto& dst_chunk = tmp[dst_c];

                        // Copy particle fields manually
                        dst_chunk.copy_from(dst_l, src_l, src_chunk);

                        id_to_index_map[dst_chunk.id[dst_l]] = current_new_idx;

                        ++current_new_idx;
                        ++dst_l;
                        if (dst_l == chunk_size) {
                            dst_l = 0;
                            ++dst_c;
                        }
                    }

                    // Sentinel pad logic
                    if (sentinel_pad) {
                        const size_t remainder = bin.size() % chunk_size;
                        if (remainder > 0) {
                            // dst_c is already pointing at the correct final chunk
                            auto& dst_chunk = tmp[dst_c];

                            for (size_t pad_l = remainder; pad_l < chunk_size; ++pad_l) {
                                dst_chunk.state[pad_l] = ParticleState::INVALID;
                                dst_chunk.pos_x[pad_l] = 1e50;
                                dst_chunk.pos_y[pad_l] = 1e50;
                                dst_chunk.pos_z[pad_l] = 1e50;
                                dst_chunk.id[pad_l]    = std::numeric_limits<ParticleID>::max();
                                dst_chunk.mass[pad_l]  = 1.0;
                            }
                        }
                    }
                }
            });

            std::swap(data, tmp);
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
        auto get_field_ptr(this auto&& self, const size_t i) {
            // locate data
            const auto [chunk_idx, lane_idx] = self.locate(i);
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
        alignas(64) packed::value_type idx_arr[packed::size()]{};

        template <ParticleField Read, ParticleField Write>
        [[nodiscard]] auto access_particle(this auto&& self, size_t chunk_idx, size_t lane_idx) {
            constexpr bool is_const = std::is_const_v<std::remove_reference_t<decltype(self)>>;

            static_assert(!(is_const && Write != ParticleField::none),
                          "[APRIL] Cannot request write permissions (WriteMask != none) on a const Container. "
                          "Either drop the write mask or ensure the container is mutable.");

            particle::internal::ParticleSource<Read, Write, ParticleAttributes> src;
            constexpr auto Mask = Read | Write;

            auto& chunk = self.ptr_chunks[chunk_idx];

            if constexpr (particle::internal::has_field_v<Mask, ParticleField::force>)
                src.force = math::Vec3Ptr{&chunk.frc_x[lane_idx], &chunk.frc_y[lane_idx], &chunk.frc_z[lane_idx]};
            if constexpr (particle::internal::has_field_v<Mask, ParticleField::position>)
                src.position = math::Vec3Ptr{&chunk.pos_x[lane_idx], &chunk.pos_y[lane_idx], &chunk.pos_z[lane_idx]};
            if constexpr (particle::internal::has_field_v<Mask, ParticleField::velocity>)
                src.velocity = math::Vec3Ptr{&chunk.vel_x[lane_idx], &chunk.vel_y[lane_idx], &chunk.vel_z[lane_idx]};
            if constexpr (particle::internal::has_field_v<Mask, ParticleField::old_position>)
                src.old_position = math::Vec3Ptr{&chunk.old_x[lane_idx], &chunk.old_y[lane_idx], &chunk.old_z[lane_idx]};
            if constexpr (particle::internal::has_field_v<Mask, ParticleField::mass>)
                src.mass = &chunk.mass[lane_idx];
            if constexpr (particle::internal::has_field_v<Mask, ParticleField::state>)
                src.state = &chunk.state[lane_idx];
            if constexpr (particle::internal::has_field_v<Mask, ParticleField::type>)
                src.type = &chunk.type[lane_idx];
            if constexpr (particle::internal::has_field_v<Mask, ParticleField::id>)
                src.id = &chunk.id[lane_idx];
            if constexpr (particle::internal::has_field_v<Mask, ParticleField::attributes>)
                src.attributes = &chunk.attributes[lane_idx];

            return src;
        }

        template <ParallelPolicy P, exec::ExecutionMode V, bool is_const, exec::IsKernel Kernel>
        void iterate_range(this auto&& self, Kernel&& kernel, const size_t start, const size_t end) {
            // route scalar/vector execution
            auto process_sub_range = [&](const size_t r_start, const size_t r_end) AP_FORCE_INLINE {
                if constexpr (V == exec::ExecutionMode::Scalar) {
                    self.template iterate_range_scalar<P, is_const>(kernel, r_start, r_end);
                } else if constexpr (V == exec::ExecutionMode::Vector || V == exec::ExecutionMode::Hybrid) {
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
        AP_FORCE_INLINE void iterate_range_scalar(this auto&& self, Kernel&& kernel, const size_t start, const size_t end) {
            using K = std::remove_cvref_t<Kernel>;
            if (start >= end) return;

            const auto [start_chunk, start_idx] = self.locate(start);
            const auto [end_chunk, end_idx] = self.locate(end);

            auto* AP_RESTRICT chunks = self.ptr_chunks;
            auto exec_scalar = [&](size_t c, size_t i) AP_FORCE_INLINE {
                const size_t physical_idx = (c << chunk_shift) | i;
                if constexpr (is_const) kernel(physical_idx, self.template view<K::Read>(c, i));
                else kernel(physical_idx, self.template at<K::Read, K::Write>(c, i));
            };

            if (start_chunk == end_chunk) {
                AP_PREFETCH(chunks + start_chunk);
                for (size_t i = start_idx; i < end_idx; ++i) {
                    exec_scalar(start_chunk, i);
                }
            }
            else {
                AP_PREFETCH(chunks + start_chunk);
                if (start_chunk + 1 < end_chunk || end_idx > 0) AP_PREFETCH(chunks + start_chunk + 1);

                // head
                for (size_t i = start_idx; i < self.chunk_size; ++i) {
                    exec_scalar(start_chunk, i);
                }

                // body
                for (size_t c = start_chunk + 1; c < end_chunk; ++c) {
                    AP_PREFETCH(chunks + c + 1);
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
        AP_FORCE_INLINE void iterate_range_vector(this auto&& self, Kernel&& kernel, const size_t start, const size_t end) {
            using K = std::remove_cvref_t<Kernel>;
            if (start >= end) return;

            // get the first and last chunks touched in the iteration
            const auto [start_chunk, start_idx] = self.locate(start);
            const auto [end_chunk, end_idx] = self.locate(end);

            constexpr size_t simd_width = packed::size();
            auto* AP_RESTRICT chunks = self.ptr_chunks;

            // Align backwards to the nearest vector boundary for head masking
            const size_t head_offset = start_idx % simd_width;
            const size_t aligned_start_idx = start_idx - head_offset;

            // kernel expects the logical base index of the vector register
            const auto lane_indices = packed::load_aligned(self.idx_arr);

            auto exec_vector = [&](size_t c, size_t i) AP_FORCE_INLINE {
                const size_t physical_idx = (c << chunk_shift) | i;
                if constexpr (is_const) kernel(physical_idx, self.template view_packed<K::Read>(c, i));
                else kernel(physical_idx, self.template at_packed<K::Read, K::Write>(c, i));
            };

            auto exec_vector_masked = [&](size_t c, size_t i, auto mask) AP_FORCE_INLINE {
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
                AP_PREFETCH(chunks + start_chunk);

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
                AP_PREFETCH(chunks + start_chunk);
                if (start_chunk + 1 < end_chunk || end_idx > 0) AP_PREFETCH(chunks + start_chunk + 1);

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
                    AP_PREFETCH(chunks + c + 1);
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


