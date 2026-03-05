#pragma once


#include <array>
#include <cstddef>
#include <bit>

#include "april/containers/container.hpp"
#include "april/base/types.hpp"
#include "april/particle/particle_types.hpp"
#include "april/exec/policy.hpp"


namespace april::container::layout {
    template <particle::IsParticleAttributes Attributes, size_t Size>
    struct alignas(64) ParticleChunk {
        // Enforce alignment requirements (Size 8 * double 8 bytes = 64 bytes = 1 AVX-512 Register)
        static_assert(std::has_single_bit(Size),
                      "Chunk Size must be a Power of 2 (e.g., 8, 16) for bitwise indexing optimizations.");
        static_assert(Size * sizeof(vec3::type) >= 64,
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

        // Old Position
        alignas(64) std::array<vec3::type, Size> old_x;
        alignas(64) std::array<vec3::type, Size> old_y;
        alignas(64) std::array<vec3::type, Size> old_z;

        // Scalars
        alignas(64) std::array<double, Size> mass;
        alignas(64) std::array<ParticleState, Size> state;
        alignas(64) std::array<ParticleType, Size> type;
        alignas(64) std::array<ParticleID, Size> id;
        alignas(64) std::array<Attributes, Size> attributes;
    };


    template <typename Config, particle::IsParticleAttributes Attributes, size_t ChunkSize>
    class AoSoA : public Container<Config, Attributes> {
    public:
        static constexpr size_t chunk_size = ChunkSize;
        using ChunkT = ParticleChunk<Attributes, chunk_size>;

        using Base = Container<Config, Attributes>;
        using Base::force_schema;
        using Base::Base; // Inherit constructors
        friend Base;

        // make inherited accessors explicit.
        // Otherwise the compiler cant find them due to existing overrides in this class
        using Base::view;
        using Base::at;
        using Base::access_particle;


        // ACCESSORS (chunk based)
        template <ParticleField Read, ParticleField Write>
        [[nodiscard]] auto at(this auto&& self, size_t chunk_idx, size_t lane_idx) {
            return particle::internal::ScalarParticleRef<Read, Write, Attributes>{
                self.template access_particle<Read, Write>(chunk_idx, lane_idx)
            };
        }

        template <ParticleField Read>
        [[nodiscard]] auto view(this const auto& self, size_t chunk_idx, size_t lane_idx) {
            return particle::internal::ScalarParticleRef<Read, ParticleField::none, Attributes>{
                self.template access_particle<Read, ParticleField::none>(chunk_idx, lane_idx)
            };
        }

        template <ParticleField Read, ParticleField Write>
        [[nodiscard]] auto at_packed(this auto&& self, size_t chunk_idx, size_t lane_idx) {
            return particle::internal::PackedParticleRef<Read, Write, Attributes>{
                self.template access_particle<Read, Write>(chunk_idx, lane_idx)
            };
        }

        template <ParticleField Read>
        [[nodiscard]] auto view_packed(this const auto& self, size_t chunk_idx, size_t lane_idx) {
            return particle::internal::PackedParticleRef<Read, ParticleField::none, Attributes>{
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

        void update_cache() {
            ptr_chunks = data.data();
        }

        void build_storage(const std::vector<particle::ParticleRecord<Attributes>>& particles) {
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

                chunk.mass[l_idx] = p.mass;
                chunk.state[l_idx] = p.state;
                chunk.type[l_idx] = p.type;
                chunk.id[l_idx] = p.id;
                chunk.attributes[l_idx] = p.attributes;

                // Map ID
                id_to_index_map[static_cast<size_t>(p.id)] = i;
            }

            for (size_t i = n_particles; i < particle_capacity; ++i) {
                const auto [c_idx, l_idx] = locate(i);
                data[c_idx].state[l_idx] = ParticleState::INVALID;
                data[c_idx].id[l_idx] = std::numeric_limits<ParticleID>::max();
                data[c_idx].type[l_idx] = std::numeric_limits<ParticleType>::max();
            }
        }

        void reorder_storage(const std::vector<std::vector<size_t>>& bins, const bool sentinel_pad = true) {
            tmp.clear();
            bin_starts.clear();
            bin_sizes.clear();

            // allocate space, calculate chunk index of first chunk in each bin
            bin_starts.reserve(bins.size() + 1);

            size_t n_chunks = 0;
            for (const auto& bin : bins) {
                bin_starts.push_back(n_chunks * chunk_size);
                bin_sizes.push_back(bin.size());
                if (bin.empty()) continue;
                const size_t bin_size = bin.size();
                n_chunks += (bin_size + chunk_size - 1) / chunk_size;
            }

            particle_capacity = n_chunks * chunk_size;
            bin_starts.push_back(particle_capacity);

            tmp.resize(n_chunks);
            id_to_index_map.assign(particle_capacity, ID_NOT_FOUND);

            // traversal cursors
            size_t dst_c = 0; // destination chunk index
            size_t dst_l = 0; // destination lane index

            for (const auto& bin : bins) {
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

                    dst_chunk.mass[dst_l] = src_chunk.mass[src_l];
                    dst_chunk.state[dst_l] = src_chunk.state[src_l];
                    dst_chunk.type[dst_l] = src_chunk.type[src_l];
                    dst_chunk.id[dst_l] = src_chunk.id[src_l];
                    dst_chunk.attributes[dst_l] = src_chunk.attributes[src_l];

                    const size_t new_physical_idx = dst_c * chunk_size + dst_l;
                    const ParticleID id = dst_chunk.id[dst_l];
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
                        dst_chunk.state[dst_l] = ParticleState::INVALID;

                        // move far away to be safe against distance checks
                        dst_chunk.pos_x[dst_l] = 1e50;
                        dst_chunk.pos_y[dst_l] = 1e50;
                        dst_chunk.pos_z[dst_l] = 1e50;

                        // set ID to max to avoid map lookups
                        dst_chunk.id[dst_l] = std::numeric_limits<ParticleID>::max();

                        // make mass safe
                        dst_chunk.mass[dst_l] = 1.0;

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
        template <ParticleField Read, ParticleField Write>
        [[nodiscard]] auto access_particle(this auto&& self, size_t chunk_idx, size_t lane_idx) {
            constexpr bool is_const = std::is_const_v<std::remove_reference_t<decltype(self)>>;

            static_assert(!(is_const && Write != ParticleField::none),
                          "APRIL ERROR: Cannot request write permissions (WriteMask != none) on a const Container. "
                          "Either drop the write mask or ensure the container is mutable.");

            particle::internal::ParticleSource<Read, Write, Attributes> src;
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

        template <ParallelPolicy P, exec::internal::ExecutionMode V, bool is_const, exec::IsKernel Kernel>
        void iterate_range(this auto&& self, Kernel&& kernel, const size_t start, const size_t end) {
            if constexpr (V == exec::internal::ExecutionMode::Scalar) {
                self.template iterate_range_scalar<P, is_const>(std::forward<Kernel>(kernel), start, end);
            }
            else if constexpr (V == exec::internal::ExecutionMode::Vector) {
                self.template iterate_range_vector<P, is_const>(std::forward<Kernel>(kernel), start, end);
            }
            else if constexpr (V == exec::internal::ExecutionMode::Hybrid) {
                self.template iterate_range_hybrid<P, is_const>(std::forward<Kernel>(kernel), start, end);
            }
        }

        template <ParallelPolicy P, bool is_const, exec::IsKernel Kernel>
        AP_FORCE_INLINE void iterate_range_scalar(this auto&& self, Kernel&& kernel, const size_t start,
                                                  const size_t end) {
            using K = std::remove_cvref_t<Kernel>;
            if (start >= end) return;

            const auto [start_chunk, start_idx] = self.locate(start);
            const auto [end_chunk, end_idx] = self.locate(end);

            size_t curr_idx = start;
            auto* AP_RESTRICT chunks = self.ptr_chunks;
            // constexpr size_t simd_width = packed::size();

            auto exec_scalar = [&](size_t c, size_t i) {
                if constexpr (is_const) kernel(curr_idx++, self.template view<K::Read>(c, i));
                else kernel(curr_idx++, self.template at<K::Read, K::Write>(c, i));
            };

            if (start_chunk == end_chunk) {
                AP_PREFETCH(chunks + start_chunk);
                for (size_t i = start_idx; i < end_idx; ++i) {
                    exec_scalar(start_chunk, i);
                }
            }
            else {
                AP_PREFETCH(chunks + start_chunk);
                if (start_chunk + 1 < end_chunk || end_idx > 0)
                    AP_PREFETCH(chunks + start_chunk + 1);

                for (size_t i = start_idx; i < self.chunk_size; ++i) {
                    exec_scalar(start_chunk, i);
                }

                for (size_t c = start_chunk + 1; c < end_chunk; ++c) {
                    AP_PREFETCH(chunks + c + 1);
                    for (size_t i = 0; i < self.chunk_size; ++i) exec_scalar(c, i);
                }

                if (end_idx > 0) {
                    for (size_t i = 0; i < end_idx; ++i) {
                        exec_scalar(end_chunk, i);
                    }
                }
            }
        }

        template <ParallelPolicy P, bool is_const, exec::IsKernel Kernel>
        AP_FORCE_INLINE void iterate_range_vector(this auto&& self, Kernel&& kernel, const size_t start,
                                                  const size_t end) {
            using K = std::remove_cvref_t<Kernel>;
            if (start >= end) return;

            const auto [start_chunk, start_idx] = self.locate(start);
            const auto [end_chunk, end_idx] = self.locate(end);

            size_t curr_idx = start;
            auto* AP_RESTRICT chunks = self.ptr_chunks;
            constexpr size_t simd_width = packed::size();

            auto exec_vector = [&](size_t c, size_t i) {
                if constexpr (is_const) kernel(curr_idx, self.template view_packed<K::Read>(c, i));
                else kernel(curr_idx, self.template at_packed<K::Read, K::Write>(c, i));
                curr_idx += simd_width;
            };

            AP_ASSERT(start_idx % simd_width == 0, "Vector execution requires an aligned start index");

            if (start_chunk == end_chunk) {
                AP_PREFETCH(chunks + start_chunk);
                for (size_t i = start_idx; i < end_idx; i += simd_width) {
                    exec_vector(start_chunk, i);
                }
            }
            else {
                AP_PREFETCH(chunks + start_chunk);
                if (start_chunk + 1 < end_chunk || end_idx > 0)
                    AP_PREFETCH(chunks + start_chunk + 1);

                for (size_t i = start_idx; i < self.chunk_size; i += simd_width) {
                    exec_vector(start_chunk, i);
                }

                for (size_t c = start_chunk + 1; c < end_chunk; ++c) {
                    AP_PREFETCH(chunks + c + 1);
                    for (size_t i = 0; i < self.chunk_size; i += simd_width) {
                        exec_vector(c, i);
                    }
                }

                if (end_idx > 0) {
                    for (size_t i = 0; i < end_idx; i += simd_width) {
                        exec_vector(end_chunk, i);
                    }
                }
            }
        }

        template <ParallelPolicy P, bool is_const, exec::IsKernel Kernel>
        AP_FORCE_INLINE void iterate_range_hybrid(this auto&& self, Kernel&& kernel, const size_t start,
                                                  const size_t end) {
            using K = std::remove_cvref_t<Kernel>;
            if (start >= end) return;

            const auto [start_chunk, start_idx] = self.locate(start);
            const auto [end_chunk, end_idx] = self.locate(end);

            size_t curr_idx = start;
            auto* AP_RESTRICT chunks = self.ptr_chunks;
            constexpr size_t simd_width = packed::size();

            // scoped macros to make code DRY but ensure inlining
            #define EXEC_SCALAR(c, i) \
                if constexpr (is_const) kernel(curr_idx++, self.template view<K::Read>(c, i)); \
                else kernel(curr_idx++, self.template at<K::Read, K::Write>(c, i))

            #define EXEC_VECTOR(c, i) \
                if constexpr (is_const) kernel(curr_idx, self.template view_packed<K::Read>(c, i)); \
                else { \
                    auto packed = self.template at_packed<K::Read, K::Write>(c, i); \
                    auto buffer = packed.load_buffer(); \
                    kernel(curr_idx, buffer.to_view()); \
                    buffer.update_into(packed); \
                } \
                curr_idx += simd_width

            // iterate a single chunk
            if (start_chunk == end_chunk) {
                AP_PREFETCH(chunks + start_chunk);

                const size_t head_end = (start_idx % simd_width == 0)
                                            ? start_idx
                                            : std::min(end_idx, start_idx + (simd_width - (start_idx % simd_width)));
                for (size_t i = start_idx; i < head_end; ++i) { EXEC_SCALAR(start_chunk, i); }

                const size_t body_start = head_end;
                const size_t body_end = body_start + ((end_idx - body_start) / simd_width) * simd_width;
                for (size_t i = body_start; i < body_end; i += simd_width) { EXEC_VECTOR(start_chunk, i); }

                for (size_t i = body_end; i < end_idx; ++i) { EXEC_SCALAR(start_chunk, i); }
            }
            // iterate multiple chunks
            else {
                AP_PREFETCH(chunks + start_chunk);
                if (start_chunk + 1 < end_chunk || end_idx > 0)
                    AP_PREFETCH(chunks + start_chunk + 1);

                // handle head chunk
                const size_t head_end = (start_idx % simd_width == 0)
                                            ? start_idx
                                            : (start_idx + (simd_width - (start_idx % simd_width)));

                for (size_t i = start_idx; i < head_end; ++i) { EXEC_SCALAR(start_chunk, i); }
                for (size_t i = head_end; i < self.chunk_size; i += simd_width) { EXEC_VECTOR(start_chunk, i); }

                // body
                for (size_t c = start_chunk + 1; c < end_chunk; ++c) {
                    AP_PREFETCH(chunks + c + 1);
                    // Guaranteed aligned inside full body chunks
                    for (size_t i = 0; i < self.chunk_size; i += simd_width) { EXEC_VECTOR(c, i); }
                }

                // handle tail
                if (end_idx > 0) {
                    const size_t tail_body_end = (end_idx / simd_width) * simd_width;
                    for (size_t i = 0; i < tail_body_end; i += simd_width) { EXEC_VECTOR(end_chunk, i); }
                    for (size_t i = tail_body_end; i < end_idx; ++i) { EXEC_SCALAR(end_chunk, i); }
                }
            }

            #undef EXEC_SCALAR
            #undef EXEC_VECTOR
        }
    };
}
