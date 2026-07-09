#pragma once


#include <array>
#include <cstddef>
#include <bit>

#include "april/base/types.hpp"
#include "april/particle/properties.hpp"
#include "april/exec/info.hpp"
#include "april/particle/attributes.hpp"


namespace april::container::layout {
    template <particle::IsParticleAttributes Attributes, size_t Size>
    struct alignas(64) ParticleChunk {
        // Enforce alignment requirements (Size 8 * double 8 bytes = 64 bytes = 1 AVX-512 Register)
        static_assert(std::has_single_bit(Size),
                      "Chunk Size must be a Power of 2 (e.g., 8, 16) for bitwise indexing optimizations.");
        static_assert(Size * sizeof(vec3::type) >= exec::CACHE_LINE_SIZE,
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

        // Attributes
        alignas(64) std::array<Attributes, Size> attributes;

        void insert_particle(size_t l_idx, const particle::ParticleRecord<Attributes> & p) {
            pos_x[l_idx] = p.position.x;
            pos_y[l_idx] = p.position.y;
            pos_z[l_idx] = p.position.z;

            vel_x[l_idx] = p.velocity.x;
            vel_y[l_idx] = p.velocity.y;
            vel_z[l_idx] = p.velocity.z;

            frc_x[l_idx] = p.force.x;
            frc_y[l_idx] = p.force.y;
            frc_z[l_idx] = p.force.z;

            old_x[l_idx] = p.old_position.x;
            old_y[l_idx] = p.old_position.y;
            old_z[l_idx] = p.old_position.z;

            mass[l_idx] = p.mass;
            state[l_idx] = p.state;
            type[l_idx] = p.type;
            id[l_idx] = p.id;

            attributes[l_idx] = p.attributes;
        }


        void copy_from(size_t dst_l, size_t src_l, const ParticleChunk & src_chunk) {
            pos_x[dst_l]      = src_chunk.pos_x[src_l];
            pos_y[dst_l]      = src_chunk.pos_y[src_l];
            pos_z[dst_l]      = src_chunk.pos_z[src_l];

            vel_x[dst_l]      = src_chunk.vel_x[src_l];
            vel_y[dst_l]      = src_chunk.vel_y[src_l];
            vel_z[dst_l]      = src_chunk.vel_z[src_l];

            frc_x[dst_l]      = src_chunk.frc_x[src_l];
            frc_y[dst_l]      = src_chunk.frc_y[src_l];
            frc_z[dst_l]      = src_chunk.frc_z[src_l];

            old_x[dst_l]      = src_chunk.old_x[src_l];
            old_y[dst_l]      = src_chunk.old_y[src_l];
            old_z[dst_l]      = src_chunk.old_z[src_l];

            mass[dst_l]       = src_chunk.mass[src_l];
            state[dst_l]      = src_chunk.state[src_l];
            type[dst_l]       = src_chunk.type[src_l];
            id[dst_l]         = src_chunk.id[src_l];

            attributes[dst_l] = src_chunk.attributes[src_l];
        }
    };
}