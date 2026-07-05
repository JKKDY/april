#pragma once
#include <memory>
#include <vector>
#include "april/particle/particle.hpp"


namespace april::container::layout {
    template<particle::IsParticleAttributes Attributes>
    struct SoAStorage {
        alignas(64) std::vector<vec3::type> pos_x, pos_y, pos_z;
        alignas(64) std::vector<vec3::type> vel_x, vel_y, vel_z;
        alignas(64) std::vector<vec3::type> frc_x, frc_y, frc_z;
        alignas(64) std::vector<vec3::type> old_x, old_y, old_z;

        alignas(64) std::vector<double> mass;
        alignas(64) std::vector<ParticleState> state;
        alignas(64) std::vector<ParticleType> type;
        alignas(64) std::vector<ParticleID> id;
        alignas(64) std::vector<Attributes> attributes;

        vec3::type * APRIL_RESTRICT ptr_pos_x = nullptr;
        vec3::type * APRIL_RESTRICT ptr_pos_y = nullptr;
        vec3::type * APRIL_RESTRICT ptr_pos_z = nullptr;

        // Velocities
        vec3::type * APRIL_RESTRICT ptr_vel_x = nullptr;
        vec3::type * APRIL_RESTRICT ptr_vel_y = nullptr;
        vec3::type * APRIL_RESTRICT ptr_vel_z = nullptr;

        // Forces
        vec3::type * APRIL_RESTRICT ptr_frc_x = nullptr;
        vec3::type * APRIL_RESTRICT ptr_frc_y = nullptr;
        vec3::type * APRIL_RESTRICT ptr_frc_z = nullptr;

        // Old Position
        vec3::type * APRIL_RESTRICT ptr_old_x = nullptr;
        vec3::type * APRIL_RESTRICT ptr_old_y = nullptr;
        vec3::type * APRIL_RESTRICT ptr_old_z = nullptr;

        // Scalars
        double * APRIL_RESTRICT ptr_mass = nullptr;
        ParticleState * APRIL_RESTRICT ptr_state = nullptr;
        ParticleType  * APRIL_RESTRICT ptr_type  = nullptr;
        ParticleID    * APRIL_RESTRICT ptr_id    = nullptr;
        Attributes * APRIL_RESTRICT ptr_attributes = nullptr;

        size_t capacity{};
        size_t size{};

        void update_pointer_cache() {
            ptr_pos_x = pos_x.data(); ptr_pos_y = pos_y.data(); ptr_pos_z = pos_z.data();
            ptr_vel_x = vel_x.data(); ptr_vel_y = vel_y.data(); ptr_vel_z = vel_z.data();
            ptr_frc_x = frc_x.data(); ptr_frc_y = frc_y.data(); ptr_frc_z = frc_z.data();
            ptr_old_x = old_x.data(); ptr_old_y = old_y.data(); ptr_old_z = old_z.data();

            ptr_mass = mass.data();
            ptr_state = state.data();
            ptr_type = type.data();
            ptr_id = id.data();
            ptr_attributes = attributes.data();
        }

        void insert_particle(size_t idx, const particle::ParticleRecord<Attributes> & p) {
            pos_x[idx] = p.position.x;
            pos_y[idx] = p.position.y;
            pos_z[idx] = p.position.z;

            vel_x[idx] = p.velocity.x;
            vel_y[idx] = p.velocity.y;
            vel_z[idx] = p.velocity.z;

            frc_x[idx] = p.force.x;
            frc_y[idx] = p.force.y;
            frc_z[idx] = p.force.z;

            old_x[idx] = p.old_position.x;
            old_y[idx] = p.old_position.y;
            old_z[idx] = p.old_position.z;

            mass[idx] = p.mass;
            state[idx] = p.state;
            type[idx] = p.type;
            id[idx] = p.id;

            attributes[idx] = p.attributes;
        }

        void resize(const size_t n) {
            // pad with packed size then round up to next multiple of packed size
            capacity = ((n + packed::size() - 1) / packed::size()) * packed::size() + packed::size();
            size = n;

            pos_x.resize(capacity); pos_y.resize(capacity); pos_z.resize(capacity);
            vel_x.resize(capacity); vel_y.resize(capacity); vel_z.resize(capacity);
            frc_x.resize(capacity); frc_y.resize(capacity); frc_z.resize(capacity);
            old_x.resize(capacity); old_y.resize(capacity); old_z.resize(capacity);

            mass.resize(capacity);
            state.resize(capacity);
            type.resize(capacity);
            id.resize(capacity);
            attributes.resize(capacity);

            capacity -= packed::size();
            update_pointer_cache();
        }

        // Copy particle data from source index to this index
        // Used for rebuilding/sorting storage
        void copy_from(const size_t dest_i, const SoAStorage& src, const size_t src_i) {
            pos_x[dest_i] = src.pos_x[src_i]; pos_y[dest_i] = src.pos_y[src_i]; pos_z[dest_i] = src.pos_z[src_i];
            vel_x[dest_i] = src.vel_x[src_i]; vel_y[dest_i] = src.vel_y[src_i]; vel_z[dest_i] = src.vel_z[src_i];
            frc_x[dest_i] = src.frc_x[src_i]; frc_y[dest_i] = src.frc_y[src_i]; frc_z[dest_i] = src.frc_z[src_i];
            old_x[dest_i] = src.old_x[src_i]; old_y[dest_i] = src.old_y[src_i]; old_z[dest_i] = src.old_z[src_i];

            mass[dest_i]      = src.mass[src_i];
            state[dest_i]     = src.state[src_i];
            type[dest_i]      = src.type[src_i];
            id[dest_i]        = src.id[src_i];
            attributes[dest_i] = src.attributes[src_i];
        }
    };
}