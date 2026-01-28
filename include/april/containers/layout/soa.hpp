#pragma once
#include <vector>
#include "april/containers/container.hpp"
#include "april/containers/batching/common.hpp"
#include "april/particle/particle.hpp"

namespace april::container::layout {

    template<typename UserData>
    struct SoAStorage {
        alignas(64) std::vector<vec3::type> pos_x, pos_y, pos_z;
        alignas(64) std::vector<vec3::type> vel_x, vel_y, vel_z;
        alignas(64) std::vector<vec3::type> frc_x, frc_y, frc_z;
        alignas(64) std::vector<vec3::type> old_x, old_y, old_z;

        alignas(64) std::vector<double> mass;
        alignas(64) std::vector<env::ParticleState> state;
        alignas(64) std::vector<env::ParticleType> type;
        alignas(64) std::vector<env::ParticleID> id;
        alignas(64) std::vector<UserData> user_data;

        vec3::type * AP_RESTRICT ptr_pos_x = nullptr;
        vec3::type * AP_RESTRICT ptr_pos_y = nullptr;
        vec3::type * AP_RESTRICT ptr_pos_z = nullptr;

        // Velocities
        vec3::type * AP_RESTRICT ptr_vel_x = nullptr;
        vec3::type * AP_RESTRICT ptr_vel_y = nullptr;
        vec3::type * AP_RESTRICT ptr_vel_z = nullptr;

        // Forces
        vec3::type * AP_RESTRICT ptr_frc_x = nullptr;
        vec3::type * AP_RESTRICT ptr_frc_y = nullptr;
        vec3::type * AP_RESTRICT ptr_frc_z = nullptr;

        // Old Position
        vec3::type * AP_RESTRICT ptr_old_x = nullptr;
        vec3::type * AP_RESTRICT ptr_old_y = nullptr;
        vec3::type * AP_RESTRICT ptr_old_z = nullptr;

        // Scalars
        double * AP_RESTRICT ptr_mass = nullptr;
        env::ParticleState * AP_RESTRICT ptr_state = nullptr;
        env::ParticleType  * AP_RESTRICT ptr_type  = nullptr;
        env::ParticleID    * AP_RESTRICT ptr_id    = nullptr;
        UserData * AP_RESTRICT ptr_user_data = nullptr;

        void update_pointer_cache() {
            ptr_pos_x = pos_x.data(); ptr_pos_y = pos_y.data(); ptr_pos_z = pos_z.data();
            ptr_vel_x = vel_x.data(); ptr_vel_y = vel_y.data(); ptr_vel_z = vel_z.data();
            ptr_frc_x = frc_x.data(); ptr_frc_y = frc_y.data(); ptr_frc_z = frc_z.data();
            ptr_old_x = old_x.data(); ptr_old_y = old_y.data(); ptr_old_z = old_z.data();

            ptr_mass = mass.data();
            ptr_state = state.data();
            ptr_type = type.data();
            ptr_id = id.data();
            ptr_user_data = user_data.data();
        }

        void resize(const size_t n) {
            pos_x.resize(n); pos_y.resize(n); pos_z.resize(n);
            vel_x.resize(n); vel_y.resize(n); vel_z.resize(n);
            frc_x.resize(n); frc_y.resize(n); frc_z.resize(n);
            old_x.resize(n); old_y.resize(n); old_z.resize(n);
            mass.resize(n); state.resize(n); type.resize(n); id.resize(n);
            user_data.resize(n);
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
            user_data[dest_i] = src.user_data[src_i];
        }

        void swap(const size_t i, const size_t j) {
            if (i == j) return;

            std::swap(pos_x[i], pos_x[j]); std::swap(pos_y[i], pos_y[j]); std::swap(pos_z[i], pos_z[j]);
            std::swap(vel_x[i], vel_x[j]); std::swap(vel_y[i], vel_y[j]); std::swap(vel_z[i], vel_z[j]);
            std::swap(frc_x[i], frc_x[j]); std::swap(frc_y[i], frc_y[j]); std::swap(frc_z[i], frc_z[j]);
            std::swap(old_x[i], old_x[j]); std::swap(old_y[i], old_y[j]); std::swap(old_z[i], old_z[j]);

            // Swap Scalars
            std::swap(mass[i], mass[j]);
            std::swap(state[i], state[j]);
            std::swap(type[i], type[j]);
            std::swap(id[i], id[j]);
            std::swap(user_data[i], user_data[j]);
        }
    };



    template<typename Config, env::IsUserData U>
    class SoA : public Container<Config, U> {
    public:
        using Base = Container<Config, U>;
        using Base::force_schema;
        using Base::Base;
        friend Base;

        SoA(const Config & config, const internal::ContainerCreateInfo & info)
            : Base(config, info)
        {
            // precompute topology batches (id based batches)
            for (size_t i = 0; i < force_schema.interactions.size(); ++i) {
                const auto& prop = force_schema.interactions[i];

                if (!prop.used_by_ids.empty() && prop.is_active) {
                    TopologyBatch batch;
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
        void iterate_range(this auto&& self, Kernel && kernel, const size_t start, const size_t end) {
            for (size_t i = start; i < end; i++) {
                if constexpr (is_const) {
                    kernel(i, self.template view<M>(i));
                } else {
                    kernel(i, self.template at<M>(i));
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
            return static_cast<env::ParticleID>(id_to_index_map.size());
        }
        [[nodiscard]] bool index_is_valid(const size_t index) const {
            return index < particle_count();
        }
        [[nodiscard]] bool contains_id(const env::ParticleID id) const {
            return id <= max_id();
        }


        // QUERIES
        [[nodiscard]] size_t capacity() const {
            return particle_count();
        }
        [[nodiscard]] size_t particle_count() const {
            return data.pos_x.size();
        }

    protected:
        SoAStorage<U> tmp;
        SoAStorage<U> data;
        std::vector<size_t> bin_starts; // first particle index of each bin
        std::vector<size_t> bin_sizes; // number of particles in each bin
        std::vector<uint32_t> id_to_index_map;

        // explode AoS input into SoA vectors
        void build_storage(const std::vector<env::internal::ParticleRecord<U>>& particles) {
            size_t n = particles.size();
            data.resize(n);
            id_to_index_map.resize(n);

            bin_starts.clear();
            bin_sizes.clear();
            bin_starts.push_back(0);
            bin_sizes.push_back(particles.size());

            for (size_t i = 0; i < n; ++i) {
                const auto& p = particles[i];

                // Vectors
                data.pos_x[i] = p.position.x;        data.pos_y[i] = p.position.y;        data.pos_z[i] = p.position.z;
                data.vel_x[i] = p.velocity.x;        data.vel_y[i] = p.velocity.y;        data.vel_z[i] = p.velocity.z;
                data.frc_x[i] = p.force.x;           data.frc_y[i] = p.force.y;           data.frc_z[i] = p.force.z;
                data.old_x[i] = p.old_position.x;    data.old_y[i] = p.old_position.y;    data.old_z[i] = p.old_position.z;

                // Scalars
                data.mass[i] = p.mass;
                data.state[i] = p.state;
                data.type[i] = p.type;
                data.id[i] = p.id;
                data.user_data[i] = p.user_data;

                // ID Map
                id_to_index_map[static_cast<size_t>(p.id)] = i;
            }

            tmp.resize(n);
        }

        void reorder_storage(const std::vector<std::vector<size_t>>& bins) {
            bin_starts.clear();
            bin_sizes.clear();

            // scatter particles into bins & update bin meta data
            size_t current_idx = 0;
            size_t current_offset = 0;

            for (const auto& bin : bins) {
                bin_starts.push_back(current_offset);
                bin_sizes.push_back(bin.size());
                current_offset += bin.size();

                for (size_t old_idx : bin) {
                    tmp.copy_from(current_idx, data, old_idx);
                    current_idx++;
                }
            }

            // swap buffers and update pointer caches
            std::swap(data, tmp);
            data.update_pointer_cache();
            tmp.update_pointer_cache();

            // rebuild iD Map
            for (size_t i = 0; i < data.id.size(); i++) {
                const auto id = data.id[i];
                if (static_cast<size_t>(id) >= id_to_index_map.size()) {
                    id_to_index_map.resize(static_cast<size_t>(id) + 1, std::numeric_limits<uint32_t>::max());
                }
                id_to_index_map[static_cast<size_t>(id)] = static_cast<uint32_t>(i);
            }
        }

        [[nodiscard]] std::pair<size_t, size_t> get_physical_bin_range(const size_t type) const {
            const size_t start = bin_starts[type];
            return {start, start + bin_sizes[type]};
        }

        template<env::Field F>
        auto get_field_ptr(this auto&& self, size_t i) {
            if constexpr (F == env::Field::position)
                return utils::Vec3Ptr { self.data.ptr_pos_x + i, self.data.ptr_pos_y + i, self.data.ptr_pos_z + i };
            else if constexpr (F == env::Field::velocity)
                return utils::Vec3Ptr { self.data.ptr_vel_x + i, self.data.ptr_vel_y + i, self.data.ptr_vel_z + i };
            else if constexpr (F == env::Field::force)
                return utils::Vec3Ptr { self.data.ptr_frc_x + i, self.data.ptr_frc_y + i, self.data.ptr_frc_z + i };
            else if constexpr (F == env::Field::old_position)
                return utils::Vec3Ptr { self.data.ptr_old_x + i, self.data.ptr_old_y + i, self.data.ptr_old_z + i };

            else if constexpr (F == env::Field::mass)      return self.data.ptr_mass + i;
            else if constexpr (F == env::Field::state)     return self.data.ptr_state + i;
            else if constexpr (F == env::Field::type)      return self.data.ptr_type + i;
            else if constexpr (F == env::Field::id)        return self.data.ptr_id + i;
            else if constexpr (F == env::Field::user_data) return self.data.ptr_user_data + i;
        }

    private:
        std::vector<TopologyBatch> topology_batches;
    };
}
