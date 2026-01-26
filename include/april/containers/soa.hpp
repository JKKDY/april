#pragma once
#include <vector>
#include "april/containers/container.hpp"
#include "april/containers/batching.hpp"
#include "april/particle/particle.hpp"

namespace april::container {

    template<typename Config, env::IsUserData U>
    class SoAContainer : public Container<Config, U> {
    public:
        using Base = Container<Config, U>;
        using Base::force_schema;
        using Base::Base;
        friend Base;

        SoAContainer(const Config & config, const internal::ContainerCreateInfo & info)
            : Base(config, info)
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
                if (!self.index_is_valid(i)) continue;
                auto raw_state = self.data.state[i];

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
        struct Storage {
            alignas(64) std::vector<vec3::type> pos_x, pos_y, pos_z;
            alignas(64) std::vector<vec3::type> vel_x, vel_y, vel_z;
            alignas(64) std::vector<vec3::type> frc_x, frc_y, frc_z;
            alignas(64) std::vector<vec3::type> old_x, old_y, old_z;

            alignas(64) std::vector<double> mass;
            alignas(64) std::vector<env::ParticleState> state;
            alignas(64) std::vector<env::ParticleType> type;
            alignas(64) std::vector<env::ParticleID> id;
            alignas(64) std::vector<U> user_data;

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
            U * AP_RESTRICT ptr_user_data = nullptr;

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

        Storage data;
        std::vector<uint32_t> id_to_index_map;

        // explode AoS input into SoA vectors
        void build_storage(const std::vector<env::internal::ParticleRecord<U>>& particles) {
            AP_ASSERT(!is_built, "storage already built");

            size_t n = particles.size();
            data.resize(n);
            id_to_index_map.resize(n);

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
            is_built = true;
        }

        void swap_particles(size_t i, size_t j) {
            data.swap(i, j);
            std::swap(id_to_index_map[data.id[i]], id_to_index_map[data.id[j]]); // update map
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
        std::vector<batching::TopologyBatch> topology_batches;
        bool is_built = false;
    };
}
