#pragma once
#include <vector>
#include "april/containers/container.hpp"
#include "april/containers/batch.hpp"
#include "april/particle/particle.hpp"

namespace april::container {

    template<typename Config, env::IsUserData U>
    class SoAContainer : public Container<Config, U> {
    public:
        using Base = Container<Config, U>;
        using Particle = env::internal::ParticleRecord<U>;
        using Base::force_schema;
        using Base::Base;

        SoAContainer(const Config & config, const internal::ContainerCreateInfo & info)
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

        // explode AoS input into SoA vectors
        void build_storage(const std::vector<Particle>& particles) {
            AP_ASSERT(!is_built, "storage already built");

            size_t n = particles.size();
            pos_x.resize(n); pos_y.resize(n); pos_z.resize(n);
            vel_x.resize(n); vel_y.resize(n); vel_z.resize(n);
            frc_x.resize(n); frc_y.resize(n); frc_z.resize(n);
            old_x.resize(n); old_y.resize(n); old_z.resize(n);
            mass.resize(n); state.resize(n); type.resize(n); id.resize(n);
            user_data.resize(n);
            id_to_index_map.resize(n);

            for (size_t i = 0; i < n; ++i) {
                const auto& p = particles[i];

                // Vectors
                pos_x[i] = p.position.x;        pos_y[i] = p.position.y;        pos_z[i] = p.position.z;
                vel_x[i] = p.velocity.x;        vel_y[i] = p.velocity.y;        vel_z[i] = p.velocity.z;
                frc_x[i] = p.force.x;           frc_y[i] = p.force.y;           frc_z[i] = p.force.z;
                old_x[i] = p.old_position.x;    old_y[i] = p.old_position.y;    old_z[i] = p.old_position.z;

                // Scalars
                mass[i] = p.mass;
                state[i] = p.state;
                type[i] = p.type;
                id[i] = p.id;
                user_data[i] = p.user_data;

                // ID Map
                id_to_index_map[static_cast<size_t>(p.id)] = i;
            }
            is_built = true;
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
            return pos_x.size();
        }

    protected:
        void swap_particles(size_t i, size_t j) {
            if (i == j) return;

            // Swap Vectors
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

            // Update Map
            std::swap(id_to_index_map[id[i]], id_to_index_map[id[j]]);
        }

        // deducing this handles const correctness automatically via CTAD
        template<env::Field F>
        auto get_field_ptr(this auto&& self, size_t i) {

            if constexpr (F == env::Field::position) {
                return utils::Vec3Ptr { &self.pos_x[i], &self.pos_y[i], &self.pos_z[i] };
            }
            else if constexpr (F == env::Field::velocity) {
                return utils::Vec3Ptr { &self.vel_x[i], &self.vel_y[i], &self.vel_z[i] };
            }
            else if constexpr (F == env::Field::force) {
                return utils::Vec3Ptr { &self.frc_x[i], &self.frc_y[i], &self.frc_z[i] };
            }
            else if constexpr (F == env::Field::old_position) {
                return utils::Vec3Ptr { &self.old_x[i], &self.old_y[i], &self.old_z[i] };
            }
            // Scalars are simple pointers
            else if constexpr (F == env::Field::mass)      return &self.mass[i];
            else if constexpr (F == env::Field::state)     return &self.state[i];
            else if constexpr (F == env::Field::type)      return &self.type[i];
            else if constexpr (F == env::Field::id)        return &self.id[i];
            else if constexpr (F == env::Field::user_data) return &self.user_data[i];
        }

        // SoA Storage
        std::vector<double> pos_x, pos_y, pos_z;
        std::vector<double> vel_x, vel_y, vel_z;
        std::vector<double> frc_x, frc_y, frc_z;
        std::vector<double> old_x, old_y, old_z;

        std::vector<double> mass;
        std::vector<env::ParticleState> state;
        std::vector<env::ParticleType> type;
        std::vector<env::ParticleID> id;
        std::vector<U> user_data; // User data remains AoS

        std::vector<uint32_t> id_to_index_map; // map id to index

    private:
        std::vector<TopologyBatch> topology_batches;
        bool is_built = false;
    };
}
