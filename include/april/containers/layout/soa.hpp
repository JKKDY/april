#pragma once
#include <vector>
#include "april/containers/container.hpp"
#include "../../exec/executors/omp_executor.hpp"
#include "april/particle/particle.hpp"
#include "april/exec/policy.hpp"


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
        ParticleState * AP_RESTRICT ptr_state = nullptr;
        ParticleType  * AP_RESTRICT ptr_type  = nullptr;
        ParticleID    * AP_RESTRICT ptr_id    = nullptr;
        Attributes * AP_RESTRICT ptr_attributes = nullptr;

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

        void resize(const size_t n) {
            capacity = n + packed::size();
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
            std::swap(attributes[i], attributes[j]);
        }
    };



    template<typename Config, particle::IsParticleAttributes Attributes>
    class SoA : public Container<Config, Attributes> {
    public:
        using Base = Container<Config, Attributes>;
        using Base::force_schema;
        using Base::Base;
        friend Base;

        // INDEXING
        [[nodiscard]] size_t id_to_index(const ParticleID id) const {
            return id_to_index_map[static_cast<size_t>(id)];
        }
        [[nodiscard]] ParticleID min_id() const {
            return 0;
        }
        [[nodiscard]] ParticleID max_id() const {
            return static_cast<ParticleID>(id_to_index_map.size());
        }
        [[nodiscard]] bool index_is_valid(const size_t index) const {
            return index < particle_count();
        }
        [[nodiscard]] bool contains_id(const ParticleID id) const {
            return id <= max_id();
        }


        // QUERIES
        [[nodiscard]] size_t capacity() const {
            return data.capacity;
        }
        [[nodiscard]] size_t particle_count() const {
            return data.size;
        }

    protected:
        SoAStorage<Attributes> tmp;
        SoAStorage<Attributes> data;
        std::vector<size_t> bin_starts; // first particle index of each bin
        std::vector<size_t> bin_sizes; // number of particles in each bin
        std::vector<uint32_t> id_to_index_map;

        // explode AoS input into SoA vectors
        void build_storage(const std::vector<particle::ParticleRecord<Attributes>>& particles) {
            size_t n = particles.size();
            data.resize(n);
            id_to_index_map.resize(n);

            bin_starts.clear();
            bin_sizes.clear();
            bin_starts.push_back(0);
            bin_sizes.push_back(particles.size());

            for (size_t i = 0; i < particle_count(); ++i) {
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
                data.attributes[i] = p.attributes;

                // ID Map
                id_to_index_map[static_cast<size_t>(p.id)] = i;
            }

             for (size_t i = particle_count(); i < capacity(); ++i) {
                 data.pos_x[i] = 1e50;
                 data.pos_y[i] = 1e50;
                 data.pos_z[i] = 1e50;
                 data.mass[i] = 1.0;
                 data.id[i] = -1;
                 data.state[i] = ParticleState::INVALID;
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
            for (size_t i = 0; i < particle_count(); i++) {
                const auto id = data.id[i];
                if (static_cast<size_t>(id) >= id_to_index_map.size()) {
                    id_to_index_map.resize(static_cast<size_t>(id) + 1, std::numeric_limits<uint32_t>::max());
                }
                id_to_index_map[static_cast<size_t>(id)] = static_cast<uint32_t>(i);
            }
        }

        [[nodiscard]] math::Range get_physical_bin_range(const size_t type) const {
            const size_t start = bin_starts[type];
            return {start, start + bin_sizes[type]};
        }

        template<ParticleField F>
        auto get_field_ptr(this auto&& self, size_t i) {
            if constexpr (F == ParticleField::position)
                return math::Vec3Ptr { self.data.ptr_pos_x + i, self.data.ptr_pos_y + i, self.data.ptr_pos_z + i };
            else if constexpr (F == ParticleField::velocity)
                return math::Vec3Ptr { self.data.ptr_vel_x + i, self.data.ptr_vel_y + i, self.data.ptr_vel_z + i };
            else if constexpr (F == ParticleField::force)
                return math::Vec3Ptr { self.data.ptr_frc_x + i, self.data.ptr_frc_y + i, self.data.ptr_frc_z + i };
            else if constexpr (F == ParticleField::old_position)
                return math::Vec3Ptr { self.data.ptr_old_x + i, self.data.ptr_old_y + i, self.data.ptr_old_z + i };

            else if constexpr (F == ParticleField::mass)      return self.data.ptr_mass + i;
            else if constexpr (F == ParticleField::state)     return self.data.ptr_state + i;
            else if constexpr (F == ParticleField::type)      return self.data.ptr_type + i;
            else if constexpr (F == ParticleField::id)        return self.data.ptr_id + i;
            else if constexpr (F == ParticleField::attributes) return self.data.ptr_attributes + i;
        }

        template<exec::ExecutionMode E, bool is_const, exec::IsKernel Kernel>
        void iterate_chunk(this auto&& self, Kernel && kernel, const size_t c_start, const size_t c_end) {
            using K = std::remove_cvref_t<Kernel>;

             if constexpr (E == exec::ExecutionMode::Scalar) {
                for (size_t i = c_start; i < c_end; i++) {
                    if constexpr (is_const) {
                        kernel(i, self.template view<K::Read>(i));
                    } else {
                        kernel(i, self.template at<K::Read, K::Write>(i));
                    }
                }
            }
            else if constexpr (E == exec::ExecutionMode::Vector) {
                for (size_t i = c_start; i < c_end; i += packed::size()) {
                    AP_ASSERT(i % packed::size() == 0, "In vectorized execution start must be aligned to the packed type");
                    if constexpr (is_const) {
                        kernel(i, self.template view_packed<K::Read>(i));
                    } else {
                        kernel(i, self.template at_packed<K::Read, K::Write>(i));
                    }
                }
            }
            else if constexpr (E == exec::ExecutionMode::Hybrid) {
                // head
                constexpr size_t vector_size = packed::size();
                const size_t remainder = c_start % vector_size;
                const size_t head_end = (remainder == 0) ? c_start : std::min(c_end, c_start + (vector_size - remainder));

                for (size_t i = c_start; i < head_end; ++i) {
                    if constexpr (is_const) {
                        kernel(i, self.template view<K::Read>(i));
                    } else {
                        kernel(i, self.template at<K::Read, K::Write>(i));
                    }
                }

                // body
                const size_t body_start = head_end;
                const size_t body_end = body_start + ((c_end - body_start) / vector_size) * vector_size;

                for (size_t i = body_start; i < body_end; i += vector_size) {
                    AP_ASSERT(i % packed::size() == 0, "In vectorized execution, index must be aligned");
                    if constexpr (is_const) {
                        kernel(i, self.template view_packed<K::Read>(i));
                    } else {
                        kernel(i, self.template at_packed<K::Read, K::Write>(i));
                    }
                }

                // tail
                for (size_t i = body_end; i < c_end; ++i) {
                    if constexpr (is_const) {
                        kernel(i, self.template view<K::Read>(i));
                    } else {
                        kernel(i, self.template at<K::Read, K::Write>(i));
                    }
                }
            }
        }

        template<ParallelPolicy P, exec::ExecutionMode E, bool is_const, exec::IsKernel Kernel>
        void iterate_range(this auto&& self, Kernel && kernel, const size_t start, const size_t end) {
            auto process_chunk = [&](const size_t c_start, const size_t c_end) {
                self.template iterate_chunk<E, is_const>(kernel, c_start, c_end);
            };

            // Dispatch based on ParallelPolicy
            if constexpr (P == ParallelPolicy::Serial) {
                process_chunk(start, end);
            }
            else if constexpr (P == ParallelPolicy::Threaded) {
                const size_t total_elements = end - start;
                if (total_elements == 0) return;

                constexpr size_t v_size = packed::size();
                const size_t total_packed = total_elements / v_size;

                // Determine block count (2x oversubscription)
                size_t B = self.thread_executor.num_threads() * 2;
                if (B == 0) B = 1;

                // number chunks should not be more then number packed
                if (total_packed < B) {
                    B = std::max<size_t>(1, total_packed);
                }

                if (B <= 1) {
                    process_chunk(start, end);
                    return;
                }

                const size_t vectors_per_block = total_packed / B;
                const size_t remainder_vectors = total_packed % B;

                std::vector<math::Range> blocks(B);
                size_t current = start;

                for (size_t i = 0; i < B; ++i) {
                    const size_t v_count = vectors_per_block + (i < remainder_vectors ? 1 : 0);
                    size_t size = v_count * v_size;

                    if (i == B - 1) {
                        size += total_elements % v_size; // attach scalar tail
                    }

                    blocks[i] = {current, current + size};
                    current += size;
                }

                self.thread_executor.execute(B, [&](const size_t i) {
                    if (blocks[i].start < blocks[i].stop) {
                        process_chunk(blocks[i].start, blocks[i].stop);
                    }
                });
            }
        }
    };
}


