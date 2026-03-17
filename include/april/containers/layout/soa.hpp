#pragma once
#include <vector>
#include "april/containers/container.hpp"
#include "april/exec/parallel_utils.hpp"
#include "april/particle/particle.hpp"
#include "april/exec/policy.hpp"

#include "april/containers/layout/internal/soa_storage.hpp"


namespace april::container::layout {

    template<typename Config, particle::IsParticleAttributes Attributes>
    class SoA : public Container<Config, Attributes> {
    public:
        using Base = Container<Config, Attributes>;
        using Base::force_schema;
        friend Base;

        SoA(const Config & config, const internal::ContainerCreateInfo & info, const exec::Executor & executor):
            Base(config, info, executor)
        {
            this->pair_schedule_config = exec::BlockConfig(executor.num_threads(), 2);
            this->linear_schedule_config = exec::BlockConfig(executor.num_threads(), 8);
            for (size_t k = 0; k < packed::size(); ++k) idx_arr[k] = static_cast<double>(k);
        }

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

            // insert particles into storage
            for (size_t i = 0; i < particle_count(); ++i) {
                const auto& p = particles[i];
                data.insert_particle(i, p);
                id_to_index_map[static_cast<size_t>(p.id)] = i;
            }

             // pad with garbage data
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


        void reorder_storage(const std::vector<std::vector<size_t>>& new_bins) {
            bin_starts.clear();
            bin_sizes.clear();

            // calculate offset of each bin sequentially
            std::vector<size_t> offsets(new_bins.size());
            size_t current_offset = 0;

            for (size_t i = 0; i < new_bins.size(); ++i) {
                bin_starts.push_back(current_offset);
                bin_sizes.push_back(new_bins[i].size());
                offsets[i] = current_offset;
                current_offset += new_bins[i].size();
            }

            // schedule tasks over the bins rather than the particles inside them
            auto blocks = exec::make_linear_schedule(math::Range{0, new_bins.size()}, this->linear_schedule_config);

            // execute thread pool exactly once
            this->thread_executor.execute(blocks.size(), [&](const size_t b_idx) {
                const auto& block = blocks[b_idx];

                for (size_t bin_idx = block.start; bin_idx < block.stop; ++bin_idx) {
                    const auto& bin = new_bins[bin_idx];
                    if (bin.empty()) continue;

                    const size_t start_offset = offsets[bin_idx];

                    for (size_t i = 0; i < bin.size(); ++i) {
                        const size_t old_idx = bin[i];
                        const size_t new_idx = start_offset + i;

                        // copy particle
                        tmp.copy_from(new_idx, data, old_idx);

                        // update id map
                        const auto id = static_cast<size_t>(data.id[old_idx]);
                        id_to_index_map[id] = static_cast<uint32_t>(new_idx);
                    }
                }
            });

            // swap old and new storage
            std::swap(data, tmp);

            // update pointer caches to reflect the swap
            data.update_pointer_cache();
            tmp.update_pointer_cache();
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


        template<ParallelPolicy P, exec::ExecutionMode E, bool is_const, exec::IsKernel Kernel>
        void iterate_range(this auto&& self, Kernel && kernel, const size_t start, const size_t end) {
            using K = std::remove_cvref_t<Kernel>;
            math::Range range = {start, end};

            auto get_scalar = [&](size_t i) AP_FORCE_INLINE {
                if constexpr (is_const) return self.template view<K::Read>(i);
                else return self.template at<K::Read, K::Write>(i);
            };

            auto get_vector_ref = [&](size_t i) AP_FORCE_INLINE {
                if constexpr (is_const) return self.template view_packed<K::Read>(i);
                else return self.template at_packed<K::Read, K::Write>(i);
            };

            // scalar/vector routing
            auto process_chunk = [&](const math::Range & chunk) {
                // process a chunk of work in scalar mode
                if constexpr (E == exec::ExecutionMode::Scalar) {
                    for (size_t i : chunk) {
                        kernel(i, get_scalar(i));
                    }
                }

                // process a chunk of work in vector mode
                else if constexpr (E == exec::ExecutionMode::Vector || E == exec::ExecutionMode::Hybrid) {
                    const size_t body = chunk.size() / packed::size();
                    const size_t tail = chunk.size() % packed::size();

                    //body (loads are unaligned so no need for a head)
                    for (size_t i = 0; i < body; i++) {
                        size_t idx = chunk.start + i * packed::size();
                        kernel(idx, get_vector_ref(idx));
                    }

                    // tail (with mask)
                    if (tail > 0) {
                        const auto lane_indices = packed::load_aligned(self.idx_arr);
                        const auto mask = lane_indices < static_cast<packed::value_type>(tail);

                        const size_t tail_idx = chunk.start + body * packed::size();
                        auto ref = get_vector_ref(tail_idx);
                        kernel(tail_idx, ref.mask_with(mask));
                    }
                } else {
                    static_assert(false, "Execution mode must be a valid value (Scalar, Vector, Hybrid)");
                }
            };

            // serial/multi-threaded routing
            if constexpr (P == ParallelPolicy::Serial) {
                process_chunk(range);
            }
            else if constexpr (P == ParallelPolicy::Threaded) {
                const auto blocks = exec::make_linear_schedule(range, self.linear_schedule_config);

                self.thread_executor.execute(blocks.size(), [&](const size_t i) {
                    process_chunk(blocks[i]);
                });
            }
        }

    private:
        alignas(64) packed::value_type idx_arr[packed::size()]{};
    };
}


