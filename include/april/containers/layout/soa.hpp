#pragma once
#include <vector>
#include "april/containers/container.hpp"
#include "april/exec/parallel_utils.hpp"
#include "april/particle/particle.hpp"
#include "april/exec/policy.hpp"

#include "april/containers/layout/internal/soa_storage.hpp"


namespace april::container::layout {

    template<typename ContainerConfig>
    class SoA : public Container<ContainerConfig> {
    public:
        using Base = Container<ContainerConfig>;
        using Base::interaction_map;
        using Base::thread_executor;
        using ParticleAttributes = Base::ParticleAttributes;
        friend Base;

        explicit SoA(const ContainerConfig & config): Base(config) {
            for (size_t k = 0; k < packed::size(); ++k) idx_arr[k] = static_cast<double>(k);
        }

        void bind_executor(Base::ThreadExecutor* raw_executor_ptr) {
            thread_executor.bind(raw_executor_ptr);
            pair_schedule_config = exec::BlockConfig(thread_executor.num_threads(), 2);
            linear_schedule_config = exec::BlockConfig(thread_executor.num_threads(), 8);
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
        SoAStorage<ParticleAttributes> tmp;
        SoAStorage<ParticleAttributes> data;
        // std::vector<size_t> bin_starts; // first particle index of each bin
        // std::vector<size_t> bin_sizes; // number of particles in each bin
        std::vector<uint32_t> id_to_index_map;

        exec::BlockConfig pair_schedule_config;
        exec::BlockConfig linear_schedule_config;

        // explode AoS input into SoA vectors
        void build_storage(const std::vector<particle::ParticleRecord<ParticleAttributes>>& particles) {
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

        SoAStorage<ParticleAttributes> movers_buffer;

        std::vector<size_t> bin_starts; // first particle index of each bin
        std::vector<size_t> bin_sizes; // number of particles in each bin


        std::vector<size_t> flat_bin_sizes;
        std::vector<size_t> flat_mover_counts;
        std::vector<size_t> flat_write_offsets;
        std::vector<size_t> global_mover_starts;
        std::vector<size_t> bin_mover_totals;


        void reorder_storage(const size_t n_bins, const std::vector<uint32_t>& target_bins) {
            // target bins maps index -> bin
            const size_t num_particles = this->particle_count();
            const size_t num_bins = n_bins;
            bin_sizes.resize(n_bins);


            // create schedules
            auto particle_blocks = exec::make_linear_schedule(math::Range{0, num_particles}, this->linear_schedule_config);
            auto bin_blocks = exec::make_linear_schedule(math::Range{0, num_bins}, this->linear_schedule_config);

            const size_t num_particle_tasks = particle_blocks.size();
            const size_t num_bin_tasks = bin_blocks.size();
            const size_t flat_size = num_bin_tasks * num_bins;


            // RESIZING
            // resize thread local accumulators if needed
            const size_t flat_particle_size = num_particle_tasks * num_bins;
            if (flat_bin_sizes.size() < flat_particle_size) {
                flat_bin_sizes.resize(flat_particle_size, 0);
            } else {
                std::ranges::fill(flat_bin_sizes, 0);
            }

            // calculate how large each bin is and where it starts
            // first accumulate sizes in local buffers
            this->thread_executor.execute(num_particle_tasks, [&](const size_t t_idx) {
                const auto& block = particle_blocks[t_idx];
                size_t* bin_sizes_buffer = &flat_bin_sizes[t_idx * num_bins];

                for (size_t i = block.start; i < block.stop; ++i) {
                    ++bin_sizes_buffer[target_bins[i]];
                }
            });

            // then sum all buffers up to get final bin sizes
            this->thread_executor.execute(num_bin_tasks, [&](const size_t t_idx) {
                const auto& bin_block = bin_blocks[t_idx];

                for (size_t i = bin_block.start; i < bin_block.stop; ++i) {
                    size_t sum = 0;
                    for (size_t t = 0; t < num_particle_tasks; ++t) {
                        sum += flat_bin_sizes[t * num_bins + i];
                    }
                    bin_sizes[i] = sum;
                }
            });

            // calculate bin starts (prefix sum) (serial)
            size_t total_size = 0;
            for (size_t i = 0; i < num_bins; i++) {
                bin_starts[i] = total_size;
                total_size += bin_sizes[i];
            }

            // find the number of movers per bin
            if (flat_mover_counts.size() < flat_size) {
                flat_mover_counts.resize(flat_size, 0);
            } else {
                std::ranges::fill(flat_mover_counts, 0);
            }

            this->thread_executor.execute(num_bin_tasks, [&](const size_t task_idx) {
                const auto & bin_block = bin_blocks[task_idx];
                size_t* local_counts = &flat_mover_counts[task_idx * num_bins];

                for (const size_t src_bin : bin_block) {
                   const math::Range bin_range = {bin_starts[src_bin], bin_starts[src_bin] + bin_sizes[src_bin]};
                    for (const auto i : bin_range) {
                        const size_t dst_bin = target_bins[i];
                        if (dst_bin != src_bin) { // particle i must be moved
                            local_counts[dst_bin]+=1;
                        }
                    }
                }
            });

            if (global_mover_starts.size() < num_bins) global_mover_starts.resize(num_bins, 0);
            if (flat_write_offsets.size() < flat_size) flat_write_offsets.resize(flat_size, 0);

            // TODO parallelize
            // 1. Parallel Reduction (Sum movers per bin across all threads)
            if (bin_mover_totals.size() < num_bins) bin_mover_totals.resize(num_bins, 0);

            this->thread_executor.execute(num_bin_tasks, [&](const size_t t_idx) {
                const auto& bin_block = bin_blocks[t_idx];
                for (size_t b = bin_block.start; b < bin_block.stop; ++b) {
                    size_t sum = 0;
                    for (size_t t = 0; t < num_bin_tasks; ++t) {
                        sum += flat_mover_counts[t * num_bins + b];
                    }
                    bin_mover_totals[b] = sum;
                }
            });

            // 2. Serial Prefix Sum (Calculate global starts per bin)
            size_t total_movers = 0;
            for (size_t b = 0; b < num_bins; ++b) {
                global_mover_starts[b] = total_movers;
                total_movers += bin_mover_totals[b];
            }

            std::cout << total_movers << std::endl;

            if (total_movers == 0) return;
            if (movers_buffer.capacity < total_movers) movers_buffer.resize(total_movers);

            // 3. Parallel Offset Distribution (Calculate write offsets for each thread)
            this->thread_executor.execute(num_bin_tasks, [&](const size_t t_idx) {
                const auto& bin_block = bin_blocks[t_idx];
                for (size_t b = bin_block.start; b < bin_block.stop; ++b) {
                    size_t current_offset = global_mover_starts[b];
                    for (size_t t = 0; t < num_bin_tasks; ++t) {
                        const size_t flat_idx = t * num_bins + b;
                        flat_write_offsets[flat_idx] = current_offset;
                        current_offset += flat_mover_counts[flat_idx];
                    }
                }
            });


            // ---------------------------------------------------------
            // PHASE 4: Extract & Sort directly to buffer
            // ---------------------------------------------------------
            this->thread_executor.execute(num_bin_tasks, [&](const size_t t_idx) {
                const auto& bin_block = bin_blocks[t_idx];
                size_t* local_offsets = &flat_write_offsets[t_idx * num_bins];

                for (size_t b = bin_block.start; b < bin_block.stop; ++b) {
                    const size_t owner_bin = b;
                    const math::Range bin_range = {bin_starts[b], bin_starts[b] + bin_sizes[b]};

                    for (const auto i : bin_range) {
                        const size_t target = target_bins[i];
                        if (target != owner_bin) {
                            const size_t dest_idx = local_offsets[target]++;
                            movers_buffer.copy_from(dest_idx, data, i);
                        }
                    }
                }
            });

            // ---------------------------------------------------------
            // PHASE 5: In-Place Scatter
            // ---------------------------------------------------------
            this->thread_executor.execute(num_bin_tasks, [&](const size_t t_idx) {
                const auto& bin_block = bin_blocks[t_idx];

                for (size_t b = bin_block.start; b < bin_block.stop; ++b) {
                    if (bin_sizes[b] == 0) continue;

                    size_t mover_read_idx = global_mover_starts[b];
                    const math::Range bin_range = {bin_starts[b], bin_starts[b] + bin_sizes[b]};

                    for (const auto i : bin_range) {
                        if (target_bins[i] != b) { // Native hole
                            data.copy_from(i, movers_buffer, mover_read_idx);

                            const auto id = static_cast<size_t>(data.id[i]);
                            id_to_index_map[id] = static_cast<uint32_t>(i);

                            ++mover_read_idx;
                        }
                    }
                }
            });



            // // find number of movers(holes) per bin
            // std::vector<std::vector<size_t>> mover_counts_tbuffers(num_bin_tasks, std::vector<size_t>(num_bins, 0));
            // std::vector<size_t> movers_prefix_sum(num_bins);
            //
            //
            // // movers_out_tbuffers[task_index][bin] = list of source indices (where the movers are that belong into this bin)
            // // the first copy iterated through this list and gathers all particles that are supposed to go into bin in
            // // the temporary buffer. It can index it via movers_prefix_sum[bin]
            // // basically group particles that need to be moved by their destination bin so we can write them into a single chunk into
            // // the buffer. This allows for the second copy phase to just linearly san that block (bin) of data and copy them into the holes
            // // into the according bin
            //
            // // example
            // // Bins 1,2,3,...                        01234 56 789A
            // // 1111|222|3333  ->  1211|231|3133  ->  12112|31|3133  ->  1X11X|XX|3X33  ->
            // //  current         new target bins  -> new bin sizes -> X need to be copied
            // // mapping to buffer: 1: [6, 8], 2: [1,4], 3:5
            // // buffer: 6,8|1,4|5
            // // mapping to data: 1: [1,4], 2: [5,6], 3: [8]
            //
            // std::vector<std::vector<std::vector<size_t>>> movers_out_tbuffers(num_bin_tasks);
            // for (auto & x: movers_out_tbuffers)  x.resize(num_bins);
            //
            // // movers_in_tbuffers[task_index][bin] = list of source indices (where the movers are that belong into this bin)
            // std::vector<std::vector<std::vector<size_t>>> movers_back_tbuffers(num_bin_tasks);
            // for (auto & x: movers_back_tbuffers) x.resize(num_bins);
            //
            //
            //
            // this->thread_executor.execute(num_bin_tasks, [&](const size_t task_idx) {
            //     const auto & bin_block = bin_blocks[task_idx];
            //     auto & mover_counts = mover_counts_tbuffers[task_idx];
            //     auto & movers_out = movers_out_tbuffers[task_idx];
            //     auto & movers_back = movers_back_tbuffers[task_idx];
            //
            //     for (const size_t src_bin : bin_block) {
            //         const math::Range bin_range = {bin_starts[src_bin], bin_starts[src_bin] + bin_sizes[src_bin]};
            //         for (const auto i : bin_range) {
            //             const size_t dst_bin = target_bins[i];
            //             if (dst_bin != src_bin) { // particle i must be moved
            //                 mover_counts[dst_bin]+=1;
            //                 movers_out[src_bin].push_back(i);
            //                 movers_back[src_bin].push_back(i);
            //             }
            //         }
            //     }
            // });
            //
            // // sum all thread local buffers together to get total number of movers per bin
            // std::vector<size_t> mover_counts (num_bins);
            // this->thread_executor.execute(num_bin_tasks, [&](const size_t t_idx) {
            //     const auto& bin_block = bin_blocks[t_idx];
            //
            //     for (size_t i = bin_block.start; i < bin_block.stop; ++i) {
            //         size_t sum = 0;
            //         for (const auto & x : mover_counts_tbuffers) sum += x[i];
            //         mover_counts[i] = sum;
            //     }
            // });
            //
            // size_t offset = 0;
            // for (size_t i = 0; i < num_bins; i++) {
            //     movers_prefix_sum[i] = offset;
            //     offset += mover_counts[i];
            // }
            //
            //
            //
            //
            // // find all hole indices per bin
            // std::vector<std::vector<size_t>> holes_per_bin(num_bin_tasks, std::vector<size_t>(num_bins));
            //
            //
            // this->thread_executor.execute(num_bin_tasks, [&](const size_t t_idx) {
            //     const auto & bins = bin_blocks[t_idx];
            //
            //     for (const auto src_bin: bins) {
            //         for (const auto i : src_bin) {
            //             const size_t dst_bin = target_bins[i];
            //             if (dst_bin != src_bin) {
            //                 mover_counts[dst_bin]+=1;
            //             }
            //         }
            //     }
            // });
            //
            //
            //
            //
            // // the to buffer copy iterates through this. Bin sorted by destination
            // // mover_copy1_indices_tbuffers[bin] =
            // std::vector<std::vector<size_t>> mover_copy1_indices_tbuffers(num_bins);
            //
            // // the to data copy iterates through this. Bin sorted by destination
            // // mover_copy2_indices_tbuffers[bin] = holes in bin
            // std::vector<std::vector<size_t>> mover_copy2_indices_tbuffers(num_bins);
            //
            //
            //
            //
            //
            // std::vector<std::vector<Mover>> movers_buffers(num_bins);
            //
            // // where to copy each mover to in buffer (bin sorting)
            // std::vector<size_t> global_mover_starts(num_bins, 0);
            // std::vector<size_t> current_bin_write_index(num_bin_tasks, 0);
            //
            // // track the holes for each bin
            // std::vector<std::vector<size_t>> bin_holes(num_bin_tasks, std::vector<size_t>(num_bins, 0));
            //
            //
            //
            // // identify which particles need to move and how many
            // this->thread_executor.execute(num_bin_tasks, [&](const size_t task_idx) {
            //     const auto & bin_block = bin_blocks[task_idx];
            //     auto & mover_counts = mover_counts_tbuffers[task_idx];
            //     auto & movers = movers_buffers[task_idx];
            //     auto & holes = holes_tbuffer[task_idx];
            //
            //     for (size_t b = bin_block.start; b < bin_block.stop; ++b) {
            //         const size_t src_bin = b; // this bin owns particle i
            //         const math::Range bin_range = {bin_starts[b], bin_starts[b] + bin_sizes[b]};
            //         for (const auto i : bin_range) {
            //             const size_t dst_bin = target_bins[i];
            //             if (dst_bin != src_bin) { // particle i must be moved
            //
            //                 // we get the index of the mover
            //                 // we get the bin it is in
            //                 // we get the destination bin of where it needs to go to
            //
            //                 holes[src_bin].push_back(i);
            //
            //                 movers.push_back(Mover{i, dst_bin});
            //                 mover_counts[dst_bin]+=1;
            //             }
            //         }
            //     }
            // });
            //
            //
            // this->thread_executor.execute(num_bin_tasks, [&](const size_t t_idx) {
            //     const auto& bin_block = bin_blocks[t_idx];
            //
            //     for (size_t i = bin_block.start; i < bin_block.stop; ++i) {
            //         size_t sum = 0;
            //         for (const auto & x : mover_counts_tbuffers) sum += x[i];
            //         global_mover_starts[i] = sum;
            //     }
            // });
            //
            //
            // size_t total_movers = 0;
            // for (size_t b = 0; b < num_bins; ++b) {
            //     global_mover_starts[b] = total_movers; // Where this bin's sorted block starts
            //     for (size_t t = 0; t < num_bin_tasks; ++t) {
            //         write_offsets[t][b] = total_movers; // Where this thread starts writing for this bin
            //         total_movers += mover_counts_tbuffers[t][b];
            //     }
            // }
            //






            // // TODO this should be made persistent to avoid reallocation
            // // Thread-local 2D buffer: local_movers[thread_idx][target_bin] = list of old indices
            // struct Hole {
            //     size_t index; // particle index
            //     size_t src_bin; // current bin index
            //     size_t dst_bin; // new bin index
            // };
            // std::vector<std::vector<size_t>> particles_to_move_buffers(num_bin_tasks);
            // std::vector<std::vector<Hole>> bin_holes(num_bin_tasks);

           //  // identify which particles need to move
           //  this->thread_executor.execute(num_bin_tasks, [&](const size_t task_idx) {
           //      const auto & bin_block = bin_blocks[task_idx];
           //      auto & movers = particles_to_move_buffers[task_idx];
           //      auto & holes = bin_holes[task_idx];
           //
           //      for (size_t b = bin_block.start; b < bin_block.stop; ++b) {
           //          const size_t src_bin = b; // this bin owns particle i
           //          const math::Range bin_range = {bin_starts[b], bin_starts[b] + bin_sizes[b]};
           //          for (const auto i : bin_range) {
           //              const size_t dst_bin = target_bins[i];
           //              if (dst_bin != src_bin) { // particle i must be moved
           //                  movers.push_back(i);
           //                  holes.push_back(Hole{i, src_bin, dst_bin});
           //              }
           //          }
           //      }
           //  });
           //
           //  // prefix sum offset for parallel copying of particles
           //  size_t offset = 0;
           //  std::vector<size_t> offsets(particles_to_move_buffers.size());
           //  for (size_t i = 0; i < offsets.size(); i++) {
           //      offsets[i] = offset;
           //      offset += particles_to_move_buffers[i].size();
           //  }
           //
           //  // copy particles in parallel to temporary buffer
           //  this->thread_executor.execute(particles_to_move_buffers.size(), [&](const size_t task_idx) {
           //      const auto& particle_idxs = particles_to_move_buffers[task_idx];
           //      const auto& write_offset = offsets[task_idx];
           //
           //      for (size_t i = 0; i < particle_idxs.size(); i++) {
           //          movers_buffer.copy_from(write_offset + i, data, particle_idxs[i]);
           //      }
           //  });
           //
           //  // copy mover particles back into correct bins
           //  const size_t num_movers = offset;
           //  auto mover_blocks = exec::make_linear_schedule(math::Range{0, num_movers}, this->linear_schedule_config);
           //
           //  this->thread_executor.execute(mover_blocks.size(), [&](const size_t task_idx) {
           //     const auto& movers = mover_blocks[task_idx];
           //
           //     for (const auto i : movers) {
           //         movers_buffer.copy_from(, data, particle_idxs[i]);
           //     }
           // });





            // std::vector<std::vector<size_t>> local_movers_to(num_tasks, std::vector<size_t>(num_bins, 0));
            //
            // this->thread_executor.execute(num_tasks, [&](const size_t t_idx) {
            //     const auto& block = particle_blocks[t_idx];
            //     auto& local_movers = local_movers_to[t_idx];
            //
            //     // O(1) jump to find the starting owner bin for this physical block
            //     auto it = std::upper_bound(bin_starts.begin(), bin_starts.end(), block.start);
            //     size_t owner_bin = std::distance(bin_starts.begin(), it) - 1;
            //
            //     for (size_t i = block.start; i < block.stop; ++i) {
            //         // Roll the owner bin forward if we cross a boundary
            //         // (The while loop safely handles empty bins)
            //         while (owner_bin + 1 < num_bins && i >= bin_starts[owner_bin + 1]) {
            //             owner_bin++;
            //         }
            //
            //         if (target_bins[i] != owner_bin) {
            //             local_movers[target_bins[i]]++;
            //         }
            //     }
            // });


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

            auto get_scalar = [&](size_t i) APRIL_FORCE_INLINE {
                if constexpr (is_const) return self.template view<K::Read>(i);
                else return self.template at<K::Read, K::Write>(i);
            };

            auto get_vector_ref = [&](size_t i) APRIL_FORCE_INLINE {
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


