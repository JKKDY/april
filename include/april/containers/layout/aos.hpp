#pragma once

#include "april/particle/particle.hpp"
#include "april/containers/container.hpp"
#include "april/math/range.hpp"
#include "april/exec/policy.hpp"
#include "april/exec/parallel_utils.hpp"

namespace april::container::layout {
    template <typename ContainerConfig>
    class AoS : public Container<ContainerConfig> {
    public:
        using Base = Container<ContainerConfig>;
        using Base::interaction_map;
        using Base::thread_executor;
        friend Base;

        using Particle = Base::ParticleRecord;

        explicit AoS(const ContainerConfig& config) :
            Base(config) {
        }

        void bind_executor(Base::ThreadExecutor* raw_executor_ptr) {
            thread_executor.bind(raw_executor_ptr);
            this->pair_schedule_config = exec::BlockConfig(thread_executor.num_threads(), 2);
            this->linear_schedule_config = exec::BlockConfig(thread_executor.num_threads(), 8);
        }


        // INDEXING
        [[nodiscard]] size_t id_to_index(const ParticleID id) const {
            return id_to_index_map[static_cast<size_t>(id)];
        }

        [[nodiscard]] ParticleID min_id() const {
            return 0;
        }

        [[nodiscard]] ParticleID max_id() const {
            return static_cast<ParticleID>(particles.size());
        }

        [[nodiscard]] bool contains_id(const ParticleID id) const {
            return id < max_id();
        }

        [[nodiscard]] bool index_is_valid(const size_t index) const {
            return index < particle_count();
        }


        // QUERIES
        [[nodiscard]] size_t capacity() const {
            return particle_count();
        }

        [[nodiscard]] size_t particle_count() const {
            return particles.size();
        }


        // DISABLE PACKED ACCESS
        template <ParticleField R, ParticleField W>
        [[nodiscard]] auto at_packed(this auto&&, size_t) {
            static_assert(false, "AoS does not support packed access");
        }

        template <ParticleField R>
        [[nodiscard]] auto view_packed(this const auto&, size_t) {
            static_assert(false, "AoS does not support packed access");
        }

    protected:
        std::vector<Particle> tmp = {};
        std::vector<Particle> particles = {};
        std::vector<size_t> bin_starts; // first particle index of each bin
        std::vector<size_t> bin_sizes; // number of particles in each bin
        std::vector<uint32_t> id_to_index_map; // map id to index

        exec::BlockConfig pair_schedule_config;
        exec::BlockConfig linear_schedule_config;

        void build_storage(const std::vector<Particle>& particles_in) {
            particles = std::vector(particles_in);
            bin_starts.clear();
            bin_sizes.clear();
            bin_starts.push_back(0);
            bin_sizes.push_back(particles.size());
            id_to_index_map.resize(particles.size());
            for (size_t i = 0; i < particles.size(); i++) {
                const auto id = static_cast<size_t>(particles[i].id);
                id_to_index_map[id] = i;
            }

            tmp.resize(particles.size());
        }


        std::vector<size_t> bin_counts;
        std::vector<size_t> cached_bins;
        std::vector<size_t> bin_counts_tls_buffers;
        std::vector<size_t> bin_particles;

        template <typename HashFunc>
        requires std::invocable<HashFunc, size_t> &&
            std::unsigned_integral<std::invoke_result_t<HashFunc, size_t>>
        void reorder_storage(const size_t n_bins, HashFunc&& calc_bin) {
            const size_t n_particles = this->particle_count();
            const unsigned n_threads = this->thread_executor.num_threads();

            if (n_particles == 0 || n_bins == 0) return;

            // Fast serial resizes
            bin_starts.resize(n_bins);
            bin_sizes.resize(n_bins);
            bin_counts.resize(n_bins);
            cached_bins.resize(n_particles);
            bin_particles.resize(n_particles);
            bin_counts_tls_buffers.assign(n_bins * n_threads, 0);

            auto particle_blocks = exec::make_linear_schedule(math::Range{0, n_particles}, this->linear_schedule_config);
            auto bin_blocks = exec::make_linear_schedule(math::Range{0, n_bins}, this->linear_schedule_config);

            // each thread counts the number of particles per bin in its assigned particle blocks
            this->thread_executor.execute(particle_blocks.size(), [&](const size_t t_idx) {
                APRIL_ASSERT(exec::thread_index() < this->thread_executor.num_threads(),
                        "[APRIL] exec::thread_index() must be in range [0, executor.num_threads()]. "
                        "Verify that the executor sets ScopedThreadContext correctly");
                const auto& block = particle_blocks[t_idx];
                auto* bin_counts_tls = &bin_counts_tls_buffers[exec::thread_index() * n_bins];

                for (size_t i = block.start; i < block.stop; ++i) {
                    const size_t bin = calc_bin(i);
                    cached_bins[i] = bin;
                    ++bin_counts_tls[bin];
                }
            });

            // reduce to a single buffer
            this->thread_executor.execute(bin_blocks.size(), [&](const size_t t_idx) {
                const auto& bins = bin_blocks[t_idx];

                for (const size_t bin : bins) bin_sizes[bin] = 0;

                for (unsigned i = 0; i < n_threads; i++) {
                    for (const size_t bin : bins) {
                        bin_sizes[bin] += bin_counts_tls_buffers[n_bins * i + bin];
                    }
                }
            });

            // compute offsets (prefix sum)
            size_t offset = 0;
            for (size_t i = 0; i < n_bins; i++) {
                bin_starts[i] = offset;
                bin_counts[i] = offset;
                offset += bin_sizes[i];
            }

            // Build the mapping: destination_idx -> source_idx
            for (size_t i = 0; i < n_particles; i++) {
                const size_t bin = cached_bins[i];
                const size_t dest_idx = bin_counts[bin]++;
                bin_particles[dest_idx] = i;
            }

            // gather copy into ping pong buffer
            this->thread_executor.execute(particle_blocks.size(), [&](const size_t t_idx) {
                 const auto& block = particle_blocks[t_idx];

                 for (size_t dest_idx = block.start; dest_idx < block.stop; ++dest_idx) {
                     const size_t src_idx = bin_particles[dest_idx];
                     tmp[dest_idx] = particles[src_idx];

                     // Update ID map
                     const auto id = static_cast<size_t>(tmp[dest_idx].id);
                     id_to_index_map[id] = static_cast<uint32_t>(dest_idx);
                 }
             });

            // ping pong swap
            std::swap(particles, tmp);
        }

        [[nodiscard]] math::Range get_physical_bin_range(const size_t type) const {
            const size_t start = bin_starts[type];
            return {start, start + bin_sizes[type]};
        }

        // Deducing 'this' automatically propagates constness to the return type
        template <ParticleField F>
        auto get_field_ptr(this auto&& self, size_t i) {
            if constexpr (F == ParticleField::force) return &self.particles[i].force;
            else if constexpr (F == ParticleField::position) return &self.particles[i].position;
            else if constexpr (F == ParticleField::velocity) return &self.particles[i].velocity;
            else if constexpr (F == ParticleField::old_position) return &self.particles[i].old_position;
            else if constexpr (F == ParticleField::mass) return &self.particles[i].mass;
            else if constexpr (F == ParticleField::state) return &self.particles[i].state;
            else if constexpr (F == ParticleField::type) return &self.particles[i].type;
            else if constexpr (F == ParticleField::id) return &self.particles[i].id;
            else if constexpr (F == ParticleField::attributes) return &self.particles[i].attributes;
        }


        template <ParallelPolicy P, exec::ExecutionMode V, bool is_const, exec::IsKernel Kernel>
        void iterate_range(this auto&& self, Kernel&& kernel, const size_t start, const size_t end) {
            static_assert(V != exec::ExecutionMode::Vector,
                          "AoS cannot be vectorized. Change the vector policy to scalar or auto.");

            auto run_kernel = [&](size_t i) APRIL_FORCE_INLINE {
                using K = std::remove_cvref_t<Kernel>;
                if constexpr (is_const) {
                    kernel(i, self.template view<K::Read>(i));
                } else {
                    kernel(i, self.template at<K::Read, K::Write>(i));
                }
            };

            math::Range full_range = {start, end};
            if constexpr (P == ParallelPolicy::Threaded) {
                auto schedule = exec::make_linear_schedule(full_range, self.linear_schedule_config);

                self.thread_executor.execute(schedule.size(), [&](size_t i) {
                    for (size_t j : schedule[i]) run_kernel(j);
                });
            } else if constexpr (P == ParallelPolicy::Serial){
                for (size_t j : full_range) run_kernel(j);
            } else {
                static_assert(false, "invalid Parallel Policy");
            }
        }
    };
}
