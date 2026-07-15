#pragma once

#include <vector>
#include <cmath>
#include <algorithm>
#include <utility>

#include "april/base/types.hpp"
#include "april/containers/batching/topology_batch.hpp"
#include "april/particle/properties.hpp"
#include "april/core/domain.hpp"
#include "april/exec/kernel.hpp"
#include "april/exec/threading/scheduling.hpp"

namespace april::container::internal {

	template <class Base>
	class DirectSumCore : public Base {
	public:
		using Base::Base;
		using Base::build_storage;
		using typename Base::ParticleRecord;
		using Base::vector_policy;
		using Base::parallel_policy;

		void build(this auto&& self, const std::vector<ParticleRecord>& particles) {
			self.build_storage(particles);
			self.build_type_batches();
			self.build_topology_batches();
		}

		template<ParallelPolicy P, typename Func>
		void for_each_topology_batch(this auto&& self, Func && f) {
			for (const auto& phase : self.topology_phases) {
				self.thread_executor.template execute<P>(phase.size(), [&](size_t i) {
					f(phase[i]);
				});
			}
		}

		template<ParallelPolicy P, typename F>
		void for_each_interaction_batch(this auto&& self, F && f) {
			// periodicity flags to jump table index
			const int mode = (self.flags.periodic_x ? 4 : 0) |
				(self.flags.periodic_y ? 2 : 0) |
				(self.flags.periodic_z ? 1 : 0);

			// trampoline lambda
			auto run_with_mode = [&] <bool PX, bool PY, bool PZ>() {
				auto bcp = [L=self.domain.extent]<typename T>(const math::Vec3<T>& dr) {
					if constexpr (!std::is_floating_point_v<T>) {
						auto res = dr;
						if constexpr (PX) res.x -= L.x * round(dr.x / L.x);
						if constexpr (PY) res.y -= L.y * round(dr.y / L.y);
						if constexpr (PZ) res.z -= L.z * round(dr.z / L.z);
						return res;
					} else {
						return minimum_image<PX, PY, PZ>(dr, L);
					}
				};

				// process symmetric groups
				for (const auto & sym_group : self.sym_groups) {
					self.thread_executor.template execute<P>(sym_group.diagonals.size(), [&](size_t i) {
						f(sym_group.diagonals[i], bcp);
					});
					for (const auto & off_diag : sym_group.off_diagonals) {
						self.thread_executor.template execute<P>(off_diag.size(), [&](size_t i) {
							f(off_diag[i], bcp);
						});
					}
				}

				// process asymmetric groups
				for (const auto & asym_group : self.asym_groups) {
					for (const auto & phase : asym_group.phases) {
						self.thread_executor.template execute<P>(phase.size(), [&](size_t i) {
							f(phase[i], bcp);
						});
					}
				}
			};

			// jump table
			switch(mode) {
				case 0: run_with_mode.template operator()<false, false, false>(); break;
				case 1: run_with_mode.template operator()<false, false, true >(); break;
				case 2: run_with_mode.template operator()<false, true , false>(); break;
				case 3: run_with_mode.template operator()<false, true , true >(); break;
				case 4: run_with_mode.template operator()<true , false, false>(); break;
				case 5: run_with_mode.template operator()<true , false, true >(); break;
				case 6: run_with_mode.template operator()<true , true , false>(); break;
				case 7: run_with_mode.template operator()<true , true , true >(); break;
				default: std::unreachable();
			}
		}

		[[nodiscard]] std::vector<size_t> collect_indices_in_region(this const auto& self, const core::Box & region) {
			std::vector<size_t> ret;
			const double domain_vol = self.domain.volume();
			const auto intersection = self.domain.intersection(region);

			// partition the entire particle range into independent blocks/tasks
			auto blocks = exec::make_linear_schedule(math::Range{0, self.capacity()}, self.linear_schedule_config);
			const size_t num_tasks = blocks.size();

			// allocate buffers for each task
			std::vector<std::vector<size_t>> local_results(self.thread_executor.num_threads());

			// preallocate storage with heuristic (assumes uniform distribution)
			if (domain_vol > 1e-9 && intersection.has_value()) {
				const double ratio = intersection->volume() / domain_vol;
				const auto est_total = static_cast<size_t>(self.particle_count() * ratio);
				const size_t est_per_task = (est_total / num_tasks) + 1;

				for (auto& loc : local_results) {
					loc.reserve(est_per_task);
				}
			}

			// process all tasks in parallel
			self.thread_executor.template execute<parallel_policy>(num_tasks, [&](const size_t t_idx) {
				const auto& block = blocks[t_idx];
				auto& local_ret = local_results[exec::thread_index()];

				// kernel checks if a particle is alive and inside the region
				auto kernel = april::universal_kernel<ParticleField::position | ParticleField::state> (
					[&]<bool is_packed>(const size_t i, const auto & particle) {
						if constexpr (is_packed) {
							// vectorized state and contains check
							auto contains_particle = region.contains(particle.position);
							auto is_alive = (particle.state == +ParticleState::ALIVE);

							// combine and export as integer bit mask
							auto valid_mask = contains_particle & is_alive;
							uint64_t bitmask = valid_mask.to_bitmask();

							// extract exact lanes without branching
							while (bitmask != 0) {
								 const uint32_t lane = std::countr_zero(bitmask);
								 local_ret.push_back(i + lane);
								 bitmask &= (bitmask - 1); // clears the lowest set bit
							}
						} else {
							if (region.contains(particle.position) && particle.state == ParticleState::ALIVE) {
								local_ret.push_back(i);
							}
						}
					}
				);

				// runt the kernel
				self.for_each_particle(block.start, block.stop, kernel);
			});

			// allocate storage for indices buffer
			size_t total_found = 0;
			for (const auto& loc : local_results) {
				total_found += loc.size();
			}

			std::vector<size_t> indices;
			ret.reserve(total_found);

			// merge all local buffers into indices buffer
			for (const auto& loc : local_results) {
				ret.insert(ret.end(), loc.begin(), loc.end());
			}

			return ret;
		}

		void rebuild_structure() { /*NoOp: nothing to rebuild for direct sum*/ }

	private:
		std::vector<std::vector<batching::TopologyBatch<DirectSumCore>>> topology_phases;

		void build_topology_batches() {

			// collect all interaction topologies into a single vector
			std::vector<utility::graph::EdgeList<ParticleID>> global_topologies;
			global_topologies.reserve(this->interaction_map.interactions.size());

			for (const auto& prop : this->interaction_map.interactions) {
				if (!prop.used_by_ids.empty() && prop.is_active) {
					global_topologies.push_back(prop.used_by_ids);
				}
			}

			// create schedule
			constexpr size_t max_partition_size = 1024;
			const size_t min_batches_threshold = this->thread_executor.num_threads();
			auto scheduled_phases = batching::build_concurrent_phases<ParticleID>(
				global_topologies,
				max_partition_size,
				min_batches_threshold
			);

			// build phases into batches
			using ContainerType = std::remove_cvref_t<decltype(*this)>;

			for (auto& phase : scheduled_phases) {
				std::vector<batching::TopologyBatch<ContainerType>> current_phase_batches;
				current_phase_batches.reserve(phase.size());

				for (auto& batch_pairs : phase) {
					if (batch_pairs.empty()) continue;

					batching::TopologyBatch<ContainerType> batch;
					batch.container_ptr = this;
					batch.representatives = batch_pairs[0];
					batch.pairs = std::move(batch_pairs);

					current_phase_batches.push_back(std::move(batch));
				}
				topology_phases.push_back(std::move(current_phase_batches));
			}
		}


		void build_type_batches(this auto&& self) {

			// reorder into type-buckets
			const size_t n_bins = self.interaction_map.types.size();
			self.reorder_storage(n_bins,[&](const size_t i) APRIL_FORCE_INLINE {
				return static_cast<size_t>(*self.template get_field_ptr<ParticleField::type>(i));
			});

			// build bucket sorted storage and create batches
			self.generate_batches();
		}

		void generate_batches(this auto && self) {
			// if (self.bin_starts.empty()) return;
			const auto n_types = static_cast<ParticleType>(self.interaction_map.types.size());

            // symmetric batches
            for (ParticleType type = 0; type < n_types; type++) {
                self.generate_symmetric_group(type);
            }

            // asymmetric batches
            for (ParticleType t1 = 0; t1 < n_types; t1++) {
                for (ParticleType t2 = t1 + 1; t2 < n_types; t2++) {
                    self.generate_asymmetric_group(t1, t2);
                }
            }
        }

        void generate_symmetric_group(this auto && self, ParticleType type) {
            using Derived = std::remove_cvref_t<decltype(self)>;
            using SymGroup = Derived::SymTaskGroup;

            auto range = self.get_physical_bin_range(type);
            if (range.size() <= 1) return;

            auto schedule = exec::make_symmetric_schedule(range, self.pair_schedule_config);
            SymGroup group;

            // diagonals
            group.diagonals.reserve(schedule.diagonals.size());
            for (const auto& r : schedule.diagonals) {
                group.diagonals.push_back(self.create_symmetric_batch(type, r));
            }

            // off diagonals
            group.off_diagonals.resize(schedule.off_diagonals.size());
            for (size_t phase = 0; phase < schedule.off_diagonals.size(); ++phase) {
                for (const auto& [r1, r2] : schedule.off_diagonals[phase]) {
                    group.off_diagonals[phase].push_back(self.create_asymmetric_batch(type, r1, type, r2));
                }
            }

            self.sym_groups.push_back(std::move(group));
        }

        void generate_asymmetric_group(this auto && self, ParticleType type1, ParticleType type2) {
            using Derived = std::remove_cvref_t<decltype(self)>;
            using AsymGroup = Derived::AsymTaskGroup;

            auto range1 = self.get_physical_bin_range(type1);
            auto range2 = self.get_physical_bin_range(type2);
            if (range1.empty() || range2.empty()) return;

            auto schedule = exec::make_bipartite_schedule(range1, range2, self.pair_schedule_config);
            AsymGroup group;
            group.phases.resize(schedule.phases.size());

            for (size_t p = 0; p < schedule.phases.size(); ++p) {
                for (const auto& [r1, r2] : schedule.phases[p]) {
                    group.phases[p].push_back(self.create_asymmetric_batch(type1, r1, type2, r2));
                }
            }
            self.asym_groups.push_back(std::move(group));
        }


		template<bool PX, bool PY, bool PZ>
		static vec3 minimum_image(vec3 dr, const vec3& L) noexcept {
			if constexpr (PX) dr.x -= L.x * std::round(dr.x / L.x);
			if constexpr (PY) dr.y -= L.y * std::round(dr.y / L.y);
			if constexpr (PZ) dr.z -= L.z * std::round(dr.z / L.z);
			return dr;
		}
	};
}


