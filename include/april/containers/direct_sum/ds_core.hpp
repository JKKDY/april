#pragma once

#include <vector>
#include <cmath>
#include <algorithm>
#include <utility>

#include "april/base/types.hpp"
#include "april/particle/particle_types.hpp"
#include "april/core/domain.hpp"
#include "april/exec/particle_kernel.hpp"
#include "april/containers/batching/common.hpp"
#include "april/exec/parallel_utils.hpp"

namespace april::container::internal {

	template <class ContainerBase>
	class DirectSumCore : public ContainerBase {
	public:
		using ContainerBase::ContainerBase;
		using ContainerBase::build_storage;
		using typename ContainerBase::ParticleRecord;

		void build(this auto&& self, const std::vector<ParticleRecord>& particles) {
			// precompute topology batches (id based batches)
			for (size_t i = 0; i < self.force_schema.interactions.size(); ++i) {
				const auto& prop = self.force_schema.interactions[i];

				if (!prop.used_by_ids.empty() && prop.is_active) {
					batching::TopologyBatch batch;
					batch.id1 = prop.used_by_ids[0].first;
					batch.id2 = prop.used_by_ids[0].second;
					batch.pairs = prop.used_by_ids;

					self.topology_batches.push_back(std::move(batch));
				}
			}

			self.build_storage(particles);
			self.build_batches();
		}

		template<typename Func>
		void for_each_topology_batch(Func && func) {
			for (const auto & batch : topology_batches) {
				func(batch);
			}
		}

		template<typename F>
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

				// subclass is responsible for populating groups via generate_batches

				// process symmetric groups
				for (const auto & sym_group : self.sym_groups) {
					self.thread_executor.execute(sym_group.diagonals.size(), [&](size_t i) {
						f(sym_group.diagonals[i], bcp);
					});
					for (const auto & off_diag : sym_group.off_diagonals) {
						self.thread_executor.execute(off_diag.size(), [&](size_t i) {
							f(off_diag[i], bcp);
						});
					}
				}

				// process asymmetric groups
				for (const auto & asym_group : self.asym_groups) {
					for (const auto & phase : asym_group.phases) {
						self.thread_executor.execute(phase.size(), [&](size_t i) {
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
			std::vector<std::vector<size_t>> local_results(num_tasks);

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
			self.thread_executor.execute(num_tasks, [&](const size_t t_idx) {
				const auto& block = blocks[t_idx];
				auto& local_ret = local_results[t_idx];

				self.for_each_particle(block.start, block.stop,
					april::universal_kernel<ParticleField::position | ParticleField::state>(
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
						})
					);
				}
			);

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

		void rebuild_structure() {}

	private:
		std::vector<batching::TopologyBatch> topology_batches;

		void build_batches(this auto&& self) {
			// calculate the buckets for bucket sorting the particles by types
			// outer vector holds buckets, inner vectors hold physical indexes to particles belonging to that bucket
			std::vector<std::vector<size_t>> buckets;

			self.for_each_particle(april::scalar_kernel<ParticleField::type>(
				[&](const size_t i, const auto& p) {
					const auto type_idx = static_cast<size_t>(p.type);
					if (type_idx >= buckets.size()) {
						buckets.resize(type_idx + 1);
					}
					buckets[type_idx].push_back(i);
				}
			));

			// reorder into type-buckets
			static_assert(requires {self.reorder_storage(buckets); }, "void rebuild_storage(bins) is not implemented");
			self.reorder_storage(buckets);

			// build bucket sorted storage and create batches
			self.generate_batches();
		}

		 void generate_batches(this auto && self) {
			// if (self.bin_starts.empty()) return;
			const auto n_types = static_cast<ParticleType>(self.force_schema.types.size());

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


