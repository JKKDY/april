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

				// subclass is responsible for populating these vectors via generate_batches
				for (const auto & batch : self.symmetric_batches) f(batch, bcp);
				for (const auto & batch : self.asymmetric_batches) f(batch, bcp);
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

			// pre allocate buffer with heuristic:  1.1x safety factor * uniform distribution
			if (domain_vol > 1e-9 && intersection.has_value()) {
				const double ratio = intersection->volume() / domain_vol;
				const auto est = static_cast<size_t>(self.particle_count() * ratio * 1.1);
				ret.reserve(std::min(est, self.particle_count()));
			}

			// gather particles
			self.for_each_particle(
				april::scalar_kernel<ParticleField::position>(
					[&](const size_t i, const auto & particle) {
						if (region.contains(particle.position)) {
							ret.push_back(i);
						}
					}),
				ParticleState::ALIVE
			);

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
			static_assert(requires {self.generate_batches(); }, "void generate_batches(ranges) is not implemented");
			self.generate_batches();
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


