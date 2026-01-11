#pragma once

#include <ranges>

#include "april/containers/soa.hpp"
#include "april/containers/aos.hpp"
#include "april/containers/batching.hpp"

namespace april::container {
	namespace internal {
		template <class ContainerBase>
		class DirectSumBase final : public ContainerBase {

			struct AsymmetricBatch : SerialBatch<BatchType::Asymmetric>{
				std::ranges::iota_view<size_t, size_t> indices1{};
				std::ranges::iota_view<size_t, size_t> indices2{};
			};

			struct SymmetricBatch : SerialBatch<BatchType::Symmetric>{
				std::ranges::iota_view<size_t, size_t> indices {};
			};

		public:
			using typename ContainerBase::ParticleRecord;
			using ContainerBase::ContainerBase;

			void build(this auto&& self, const std::vector<ParticleRecord> & particlesIn) {
				self.build_storage(particlesIn);
				self.sort_storage_by_type();
				self.build_batches();
			}

			template<typename F>
				void for_each_interaction_batch(this auto&& self, F && f) {
				// periodicity flags to jump table index
				const int mode = (self.flags.periodic_x ? 4 : 0) |
							  (self.flags.periodic_y ? 2 : 0) |
							  (self.flags.periodic_z ? 1 : 0);

				// trampoline lambda
				auto run_with_mode = [&] <bool PX, bool PY, bool PZ>() {
					auto bcp = [L=self.domain.extent](const vec3& dr) {
						return minimum_image<PX, PY, PZ>(dr, L);
					};

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

			[[nodiscard]] std::vector<size_t> collect_indices_in_region(this const auto& self, const env::Box & region) {
				std::vector<size_t> ret;

				const double domain_vol = self.domain.volume();
				const auto intersection = self.domain.intersection(region);

				if (domain_vol > 1e-9 && intersection.has_value()) { // avoid division by zero
					const double ratio = intersection->volume() / domain_vol;

					// apply 1.1x safety factor because distributions are rarely perfectly uniform
					// clamp to particles.size() to avoid over-reservation
					const auto est = static_cast<size_t>(self.particle_count() * ratio * 1.1);
					ret.reserve(std::min(est, self.particle_count()));
				}

				for (size_t i = 0; i < self.particle_count(); i++) {
					const auto & particle = self.template view<env::Field::position | env::Field::state>(i);

					if (region.contains(particle.position) && particle.state != env::ParticleState::DEAD) {
						ret.push_back(i);
					}
				}

				return ret;
			}

			void rebuild_structure() {}

		private:
			std::vector<SymmetricBatch> symmetric_batches;
			std::vector<AsymmetricBatch> asymmetric_batches;

			void sort_storage_by_type(this auto&& self) {
				const size_t n = self.particle_count();
				if (n == 0) return;

				std::vector<size_t> p(n);
				std::ranges::iota(p, size_t{0});

				std::ranges::sort(p, [&](size_t i, size_t j) {
					auto t1 = *self.template get_field_ptr<env::Field::type>(i);
					auto t2 = *self.template get_field_ptr<env::Field::type>(j);
					return t1 < t2;
				});

				std::vector done(n, false);
				for (size_t i = 0; i < n; ++i) {
					if (done[i]) continue;
					size_t current = i;
					while (i != p[current]) {
						size_t next = p[current];
						self.swap_particles(current, next);
						p[current] = current;
						current = next;
					}
					done[current] = true;
				}

				// Rebuild ID map after sort
				for (size_t i = 0; i < self.particle_count(); i++) {
					auto id_ptr = self.template get_field_ptr<env::Field::id>(i);
					self.id_to_index_map[static_cast<size_t>(*id_ptr)] = i;
				}
			}

			void build_batches(this auto&& self) {
				std::unordered_map<env::ParticleType, std::ranges::iota_view<size_t, size_t>> type_ranges;
				if (self.particle_count() == 0) return;

				size_t start = 0;
				auto current_type = *self.template get_field_ptr<env::Field::type>(0);

				for (size_t i = 0; i < self.particle_count(); i++) {
					auto type = *self.template get_field_ptr<env::Field::type>(i);

					if (current_type != type) {
						type_ranges[current_type] = std::ranges::iota_view(start, i);;
						start = i;
						current_type = type;
					}
				}
				type_ranges[current_type] = std::ranges::iota_view(start, self.particle_count()); // last batch

				const env::ParticleType n_types = current_type + 1;

				// create a batch for each particle type
				for (env::ParticleType type = 0; type < n_types; type++) {
					SymmetricBatch batch;
					batch.types = {type, type};
					batch.indices = type_ranges[type];
					self.symmetric_batches.push_back(batch);
				}

				// create a batch for each distinct pair of particle types
				for (env::ParticleType t1 = 0; t1 < n_types; t1++) {
					for (env::ParticleType t2 = t1 + 1; t2 < n_types; t2++) {
						AsymmetricBatch batch;
						batch.types = {t1, t2};
						batch.indices1 = type_ranges[t1];
						batch.indices2 = type_ranges[t2];
						self.asymmetric_batches.push_back(batch);
					}
				}
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

	struct DirectSumAoS {
		template <class U>
		using impl = internal::DirectSumBase<AoSContainer<DirectSumAoS, U>>;
	};

	struct DirectSumSoA {
		template <class U>
		using impl = internal::DirectSumBase<SoAContainer<DirectSumSoA, U>>;
	};
}
