#pragma once

#include <ranges>

#include "soa_container.h"
#include "april/containers/aos_container.hpp"
#include "april/containers/batch.hpp"

namespace april::container {
	namespace internal {
		template <class U> class DirectSumAoS;
		template <class U> class DirectSumSoA;
		template <class U> class Test;

	}

	struct DirectSumAoS {
		template<class U> using impl = internal::DirectSumAoS<U>;
	};

	struct DirectSumSoA {
		template<class U> using impl = internal::DirectSumSoA<U>;
	};



	namespace internal {
		template <class U>
		class Test final : public SoAContainer<container::DirectSumSoA, U> {
			using Base = SoAContainer<container::DirectSumSoA, U>;
			using typename Base::ParticleRecord;
			friend Base;
			using Base::domain;
			using Base::flags;
			using Base::id_to_index_map;

		public:
			using Base::particle_count;
			using Base::Base;

			void build(const std::vector<ParticleRecord> & particlesIn) {
				this->build_storage(particlesIn);
			}


			template<typename F>
			void for_each_interaction_batch(F &&) {

			}


			[[nodiscard]] std::vector<size_t> collect_indices_in_region(const env::Box &) const {
				return {};
			}

			void rebuild_structure() {}

		};


		template <class U>
		class DirectSumSoA final : public SoAContainer<container::DirectSumSoA, U> {
			using Base = SoAContainer<container::DirectSumSoA, U>;
			using typename Base::ParticleRecord;
			friend Base;
			using Base::domain;
			using Base::flags;
			using Base::id_to_index_map;
			using Base::swap_particles;

			struct AsymmetricBatch : SerialBatch<BatchSymmetry::Asymmetric>{
				std::ranges::iota_view<size_t, size_t> indices1{};
				std::ranges::iota_view<size_t, size_t> indices2{};
			};

			struct SymmetricBatch : SerialBatch<BatchSymmetry::Symmetric>{
				std::ranges::iota_view<size_t, size_t> indices {};
			};

		public:
			using Base::particle_count;
			using Base::Base;

			void build(const std::vector<ParticleRecord> & particlesIn) {
				this->build_storage(particlesIn);

				sort_storage_by_type();

				// Rebuild ID map after sort
				for (size_t i = 0; i < particle_count(); i++) {
					auto id_ptr = this->template get_field_ptr<env::Field::id>(i);
					id_to_index_map[static_cast<size_t>(*id_ptr)] = i;
				}

				build_batches();
			}

			template<typename F>
				void for_each_interaction_batch(F && f) {
				// periodicity flags to jump table index
				const int mode = (flags.periodic_x ? 4 : 0) |
							  (flags.periodic_y ? 2 : 0) |
							  (flags.periodic_z ? 1 : 0);

				// trampoline lambda
				auto run_with_mode = [&]<bool PX, bool PY, bool PZ>() {
					auto bcp = [L=domain.extent](const vec3& dr) {
						return minimum_image<PX, PY, PZ>(dr, L);
					};

					for (const auto & batch : symmetric_batches) f(batch, bcp);
					for (const auto & batch : asymmetric_batches) f(batch, bcp);
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

			[[nodiscard]] std::vector<size_t> collect_indices_in_region(const env::Box & region) const {
				std::vector<size_t> ret;

				const double domain_vol = domain.volume();
				const auto intersection = domain.intersection(region);

				if (domain_vol > 1e-9 && intersection.has_value()) { // avoid division by zero
					const double ratio = intersection->volume() / domain_vol;

					// apply 1.1x safety factor because distributions are rarely perfectly uniform
					// clamp to particles.size() to avoid over-reservation
					const auto est = static_cast<size_t>(particle_count() * ratio * 1.1);
					ret.reserve(std::min(est, particle_count()));
				}

				for (size_t i = 0; i < particle_count(); i++) {
					auto particle = this->template view<env::Field::position | env::Field::state>(i);

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

			void sort_storage_by_type() {
				size_t n = this->particle_count();
				if (n == 0) return;

				std::vector<size_t> p(n);
				std::iota(p.begin(), p.end(), 0);

				std::sort(p.begin(), p.end(), [&](size_t i, size_t j) {
					auto t1 = *this->template get_field_ptr<env::Field::type>(i);
					auto t2 = *this->template get_field_ptr<env::Field::type>(j);
					return t1 < t2;
				});

				std::vector done(n, false);
				for (size_t i = 0; i < n; ++i) {
					if (done[i]) continue;
					size_t current = i;
					while (i != p[current]) {
						size_t next = p[current];
						swap_particles(current, next); // <--- Calls AoS or SoA swap!
						p[current] = current;
						current = next;
					}
					done[current] = true;
				}
			}

			void build_batches() {
				std::unordered_map<env::ParticleType, std::ranges::iota_view<size_t, size_t>> type_ranges;
				if (particle_count() == 0) return;

				size_t start = 0;
				auto current_type = *this->template get_field_ptr<env::Field::type>(0);

				for (size_t i = 0; i < particle_count(); i++) {
					auto type = *this->template get_field_ptr<env::Field::type>(i);

					if (current_type != type) {
						type_ranges[current_type] = std::ranges::iota_view(start, i);;
						start = i;
						current_type = type;
					}
				}
				type_ranges[current_type] = std::ranges::iota_view(start, particle_count()); // last batch

				const env::ParticleType n_types = current_type + 1;

				// create a batch for each particle type
				for (env::ParticleType type = 0; type < n_types; type++) {
					SymmetricBatch batch;
					batch.types = {type, type};
					batch.indices = type_ranges[type];
					symmetric_batches.push_back(batch);
				}

				// create a batch for each distinct pair of particle types
				for (env::ParticleType t1 = 0; t1 < n_types; t1++) {
					for (env::ParticleType t2 = t1 + 1; t2 < n_types; t2++) {
						AsymmetricBatch batch;
						batch.types = {t1, t2};
						batch.indices1 = type_ranges[t1];
						batch.indices2 = type_ranges[t2];
						asymmetric_batches.push_back(batch);
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

	namespace internal {
		template <class U>
		class DirectSumAoS final : public AoSContainer<container::DirectSumAoS, U> {
			using Base = AoSContainer<container::DirectSumAoS, U>;
			using typename Base::ParticleRecord;
			friend Base;
			using Base::particles;
			using Base::domain;
			using Base::flags;
			using Base::id_to_index_map;

			// TODO: once parallelization has been added change the policy to BatchParallelPolicy::inner
			struct AsymmetricBatch : SerialBatch<BatchSymmetry::Asymmetric>{
				std::ranges::iota_view<size_t, size_t> indices1{};
				std::ranges::iota_view<size_t, size_t> indices2{};
			};

			struct SymmetricBatch : SerialBatch<BatchSymmetry::Symmetric>{
				std::ranges::iota_view<size_t, size_t> indices {};
			};

		public:
			using Base::Base;

			void build(const std::vector<ParticleRecord> & particlesIn) {
				this->build_storage(particlesIn);

				// sort particles by type to create contiguous blocks
				std::sort(particles.begin(),particles.end(),
					[](const auto& a, const auto& b) {
						return a.type < b.type;
					});

				// Rebuild ID map after sort
				for (size_t i = 0; i < particles.size(); i++) {
					const auto id = static_cast<size_t>(particles[i].id);
					id_to_index_map[id] = i;
				}

				build_batches();
			}


			template<typename F>
			void for_each_interaction_batch(F && f) {
				// periodicity flags to jump table index
				const int mode = (flags.periodic_x ? 4 : 0) |
							  (flags.periodic_y ? 2 : 0) |
							  (flags.periodic_z ? 1 : 0);

				// trampoline lambda
				auto run_with_mode = [&]<bool PX, bool PY, bool PZ>() {
					auto bcp = [L=domain.extent](const vec3& dr) {
						return minimum_image<PX, PY, PZ>(dr, L);
					};

					for (const auto & batch : symmetric_batches) f(batch, bcp);
					for (const auto & batch : asymmetric_batches) f(batch, bcp);
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


			[[nodiscard]] std::vector<size_t> collect_indices_in_region(const env::Box & region) const {
				std::vector<size_t> ret;

				const double domain_vol = domain.volume();
				const auto intersection = domain.intersection(region);

				if (domain_vol > 1e-9 && intersection.has_value()) { // avoid division by zero
					const double ratio = intersection->volume() / domain_vol;

					// apply 1.1x safety factor because distributions are rarely perfectly uniform
					// clamp to particles.size() to avoid over-reservation
					const auto est = static_cast<size_t>(particles.size() * ratio * 1.1);
					ret.reserve(std::min(est, particles.size()));
				}

				for (size_t i = 0; i < particles.size(); i++) {
					if (region.contains(particles[i].position) && particles[i].state != env::ParticleState::DEAD) {
						ret.push_back(i);
					}
				}

				return ret;
			}

			void rebuild_structure() {}

		private:

			std::vector<SymmetricBatch> symmetric_batches;
			std::vector<AsymmetricBatch> asymmetric_batches;

			void build_batches() {
				std::unordered_map<env::ParticleType, std::ranges::iota_view<size_t, size_t>> type_ranges;
				if (particles.empty()) return;

				size_t start = 0;
				auto current_type = particles[0].type;

				for (size_t i = 0; i < particles.size(); i++) {
					auto & p = particles[i];
					if (current_type != p.type) {
						type_ranges[current_type] = std::ranges::iota_view(start, i);;

						start = i;
						current_type = p.type;
					}
				}
				type_ranges[current_type] = std::ranges::iota_view(start, particles.size()); // last batch

				const env::ParticleType n_types = current_type + 1;

				// create a batch for each particle type
				for (env::ParticleType type = 0; type < n_types; type++) {
					SymmetricBatch batch;
					batch.types = {type, type};
					batch.indices = type_ranges[type];
					symmetric_batches.push_back(batch);
				}

				// create a batch for each distinct pair of particle types
				for (env::ParticleType t1 = 0; t1 < n_types; t1++) {
					for (env::ParticleType t2 = t1 + 1; t2 < n_types; t2++) {
						AsymmetricBatch batch;
						batch.types = {t1, t2};
						batch.indices1 = type_ranges[t1];
						batch.indices2 = type_ranges[t2];
						asymmetric_batches.push_back(batch);
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
}
