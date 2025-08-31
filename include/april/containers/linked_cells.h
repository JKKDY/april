#pragma once

#include "april/containers/contiguous_container.h"
#include "april/utils/set.hpp"
#include "april/env/particle.h"



namespace april::cont {
	namespace impl {
		template <class Env> class LinkedCells;
	}

	template<typename T> concept IsParticleSet = requires(T set, size_t id) {
		// Construction with maxId
		{ T(std::declval<size_t>()) };

		// Insert method
		{ set.insert(id) } -> std::same_as<void>;

		// Erase method
		{ set.erase(id) } -> std::same_as<void>;

		// Contains method
		{ set.contains(id) } -> std::same_as<bool>;

		// Iteration methods
		{ set.begin() } -> std::input_iterator;
		{ set.end() } -> std::input_iterator;

		// Size method
		{ set.size() } -> std::same_as<size_t>;
	};


	struct LinkedCells {
		template<typename  Env> using impl = impl::LinkedCells<Env>;
		double cell_size_hint;
	};


	namespace impl {
		template <class Env>
		class LinkedCells final : public ContiguousContainer<cont::LinkedCells, Env> {

			static constexpr unsigned max_uint = std::numeric_limits<unsigned>::max();

			struct Cell {
				using ParticleSet = utils::IndexSet<env::impl::ParticleID>;
				static_assert(IsParticleSet<ParticleSet>, "ParticleSet must implement IsParticleSet interface");

				ParticleSet particles;
				uint3 idx;
				unsigned int id {};
			};

			struct CellPair {
				Cell& first;
				Cell& second;
			};

			using Base = ContiguousContainer<cont::LinkedCells, Env>;
			using typename Base::Particle;
			using typename Base::ParticleID;
			using Base::interactions;
			using Base::cfg;
			using Base::domain;
			using Base::particles;
		public:
			explicit LinkedCells(const cont::LinkedCells& config)
			  : Base(config)
				// outside_cell{ Cell::ParticleSet, uint3{max_uint, max_uint, max_uint}, 0u }
			{
				outside_cell.id = 0;
				outside_cell.idx = uint3{max_uint, max_uint, max_uint};
			}

			void build(const std::vector<Particle>& particles) {
				this->build_storage(particles);
				std::ranges::sort(this->particles,
								  [](const Particle& a, const Particle& b) {
									  return a.id < b.id;
								  });
				build_cells();
				build_cell_pairs();
			}

			void calculate_forces() {
				for (auto & p : particles) {
					p.reset_force();

					Cell & old_cell = get_cell(p.old_position);
					Cell & new_cell = get_cell(p.position);

					if (old_cell.idx != new_cell.idx) {
						old_cell.particles.erase(p.id);
						new_cell.particles.insert(p.id);
					}
				}

				// go through every cell and apply direct sum
				for (Cell & cell : cells) {
					for (size_t i = 0; i < cell.particles.size(); i++) {
						for (size_t j = i + 1; j < cell.particles.size(); j++) {
							auto & p1 = particles[cell.particles[i]];
							auto & p2 = particles[cell.particles[j]];

							const vec3 force = interactions->evaluate(p1, p2);
							p1.force += force;
							p2.force -= force;
						}
					}
				}

				//go through all cell pairs
				for (auto & [c1, c2] : cell_pairs) {
					for (const unsigned int i : c1.particles) {
						for (const unsigned int j : c2.particles) {
							auto & p1 = particles[i];
							auto & p2 = particles[j];

							const vec3 force = interactions->evaluate(p1, p2);
							p1.force += force;
							p2.force -= force;
						}
					}
				}
			}
		private:
			void build_cells() {
				cfg.cell_size_hint = std::max(interactions->get_max_cutoff(), cfg.cell_size_hint);
				if (cfg.cell_size_hint <= 0) {
					cfg.cell_size_hint = domain.extent.max();
				}

				const auto num_x = static_cast<unsigned int>(std::max(1.0, floor(domain.extent[0] / cfg.cell_size_hint)));
				const auto num_y = static_cast<unsigned int>(std::max(1.0, floor(domain.extent[1] / cfg.cell_size_hint)));
				const auto num_z = static_cast<unsigned int>(std::max(1.0, floor(domain.extent[2] / cfg.cell_size_hint)));

				cell_size = {domain.extent[0] / num_x, domain.extent[1] / num_y, domain.extent[2] / num_z};
				inv_cell_size = {1/cell_size[0], 1/cell_size[1], 1/cell_size[2]};
				cell_count = uint3{num_x, num_y, num_z};
				cells.reserve(num_x * num_y * num_z);

				// create cells
				const auto N = static_cast<env::impl::ParticleID>(particles.size());
				for (unsigned int x = 0; x < num_x; x++) {
					for (unsigned int y = 0; y < num_y; y++) {
						for (unsigned int z = 0; z < num_z; z++) {
							Cell cell = {
								.particles =  utils::IndexSet(N),
								.idx = {x,y,z},
								.id = x + y * num_x + z * num_x * num_y
							};
							cells.push_back(cell);
						}
					}
				}

				outside_cell.particles.set_capacity(N);

				// fill cells with particles
				for (auto & p : particles) {
					if (p.state == env::ParticleState::DEAD) continue;
					get_cell(p.position).particles.insert(p.id);
				}
			}

			void build_cell_pairs() {
				static const int3 displacements[13] = {
					{ 1, 0, 0}, { 0, 1, 0}, { 0, 0, 1},
					{ 1, 1, 0}, { 1,-1, 0}, { 1, 0, 1},
					{-1, 0, 1}, { 0, 1, 1}, { 0,-1, 1},
					{ 1, 1, 1}, { 1,-1, 1}, {-1, 1, 1},
					{-1,-1, 1}
				};

				for (const auto d : displacements) {
					for (unsigned int z = 0; z < cell_count[2]; z++) {
						for (unsigned int y = 0; y < cell_count[1]; y++) {
							for (unsigned int x = 0; x < cell_count[0]; x++) {
								const int3 base{static_cast<int>(x), static_cast<int>(y), static_cast<int>(z)};
								const int3 n = base + d;

								if (n.x < 0 || n.y < 0 || n.z < 0)
									continue;
								if (n.x >= static_cast<int>(cell_count[0]) || n.y >= static_cast<int>(cell_count[1]) || n.z >= static_cast<int>(cell_count[2]))
									continue;

								uint3 idx1 = {x,y,z};
								uint3 idx2{static_cast<uint32_t>(n.x), static_cast<uint32_t>(n.y), static_cast<uint32_t>(n.z)};

								// if (idx1 >= idx2) // ensure only unique pairs are added
								// 	continue;

								cell_pairs.emplace_back(get_cell(idx1), get_cell(idx2));
							}
						}
					}
				}
			}

			Cell& get_cell(const vec3& position) noexcept {
				const vec3 pos = position - domain.origin;
				if (pos[0] < 0 || pos[1] < 0 || pos[2] < 0) {
					return outside_cell;
				}

				// maybe we can somehow get around the division?
				const auto x = static_cast<uint32_t>(pos[0] * inv_cell_size[0]);
				const auto y = static_cast<uint32_t>(pos[1] * inv_cell_size[1]);
				const auto z = static_cast<uint32_t>(pos[2] * inv_cell_size[2]);

				if (x >= cell_count[0] || y >= cell_count[1] || z >= cell_count[2]) {
					return outside_cell;
				}

				return get_cell(uint3{x,y,z});
			}

			Cell& get_cell(const uint3& idx) {
				const size_t id = idx[0] + idx[1] * cell_count[0] + idx[2] * cell_count[0] * cell_count[1];
				return cells[id];
			}

			vec3 cell_size;
			vec3 inv_cell_size;
			uint3 cell_count{};
			Cell outside_cell;
			std::vector<Cell> cells;
			std::vector<CellPair> cell_pairs;
		};
	}
}