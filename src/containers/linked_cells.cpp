

#include "april/containers/linked_cells.h"

#include <iostream>
#include <limits>
#include "april/env/interaction.h"
#include "april/env/particle.h"

namespace april::core {
	constexpr unsigned max_uint = std::numeric_limits<unsigned>::max();
	LinkedCells::LinkedCells(const double cell_size):
	grid_constant(cell_size),
	outside_cell(Cell::ParticleSet(0), {max_uint, max_uint, max_uint}, 0)
	{}

	void LinkedCells::build() {
		build_cells();
		build_cell_pairs();
	}

	void LinkedCells::build_cells() {
		grid_constant = std::max(interaction_manager->get_max_cutoff(), grid_constant);
		if (grid_constant <= 0) {
			grid_constant = std::max(extent[0], std::max(extent[1], extent[2]));
		}

		const auto num_x = static_cast<unsigned int>(std::max(1.0, floor(extent[0] / grid_constant)));
		const auto num_y = static_cast<unsigned int>(std::max(1.0, floor(extent[1] / grid_constant)));
		const auto num_z = static_cast<unsigned int>(std::max(1.0, floor(extent[2] / grid_constant)));

		cell_size = {extent[0] / num_x, extent[1] / num_y, extent[2] / num_z};
		cell_count = uint3{num_x, num_y, num_z};
		cells.reserve(num_x * num_y * num_z);

		// create cells
		auto N = static_cast<env::impl::ParticleID>(particles->size());
		for (unsigned int x = 0; x < num_x; x++) {
			for (unsigned int y = 0; y < num_y; y++) {
				for (unsigned int z = 0; z < num_z; z++) {
					Cell cell = {
						.particles = Cell::ParticleSet(N),
						.idx = {x,y,z},
						.id = x + y * num_x + z * num_x * num_y
					};
					cells.push_back(cell);
				}
			}
		}

		outside_cell.particles.set_capacity(N);

		// fill cells with particles
		for (auto & p : *particles) {
			if (p.state == env::ParticleState::DEAD) continue;
			get_cell(p.position).particles.insert(p.index);
		}
	}

	void LinkedCells::build_cell_pairs() {
		std::vector<int3> displacements;

		for (int dx = -1; dx <= 1; dx++) {
			for (int dy = -1; dy <= 1; dy++) {
				for (int dz = -1; dz <= 1; dz++) {
					if (dx == 0 && dy == 0 && dz == 0) continue;
					displacements.push_back({dx, dy, dz});
				}
			}
		}

		for (const auto d : displacements) {
			for (unsigned int z = 0; z < cell_count[2]; z++) {
				for (unsigned int y = 0; y < cell_count[1]; y++) {
					for (unsigned int x = 0; x < cell_count[0]; x++) {
						uint3 idx1 = {x,y,z};
						uint3 idx2 = {idx1[0] + d[0], idx1[1] + d[1], idx1[2] + d[2]};

						if (idx1 >= idx2 || // ensure only unique pairs are added
							idx2[0] < 0 || idx2[1] < 0 || idx2[2] < 0 ||
							idx2[0] >= cell_count[0] ||idx2[1] >= cell_count[1] || idx2[2] >= cell_count[2]
						) continue;
						
						cell_pairs.emplace_back(get_cell(idx1), get_cell(idx2));
					}
				}
			}
		}
	}

	LinkedCells::Cell & LinkedCells::get_cell(const vec3& position) noexcept {
		const vec3 pos = position - origin;
		if (pos[0] < 0 || pos[1] < 0 || pos[2] < 0) {
			return outside_cell;
		}

		// maybe we can somehow get around the division?
		const auto x = static_cast<uint32_t>(pos[0] / cell_size[0]);
		const auto y = static_cast<uint32_t>(pos[1] / cell_size[1]);
		const auto z = static_cast<uint32_t>(pos[2] / cell_size[2]);

		if (x >= cell_count[0] || y >= cell_count[1] || z >= cell_count[2]) {
			return outside_cell;
		}

		return get_cell(uint3{x,y,z});
	}

	LinkedCells::Cell& LinkedCells::get_cell(const uint3& idx) {
		const size_t id = idx[0] + idx[1] * cell_count[0] + idx[2] * cell_count[0] * cell_count[1];
		return cells[id];
	}

	void LinkedCells::calculate_forces() {

		for (auto & p : *particles) {
			p.reset_force();

			Cell & old_cell = get_cell(p.old_position);
			Cell & new_cell = get_cell(p.position);

			if (old_cell.idx != new_cell.idx) {
				old_cell.particles.erase(p.index);
				new_cell.particles.insert(p.index);
			}
		}

		// go through every cell and apply direct sum
		for (Cell & cell : cells) {
			for (size_t i = 0; i < cell.particles.size(); i++) {
				for (size_t j = i + 1; j < cell.particles.size(); j++) {
					auto & p1 = (*particles)[cell.particles[i]];
					auto & p2 = (*particles)[cell.particles[j]];

					const vec3 force = interaction_manager->evaluate(p1, p2);
					p1.force += force;
					p2.force -= force;
				}
			}
		}

		//go through all cell pairs
		for (auto & [c1, c2] : cell_pairs) {
			for (const unsigned int i : c1.particles) {
				for (const unsigned int j : c2.particles) {
					auto & p1 = (*particles)[i];
					auto & p2 = (*particles)[j];

					const vec3 force = interaction_manager->evaluate(p1, p2);
					p1.force += force;
					p2.force -= force;
				}
			}
		}
	}

}
