#pragma once

#include <vector>
#include <cmath>
#include <algorithm>
#include <utility>
#include <functional>
#include <limits>

#include "lc_batching.hpp"
#include "april/base/types.hpp"
#include "april/containers/batching/topology_batch.hpp"
#include "april/particle/properties.hpp"
#include "april/core/domain.hpp"
#include "april/exec/particle_kernel.hpp"
#include "april/exec/parallel_utils.hpp"
#include "april/exec/executors/context.hpp"

namespace april::math {
	struct Range;
}

namespace april::container::internal {
	template <class Base>
	class LinkedCellsCore : public Base {
		using cell_index_t = uint32_t;

		enum CellWrapFlag : uint8_t {
			NO_WRAP = 0,
			WRAP_X=1,
			WRAP_Y=2,
			WRAP_Z=4,
		};

		struct CellPair {
			const cell_index_t c1 = {};
			const cell_index_t c2 = {};
		};

		struct WrappedCellPair {
			const cell_index_t c1 = {};
			const cell_index_t c2 = {};
			const CellWrapFlag force_wrap = {};
			const vec3 shift = {};
		};

	public:
		using Base::Base;
		using Base::build_storage;
		using typename Base::ParticleRecord;
		using Base::vector_policy;
		using Base::parallel_policy;

		void build (this auto&& self, const std::vector<ParticleRecord>& particles) {
			self.setup_topology_batches();
			self.setup_cell_grid();
			self.init_cell_order();
			self.create_neighbor_stencil();
			self.compute_wrapped_cell_pairs();
			self.build_storage(particles);
			self.pre_allocate_assignment_bins();
			self.rebuild_structure_impl();
			self.schedule_phases();

			// TODO once we can use reflection automatically serialize add a vec3 last_rebuild_position member to particle attributes
			self.last_x.resize(particles.size());
			self.last_y.resize(particles.size());
			self.last_z.resize(particles.size());

			self.template for_each_particle<ParallelPolicy::Serial>(
				scalar_kernel<ParticleField::position | ParticleField::id>([&](auto && p) {
					self.last_x[p.id] = p.position.x;
					self.last_y[p.id] = p.position.y;
					self.last_z[p.id] = p.position.z;
				})
			);



		}

		template<ParallelPolicy P, typename Func>
		void for_each_topology_batch(this auto&& self, Func && f) {
			for (const auto& phase : self.topology_phases) {
				self.thread_executor.template execute<P>(phase.size(), [&](size_t i) {
					f(phase[i]);
				});
			}
		}


		struct alignas(64) PaddedThreadBuffer {
			// Stores {bin_index, particle_index}
			std::vector<std::pair<uint32_t, uint32_t>> records;
		};

		// Persistent member variable
		std::vector<PaddedThreadBuffer> thread_local_buffers;

		std::vector<vec3::type> last_x;
		std::vector<vec3::type> last_y;
		std::vector<vec3::type> last_z;

		void rebuild_structure(this auto && self) {
			std::atomic rebuild = false;

			// check if a particle has moved further than the  skin thickness
			self.template for_each_particle<parallel_policy>(
				scalar_kernel<ParticleField::position| ParticleField::id>([&](auto && p) {
					if (std::abs(self.last_x[p.id] - p.position.x) > self.verlet_skin / 2) rebuild = true;
					if (std::abs(self.last_y[p.id] - p.position.y) > self.verlet_skin / 2) rebuild = true;
					if (std::abs(self.last_z[p.id] - p.position.z) > self.verlet_skin / 2) rebuild = true;
				})
			);

			if (rebuild) {
				self.rebuild_structure_impl();

				// cache current particles positions
				self.template for_each_particle<parallel_policy>(
					scalar_kernel<ParticleField::position | ParticleField::id>([&](auto && p) {
						self.last_x[p.id] = p.position.x;
						self.last_y[p.id] = p.position.y;
						self.last_z[p.id] = p.position.z;
					})
				);
			}
		}

		void rebuild_structure_impl(this auto&& self) {
			self.reorder_storage(self.n_bins, [&](const size_t i) {
				const auto p = self.template view<ParticleField::position | ParticleField::type>(i);
				const size_t cid = self.cell_index_from_position( p.position);
				return self.bin_index(cid, p.type);
			});
		}

		[[nodiscard]] std::vector<size_t> collect_indices_in_region(this const auto& self, const core::Box & region) {
		    const std::vector<cell_index_t> cells = self.get_cells_in_region(region);

		    if (cells.empty()) return {};

			// partition the entire particle range into independent blocks/tasks
		    auto blocks = exec::make_linear_schedule(math::Range{0, cells.size()}, self.linear_schedule_config);
		    const size_t num_tasks = blocks.size();

			// allocate buffers for each task
		    std::vector<std::vector<size_t>> local_results(self.thread_executor.num_threads());

			// preallocate storage with heuristic (assumes uniform distribution)
		    const size_t est_count_per_task = (self.particle_count() * cells.size() / self.n_cells) / num_tasks + 1;
		    for (auto& loc : local_results) {
		        loc.reserve(est_count_per_task);
		    }

			// process all tasks in parallel
		    self.thread_executor.template execute<parallel_policy>(num_tasks, [&](const size_t t_idx) {
		        const auto& block = blocks[t_idx];
		        auto& local_ret = local_results[exec::thread_index()];

		        // each task loops over its assigned chunk of cells
		        for (size_t c = block.start; c < block.stop; ++c) {
		        	// get the search range (slice of particle data)
		            const size_t cid = cells[c];
		            const auto [start_idx, end_idx] = self.cell_index_range(cid);

		        	// skip empty ranges
		            if (start_idx == end_idx) continue;

		        	// kernel checks if a particle is alive and inside the region
		        	// this is scalar because SoA layout is tightly packed resulting in spilling
		        	// TODO: make this a universal kernel by propagating masks along with particle data
		        	auto kernel = april::scalar_kernel<ParticleField::position | ParticleField::state> (
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

		            // run the kernel
		            self.for_each_particle(start_idx, end_idx, kernel);
		        }
		    });

			// allocate storage for indices buffer
		    size_t total_found = 0;
		    for (const auto& loc : local_results) {
		        total_found += loc.size();
		    }

		    std::vector<size_t> ret;
		    ret.reserve(total_found);

			// merge all local buffers into indices buffer
		    for (const auto& loc : local_results) {
		        ret.insert(ret.end(), loc.begin(), loc.end());
		    }

		    return ret;
		}


	protected:
		size_t outside_cell_id {};
		size_t n_grid_cells {};
		size_t n_cells {}; // total cells = grid + outside
		size_t n_types {}; // types range from 0 ... n_types-1
		size_t n_bins {}; // number of bins (cells * types)
		double global_cutoff {}; // maximum force cutoff
		double verlet_skin {};

		vec3d cell_size; // side lengths of each cell
		vec3d inv_cell_size; // cache the inverse of each size component to avoid divisions
		uint3 cells_per_axis{}; // number of cells along each axis
		uint3::type cell_per_axis_xy{}; // = cells_per_axis.x * cells_per_axis.y

		std::vector<std::vector<size_t>> bin_assignments;
		std::vector<cell_index_t> cell_ordering; // map x,y,z flat index (Nx*Ny*z+Nx*y+x) to ordering index

		// cell pair info
		std::vector<int3> neighbor_stencil;
		std::vector<WrappedCellPair> wrapped_cell_pairs;
		std::vector<std::vector<uint3>> phase_schedule; // for user defined coloring scheme
		std::vector<std::vector<WrappedCellPair>> wrapped_phase_schedule;


		//------
		// SETUP
		//------
		void setup_topology_batches() {
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

		void setup_cell_grid(this auto&& self) {
			// determine the physical cutoff (max_rc) from interactions
			double max_cutoff = 0;
			for (const auto & interaction : self.interaction_map.interactions) {
				if (interaction.is_active && !interaction.used_by_types.empty() && interaction.cutoff > max_cutoff) {
					max_cutoff = interaction.cutoff;
				}
			}

			// handle edge cases (no interactions or cutoff > domain)
			// ensure at least a 2x2x2 grid (before applying cell size factor (CSF))
			const vec3::type min_dim = vec3 {
				self.domain.extent.x == 0 ? std::numeric_limits<double>::max() : self.domain.extent.x,
				self.domain.extent.y == 0 ? std::numeric_limits<double>::max() : self.domain.extent.y,
				self.domain.extent.z == 0 ? std::numeric_limits<double>::max() : self.domain.extent.z
			}.min();
			if (max_cutoff <= 0 || max_cutoff > min_dim / 2.0) {
				// TODO log warning
				max_cutoff = min_dim / 2.0;
			}

			// set target cell size and verlet skin
			double target_cell_size = self.config.get_width(max_cutoff);
			APRIL_ASSERT(target_cell_size > 0, "Calculated cell size must be > 0");

			self.verlet_skin = self.config.get_skin(target_cell_size);
			APRIL_ASSERT(target_cell_size >= 0, "Calculated cell size must be > 0");

			target_cell_size += self.verlet_skin;

			// compute number of cells along each axis
			// std::floor ensures that the resulting cells are larger than or equal to target_cell_size
			const auto num_x = static_cast<cell_index_t>(std::max(2.0, std::floor(self.domain.extent.x / target_cell_size)));
			const auto num_y = static_cast<cell_index_t>(std::max(2.0, std::floor(self.domain.extent.y / target_cell_size)));
			const auto num_z = static_cast<cell_index_t>(std::max(2.0, std::floor(self.domain.extent.z / target_cell_size)));

			// calculate cell size along (stretches to fit domain exactly)
			self.cell_size = {
				self.domain.extent.x / num_x,
				self.domain.extent.y / num_y,
				self.domain.extent.z / num_z
			};

			// cache inverse (useful for fast binning: index = coord * inv_cell_size)
			self.inv_cell_size = {
				self.cell_size.x > 0 ? 1.0/self.cell_size.x : 0.0,
				self.cell_size.y > 0 ? 1.0/self.cell_size.y : 0.0,
				self.cell_size.z > 0 ? 1.0/self.cell_size.z : 0.0
			  };

			self.cells_per_axis = uint3{num_x, num_y, num_z};
			self.cell_per_axis_xy = self.cells_per_axis.x * self.cells_per_axis.y;

			// set scalars
			self.n_types = self.interaction_map.types.size();
			self.n_grid_cells = num_x * num_y * num_z;
			self.n_cells = self.n_grid_cells + 1;
			self.outside_cell_id = self.n_grid_cells;
			self.global_cutoff = max_cutoff;

			// allocate buffers
			self.n_bins = self.n_cells * self.n_types;
			self.bin_starts.resize(self.n_bins);
			self.bin_assignments.resize(self.n_bins);
		}

		void init_cell_order(this auto && self) {
			if (self.config.cell_ordering_fn.has_value()) {
				const std::function<std::vector<uint32_t>(uint3)> & cell_ordering_fn = self.config.cell_ordering_fn.value();
				self.cell_ordering = cell_ordering_fn(self.cells_per_axis);
			}
		}

		void create_neighbor_stencil(this auto && self) {
			// number of cells within a box of sidelengths global_cutoff
			// ceil ensures we really check every relevant cell
			const auto nx = static_cast<int>(std::ceil(self.global_cutoff * self.inv_cell_size.x));
			const auto ny = static_cast<int>(std::ceil(self.global_cutoff * self.inv_cell_size.y));
			const auto nz = static_cast<int>(std::ceil(self.global_cutoff * self.inv_cell_size.z));

			const double cutoff_sq = self.global_cutoff * self.global_cutoff;

			// loop over every cell within a 2*rc sized box
			for (int z = 0; z <= nz; ++z) {  // we only need a half sphere so we can exclude cells with z<0
				for (int y = -ny; y <= ny; ++y) {
					for (int x = -nx; x <= nx; ++x) {

						// half sphere filter: only "forward" cells (use tuple ordering for filtering)
						if (std::make_tuple(z, y, x) <= std::make_tuple(0, 0, 0)) {
							continue;
						}

						// calculate distance between cell at x,y,z and center cell
						vec3d dist_vec = {
							std::abs(x) > 1 ? (std::abs(x) - 1) * self.cell_size.x : 0,
							std::abs(y) > 1 ? (std::abs(y) - 1) * self.cell_size.y : 0,
							std::abs(z) > 1 ? (std::abs(z) - 1) * self.cell_size.z : 0,
						};

						if (dist_vec.norm_squared() <= cutoff_sq) {
							self.neighbor_stencil.emplace_back(x,y,z);
						}
					}
				}
			}
		}

		void compute_wrapped_cell_pairs(this auto && self) {
			// check if cell index is out of bounds, wrap it to the other side and compute the spatial shift vector
			auto try_wrap_cell = [&](int3 & n, vec3 & shift, int ax) -> CellWrapFlag {
				const int dim_cells = static_cast<int>(self.cells_per_axis[ax]);

				if (n[ax] < 0) {
					// Wrapped "Left": Shift index by +Size
					n[ax] += dim_cells;
					shift[ax] = -self.domain.extent[ax];
				}
				else if (n[ax] >= dim_cells) {
					// Wrapped "Right": Shift index by -Size
					n[ax] -= dim_cells;
					shift[ax] = self.domain.extent[ax];
				}
				else {
					return NO_WRAP;
				}

				return static_cast<CellWrapFlag>(1 << ax); // maps axis to appropriate CellWrap flag
			};

			for (unsigned int z = 0; z < self.cells_per_axis.z; z++) {
				for (unsigned int y = 0; y < self.cells_per_axis.y; y++) {
					for (unsigned int x = 0; x < self.cells_per_axis.x; x++) {
						for (const auto displacement : self.neighbor_stencil) {
							if (displacement == int3{0,0,0}) continue;

							const int3 base {  // coordinates of cell 1
								static_cast<int>(x),
								static_cast<int>(y),
								static_cast<int>(z)
							};

							int3 n = base + displacement; // coordinates of cell 2
							vec3 shift = {};
							int8_t wrap_flags = {};

							// if periodic check if cell pair needs to be wrapped
							if (self.flags.periodic_x) wrap_flags |= try_wrap_cell(n, shift, 0);
							if (self.flags.periodic_y) wrap_flags |= try_wrap_cell(n, shift, 1);
							if (self.flags.periodic_z) wrap_flags |= try_wrap_cell(n, shift, 2);

							if (shift == vec3{}) continue;

							// if cell pair is not wrapped and cell 2 is outside of domain: skip pair
							if (n.x < 0 || n.y < 0 || n.z < 0 ||
								n.x >= static_cast<int>(self.cells_per_axis.x) ||
								n.y >= static_cast<int>(self.cells_per_axis.y) ||
								n.z >= static_cast<int>(self.cells_per_axis.z)
							) continue;

							self.wrapped_cell_pairs.emplace_back(
								self.cell_pos_to_idx(x,y,z),
								self.cell_pos_to_idx(n.x, n.y, n.z),
								static_cast<CellWrapFlag>(wrap_flags),
								shift
							);
						}
					}
				}
			}
		}

		void pre_allocate_assignment_bins() {
			const size_t num_bins = n_types * n_cells;
			bin_assignments.resize(num_bins);

			// for each bin assume somewhat uniform distribution + 50% buffer
			const size_t est_per_bin = static_cast<size_t>((this->particle_count() / std::max<size_t>(1, num_bins)) * 1.5);
			for(auto& bin : bin_assignments) {
				bin.reserve(est_per_bin);
			}
		}

		void schedule_phases() {
			const auto& batch_dim = this->config.block_size;

			// schedule phases for neighboring cells according to user defined color scheme
			for_each_block([&](size_t bx, size_t by, size_t bz) {
				const size_t logical_x = bx / batch_dim.x;
				const size_t logical_y = by / batch_dim.y;
				const size_t logical_z = bz / batch_dim.z;

				const size_t color = this->config.schedule_phases(logical_x, logical_y, logical_z, batch_dim);
				if (color >= phase_schedule.size()) {
					phase_schedule.resize(color + 1);
				}

				phase_schedule[color].emplace_back(bx, by, bz);
			});

			// schedule wrapped neighbor cells (if force wrapping is enabled)
			// create an edge list graph of interacting wrapped pairs
			std::vector<std::vector<uint32_t>> touched_cells;
			touched_cells.reserve(wrapped_cell_pairs.size());
			for (const auto& pair : wrapped_cell_pairs) {
				touched_cells.push_back({pair.c1, pair.c2});
			}

			// and build an intersection graph (Adjacency list of conflicting cell pairs)
			const auto conflict_graph = utility::graph::build_intersection_graph(touched_cells);

			// create a sequential processing order (indices 0 to N-1)
			std::vector<size_t> processing_order(wrapped_cell_pairs.size());
			std::iota(processing_order.begin(), processing_order.end(), 0);

			// partition into independent sets (Conflict-free phases)
			const auto independent_phases = utility::graph::greedy_independent_partitions(
				conflict_graph,
				processing_order
			);

			// Build the final schedule for the execution loop
			wrapped_phase_schedule.clear();
			wrapped_phase_schedule.reserve(independent_phases.size());

			for (const auto& phase_indices : independent_phases) {
				std::vector<WrappedCellPair> phase_pairs;
				phase_pairs.reserve(phase_indices.size());

				// Map the indices back to the actual WrappedCellPair objects
				for (size_t idx : phase_indices) {
					phase_pairs.push_back(wrapped_cell_pairs[idx]);
				}

				wrapped_phase_schedule.push_back(std::move(phase_pairs));
			}
		}


        // -----------------
        // LOOP ABSTRACTIONS
        // -----------------
        // Iterates over spatial blocks (cache blocking)
        template <typename Func>
        APRIL_FORCE_INLINE void for_each_block(Func&& fn) const {
            const auto& batch_dim = this->config.block_size;
            for (size_t bz = 0; bz < cells_per_axis.z; bz += batch_dim.z)
                for (size_t by = 0; by < cells_per_axis.y; by += batch_dim.y)
                    for (size_t bx = 0; bx < cells_per_axis.x; bx += batch_dim.x)
                        fn(bx, by, bz);
        }

		// Iterates over all unique pairs of types (T1, T2) where T2 >= T1
		template <typename Func>
		APRIL_FORCE_INLINE void for_each_type_pair(Func&& fn) const {
			for (size_t t1 = 0; t1 < this->n_types; ++t1)
				for (size_t t2 = t1; t2 < this->n_types; ++t2)
					fn(t1, t2);
		}

        // Iterates over the cells inside a specific block
        template <typename Func>
        APRIL_FORCE_INLINE void for_each_cell_in_block(size_t bx, size_t by, size_t bz, Func&& fn) const {
            const auto& bdim = this->config.block_size;

            // Calculate limits (handling edge blocks that might be smaller)
            const size_t z_end = std::min(bz + bdim.z, static_cast<size_t>(cells_per_axis.z));
            const size_t y_end = std::min(by + bdim.y, static_cast<size_t>(cells_per_axis.y));
            const size_t x_end = std::min(bx + bdim.x, static_cast<size_t>(cells_per_axis.x));

            for (size_t z = bz; z < z_end; ++z)
                for (size_t y = by; y < y_end; ++y)
                    for (size_t x = bx; x < x_end; ++x)
                        fn(x, y, z);
        }



		// ----------------
		// BATCH ITERATIONS
		// ----------------
		template <typename GetRange, typename AddSym, typename AddAsym>
		APRIL_FORCE_INLINE void process_cell_interactions(
			size_t x, size_t y, size_t z,
			size_t t1, size_t t2,
			GetRange&& get_range,
			AddSym&& add_sym,
			AddAsym&& add_asym
		) const {
			const size_t c = this->cell_pos_to_idx(x, y, z);
			auto range1 = get_range(c, t1);

			// intra-cell: process forces between particles inside the cell
			if (t1 == t2) {
				if (range1.size() > 1) add_sym(range1);
			} else {
				auto range2 = get_range(c, t2);
				if (range1.size() > 0 && range2.size() > 0) {
					add_asym(range1, range2);
				}
			}

			// If T1==T2 and Cell(T1) is empty, neighbors don't matter.
			// For mixed types (T1!=T2), we continue because we need the reverse check (Neighbor(T1) vs Cell(T2)).
			if (range1.empty() && (t1 == t2)) return;

			// inter-cell: process forces between particles of neighboring cells
			for (auto offset : this->neighbor_stencil) {
				size_t c_n = this->get_neighbor_idx(x, y, z, offset);
				if (c_n == this->outside_cell_id) continue;

				// Interaction 1: Cell(T1) -> Neighbor(T2)
				auto range_n2 = get_range(c_n, t2);
				if (!range1.empty() && !range_n2.empty()) {
					add_asym(range1, range_n2);
				}

				// Interaction 2: Neighbor(T1) -> Cell(T2)
				// (due to the half stencil we would otherwise never iterate over this combination)
				if (t1 != t2) {
					auto range2 = get_range(c, t2);
					auto range_n1 = get_range(c_n, t1);

					if (!range2.empty() && !range_n1.empty()) {
						add_asym(range_n1, range2);
					}
				}
			}
		}

		template <ParallelPolicy P, typename Func, typename GetIndices, typename ProcessBatch>
		APRIL_FORCE_INLINE void for_each_wrapped_interaction(
			this const auto& self,
			Func&& func,
			GetIndices&& get_indices,
			ProcessBatch&& process_batch
		) {
			if (self.wrapped_phase_schedule.empty()) return;

			// Iterate sequentially over the independent phases
			for (const auto& phase : self.wrapped_phase_schedule) {

				// Execute the pairs within the phase in parallel
				self.thread_executor.template execute<P>(phase.size(), [&](size_t phase_idx) {
					const auto& pair = phase[phase_idx];
					auto bcp = [&pair](const auto& diff) { return diff + pair.shift; };

					for (size_t t1 = 0; t1 < self.n_types; ++t1) {
						auto range1 = get_indices(pair.c1, t1);
						if (range1.empty()) continue;

						for (size_t t2 = 0; t2 < self.n_types; ++t2) {
							auto range2 = get_indices(pair.c2, t2);
							if (range2.empty()) continue;

							process_batch(func, range1, range2, t1, t2, bcp);
						}
					}
				});
			}
		}



		//----------
		// UTILITIES
		//----------
		[[nodiscard]] APRIL_FORCE_INLINE size_t get_neighbor_idx(const size_t x, const size_t y, const size_t z, const int3 offset) const {
			const int nx = static_cast<int>(x) + offset.x;
			const int ny = static_cast<int>(y) + offset.y;
			const int nz = static_cast<int>(z) + offset.z;

			// Note: If you implement periodic BCs via ghost cells, this check might change.
			// For now, it clamps to "Outside".
			if (nx < 0 || ny < 0 || nz < 0 ||
				nx >= static_cast<int>(this->cells_per_axis.x) ||
				ny >= static_cast<int>(this->cells_per_axis.y) ||
				nz >= static_cast<int>(this->cells_per_axis.z)) {
				return this->outside_cell_id;
				}
			return this->cell_pos_to_idx(nx, ny, nz);
		}

		// gather all cell ids whose cells have an intersection with the box region
		[[nodiscard]] std::vector<cell_index_t> get_cells_in_region(this const auto& self, const core::Box & box) {
			//  Convert world coords to cell coords (relative to domain origin)
			const vec3d min = (box.min - self.domain.min) * self.inv_cell_size;
			const vec3d max = (box.max - self.domain.min) * self.inv_cell_size;

			// clamp cell coordinates to valid ranges
			const vec3d min_clamped = {
				std::clamp(std::floor(min.x), 0.0, static_cast<double>(self.cells_per_axis.x - 1)),
				std::clamp(std::floor(min.y), 0.0, static_cast<double>(self.cells_per_axis.y - 1)),
				std::clamp(std::floor(min.z), 0.0, static_cast<double>(self.cells_per_axis.z - 1))
			};

			const vec3d max_clamped = {
				std::clamp(std::ceil(max.x), 0.0, static_cast<double>(self.cells_per_axis.x - 1)),
				std::clamp(std::ceil(max.y), 0.0, static_cast<double>(self.cells_per_axis.y - 1)),
				std::clamp(std::ceil(max.z), 0.0, static_cast<double>(self.cells_per_axis.z - 1))
			};

			// find the lowest left cell
			const uint3 min_cell = {
				static_cast<uint3::type>(min_clamped.x),
				static_cast<uint3::type>(min_clamped.y),
				static_cast<uint3::type>(min_clamped.z)
			};

			// find the highest right cell
			const uint3 max_cell = {
				static_cast<uint3::type>(max_clamped.x),
				static_cast<uint3::type>(max_clamped.y),
				static_cast<uint3::type>(max_clamped.z)
			};

			// add all cells in that region
			const auto cell_counts = max_cell - min_cell;
			std::vector<cell_index_t> cells;
			cells.reserve(cell_counts.x * cell_counts.y * cell_counts.z);
			for (size_t x = min_cell.x; x <= max_cell.x; ++x) {
				for (size_t y = min_cell.y; y <= max_cell.y; ++y) {
					for (size_t z = min_cell.z; z <= max_cell.z; ++z) {
						cells.push_back(self.cell_pos_to_idx(x,y,z));
					}
				}
			}

			if (!(box.min>= self.domain.min && box.max <= self.domain.max)) {
				cells.push_back(self.outside_cell_id);
			}

			return cells;
		}

		[[nodiscard]] size_t bin_index(const size_t cell_id, const ParticleType type = 0) const {
			return cell_id * n_types + static_cast<size_t>(type);
		}

		[[nodiscard]] std::pair<size_t, size_t> cell_index_range(this const auto & self, const uint32_t cid) {
			const size_t start_bin_idx = self.bin_index(cid);
			size_t start = self.bin_starts[start_bin_idx];
			size_t end = start_bin_idx + self.n_types >= self.bin_starts.size() ? self.capacity() : self.bin_starts[start_bin_idx + self.n_types];
			return {start, end};

		}

		[[nodiscard]] uint32_t cell_pos_to_idx(this const auto & self, const uint32_t x, const uint32_t y, const uint32_t z) noexcept{
			const uint32_t flat_idx = z * self.cell_per_axis_xy + y * self.cells_per_axis.x + x;
			return self.cell_ordering.empty() ? flat_idx : self.cell_ordering[flat_idx];
		}

		uint32_t cell_index_from_position(this const auto & self, const vec3 & position) noexcept {
			const vec3 pos = position - self.domain.min;

			if ((pos.x < 0.0) | (pos.y < 0.0) | (pos.z < 0.0)) {
				return self.outside_cell_id;
			}

			const auto x = static_cast<uint32_t>(pos.x * self.inv_cell_size.x);
			const auto y = static_cast<uint32_t>(pos.y * self.inv_cell_size.y);
			const auto z = static_cast<uint32_t>(pos.z * self.inv_cell_size.z);

			if ((x >= self.cells_per_axis.x) | (y >= self.cells_per_axis.y) | (z >= self.cells_per_axis.z)) {
				return self.outside_cell_id;
			}

			return self.cell_pos_to_idx(x, y, z);
		}
	private:
		std::vector<std::vector<batching::TopologyBatch<LinkedCellsCore>>> topology_phases;
	};
}

