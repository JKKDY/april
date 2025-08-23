#pragma once

#include "april/algo/algorithm.h"
#include "april/utils/set.hpp"
#include "april/env/particle.h"



namespace april::algo {
	namespace impl {
		class LinkedCells;
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
		using impl = impl::LinkedCells;
		double cell_size_hint;
	};


	namespace impl {
		class LinkedCells final : public Algorithm<algo::LinkedCells> {
			struct Cell {
				using ParticleSet = utils::IndexSet<env::impl::ParticleID>;
				static_assert(IsParticleSet<ParticleSet>, "ParticleSet must implement IsParticleSet interface");

				ParticleSet particles;
				uint3 idx;
				unsigned int id;
			};

			struct CellPair {
				Cell& first;
				Cell& second;
			};
		public:
			using Algorithm::Algorithm;

			void build(const std::vector<Particle>& particles) override;
			void calculate_forces() override;

			Particle & get_particle_by_id(ParticleID id) override;
			ParticleID id_start() override;
			ParticleID id_end() override;

			Particle & get_particle_by_index(size_t index) noexcept override;
			size_t index_start() override;
			size_t index_end() override;

		private:
			void build_cells();
			void build_cell_pairs();

			Cell & get_cell(const vec3& position) noexcept;
			Cell & get_cell(const uint3 & idx);

			vec3 cell_size;
			vec3 inv_cell_size;
			uint3 cell_count{};
			Cell outside_cell;
			std::vector<Cell> cells;
			std::vector<CellPair> cell_pairs;
			std::vector<Particle> particles;
		};
	}
}